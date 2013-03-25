/*
 * File:   libhuff.c
 * Author: Matthias Petri
 *
 * Canonical huffman encoding
 *
 */


#include "libutil.h"
#include "libhuff.h"
#include "libpqueue.h"

#define CHECK_BIT(var,pos) !!((var) & (1 << (pos)))

void
calc_code_len(hnode_t* n,uint32_t* clength,int32_t len)
{
    if (n != NULL) {
        if (n->sym != -1) { /* leaf */
            clength[n->sym] = len;
        } else {
            calc_code_len(n->left,clength,len+1);
            calc_code_len(n->right,clength,len+1);
        }
    }
}

void
calc_code_values(uint32_t* len,uint32_t* val)
{
    uint32_t min,max,l,i,j;
    int32_t limit[ALPHABET_SIZE] = {0};
    int32_t code_of_len[ALPHABET_SIZE] = {0};

    /* compute min max and code_of_len */
    min=255;
    max=0;
    for (i=0; i<ALPHABET_SIZE; i++) {
        j=len[i];
        if (j>0) {
            code_of_len[j]++;
            if (j>max) max = j;
            if (j<min) min = j;
        }
    }


    /* init limit: limit[i] is the largest code of a symbol of length i */
    limit[min]=code_of_len[min]-1;
    for (i=min+1; i<=max; i++) {
        limit[i] = ((limit[i-1]+1)<<1)+code_of_len[i]-1;
    }

    /* compute codes */
    for (i=0; i<ALPHABET_SIZE; i++) {
        l=len[i];
        val[i] = limit[l]-(--code_of_len[l]);
    }
}

void
encode_htree(uint32_t* clen,bit_file_t* of)
{
    uint32_t i,j;
    int32_t code_of_len[ALPHABET_SIZE] = {0};
    int32_t cstart[ALPHABET_SIZE] = {0};
    uint8_t syms[ALPHABET_SIZE] = {0};
    uint8_t n;


    /* sort by code length */
    n = 0;
    for (i=0; i<ALPHABET_SIZE; i++) {
        j=clen[i];
        code_of_len[j]++;
        if (j!=0) n++;
    }
    cstart[0] = 0;
    code_of_len[0] = 0;
    for (i=1; i<ALPHABET_SIZE; i++) {
        cstart[i] = cstart[i-1] + code_of_len[i-1];
    }

    for (i=0; i<ALPHABET_SIZE; i++) {
        if (clen[i]>0) {
            syms[ cstart[clen[i]] ] = i;
            cstart[clen[i]]++;
        }
    }

    /* output */
    BitFilePutChar(n-1,of); /* store n-1 */
    for (i=0; i<n; i++) BitFilePutChar(syms[i],of); /* output sym */
    for (i=0; i<n; i++) BitFilePutChar(clen[syms[i]],of); /* output len */

}


void
write_symbol(uint8_t sym,uint32_t* code_table,uint32_t* code_len,bit_file_t* of)
{
    int32_t i,n;
    uint32_t cw;

    n = code_len[sym];
    cw = code_table[sym];

    for (i=0; i<n; i++) {
        BitFilePutBit(CHECK_BIT(cw,n-i-1),of);
    }
}

void
encode_text(uint32_t* code_table,uint32_t* code_len,uint8_t* text,uint32_t n,bit_file_t* of)
{
    uint32_t i;

    /* write number of symbols */
    BitFilePutBitsInt(of,&n,32,sizeof(uint32_t));
    fprintf(stderr,"n = %d\n",n);

    /* encode the text */
    for (i=0; i<n; i++) {
        /* write symbol */
        write_symbol(text[i],code_table,code_len,of);
    }
}

/*
 * Main huffman encoding function. Use a heap based priority queue to create
 * the code words.
 */
void
encode_huffman(uint8_t* text,uint32_t n,bit_file_t* of)
{
    uint32_t i;
    pqueue_t* pq;
    hnode_t* left,*right,*new,*root;
    uint32_t code_len[ALPHABET_SIZE] = {0};
    uint32_t code_table[ALPHABET_SIZE] = {0};
    uint32_t freqs[ALPHABET_SIZE] = {0};

    /* count frequencies */
    for (i=0; i<n; i++) freqs[text[i]]++;

    /* create priority queue / heap */
    pq = pqueue_create();

    /* insert freqencies */
    for (i=0; i<ALPHABET_SIZE; i++) {
        if (freqs[i]>0) {
            new = (hnode_t*) safe_malloc(sizeof(hnode_t));
            new->sym = i;
            new->freq = freqs[i];
            new->left = NULL;
            new->right = NULL;

            pqueue_enqueue(pq,new,new->freq);
        }
    }

    /* create tree */
    root = NULL;
    while (pqueue_isempty(pq)==0) {
        left = pqueue_dequeue(pq);
        right = pqueue_dequeue(pq);

        if (left != NULL && right != NULL) {
            new = (hnode_t*) safe_malloc(sizeof(hnode_t));
            new->sym = -1;
            new->freq = left->freq + right->freq;
            new->left = (hnode_t*) left;
            new->right = (hnode_t*) right;
            pqueue_enqueue(pq,new,new->freq);
        } else {
            root = left;
        }
    }

    /* calculate code lengths + values */
    if (root) {
        calc_code_len(root,code_len,0);

        calc_code_values(code_len,code_table);
    }

    /* encode the tree info */
    encode_htree(code_len,of);

    /* encode the input */
    encode_text(code_table,code_len,text,n,of);
}
