/*
 * File:   libbwt.c
 * Author: Matthias Petri
 *
 *  BWT implementation based on
 *  "An efficient method for in memory construction of suffix arrays"
 *  by Itoh, Hideo and Tanaka, H.
 */
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   ds.c
   deep-shallow sorting algorithm. these routines are taken mainly by
   Seward's d_copyEQ_u12.c
   Buckets are sorted one at a time; when  bucket Y is completely sorted
   (except for YY) "pointer copying" is used to sort small buckets of
   the form XY, X=A..Z
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

#include "libutil.h"
#include "libbwt.h"
/* *******************************************************************
   globals.c
   Ver 1.0   14-oct-02
   This file contains the definition of the global variables
   which can be defined by the user + some relate procedures
   ******************************************************************* */
#include <assert.h>
#include <limits.h>

/* ---- global variables (modifiable by command line options) ----- */
int Anchor_dist;
int Shallow_limit;
int _ds_Verbose;
int _ds_Word_size;
int Mk_qs_thresh;
int32_t Max_pseudo_anchor_offset;

int32_t B2g_ratio;

int32_t Update_anchor_ranks;

int32_t Blind_sort_ratio;


int32_t  Text_size;
uint8_t*  Upper_text_limit;

int check_global_variables(void);
void set_global_variables(void);
int compute_overshoot(void);


uint8_t* Text;
int32_t* Sa;
int32_t ftab [65537];
uint16_t*  Anchor_offset;
int32_t*  Anchor_rank;

int32_t Calls_helped_sort=0;
int32_t Calls_anchor_sort_forw=0;
int32_t Calls_anchor_sort_backw=0;
int32_t Calls_pseudo_anchor_sort_forw=0;
int32_t Calls_deep_sort=0;



void deep_sort(int32_t* a, int32_t n, int32_t depth);




#define Get_small_bucket(pos) ((Text[pos]<<8) + Text[pos+1])



void general_anchor_sort(int32_t* a,int32_t n,int32_t pos,int32_t rank,int32_t off);

void pseudo_anchor_sort(int32_t* a, int32_t n,int32_t pseudo_an,int32_t offset);
int32_t split_group(int32_t* a, int n, int,int,int32_t,int*);
void update_anchors(int32_t* a, int32_t n);
void pseudo_or_deep_sort(int32_t* a, int32_t n, int32_t depth);










/* ------- node of blind trie -------- */
typedef struct nodex {
    int32_t skip;
    uint8_t key;
    struct nodex*  down;
    struct nodex* right;
} node;



#define BUFSIZE 1000
#define FREESIZE 5000
void* freearr[FREESIZE];
static node* bufn;
static int bufn_num=0, free_num=0;
static int32_t* Aux, Aux_written;
node** Stack;
int Stack_size;


int neg_integer_cmp(const void*, const void*);
node* find_companion(node* head, uint8_t* s);
void insert_suffix(node* h, int32_t suf, int n, uint8_t mmchar);
void traverse_trie(node* h);
int32_t compare_suffixes(int32_t suf1, int32_t suf2, int32_t depth);
void free_node_mem();
node* get_leaf(node* head);



/* ****************************************************************
   routine for deep-sorting the suffixes a[0] ... a[n-1]
   knowing that they have a common prefix of length "depth"
  **************************************************************** */
void blind_ssort(int32_t* a, int32_t n, int32_t depth)
{
    int32_t i,j,aj,lcp;
    node nh, *root, *h;


    qsort(a,n, sizeof(int32_t), neg_integer_cmp);


    for (j=0; j<n; j++)
        if (a[j]+depth < Text_size)
            break;
    if (j>=n-1) return;


    Stack = (node**) malloc(n*sizeof(node*));
    if (Stack==NULL) {
        fprintf(stderr,"Out of memory! (blind_ssort)\n");
        exit(1);
    }


    nh.skip = -1;   nh.right = NULL; nh.down = (void*) a[j];
    root = &nh;


    for (i=j+1; i<n; i++) {
        h=find_companion(root, Text+a[i]);
        assert(h->skip==-1);
        assert(Stack_size<=i-j);
        aj=(int32_t) h->down;
        assert(aj>a[i]);
        lcp = compare_suffixes(aj,a[i],depth);
        insert_suffix(root, a[i], lcp, Text[aj+lcp]);
    }


    Aux=a;  Aux_written = j;
    traverse_trie(root);
    assert(Aux_written==n);

    free_node_mem();
    free(Stack);
}

/* ***********************************************************************
   this function traverses the trie rooted at head following the string s.
   Returns the leaf "corresponding" to the string s
   *********************************************************************** */
node* find_companion(node* head, uint8_t* s)
{
    uint8_t c;
    node* p;
    int t;

    Stack_size = 0;
    while (head->skip >= 0) {
        Stack[Stack_size++] = head;
        t = head->skip;
        if (s+t>=Upper_text_limit)
            return get_leaf(head);
        c = s[t]; p = head->down;
repeat:
        if (c==p->key) {
            head = p;
            continue;
        } else if (c<p->key)
            return get_leaf(head);
        if ((p=(p->right))==NULL)
            return get_leaf(head);
        goto repeat;
    }
    Stack[Stack_size++] = head;
    return head;
}




node* get_leaf(node* head)
{
    assert(head->skip>=0);

    do {
        head = head->down;
    } while (head->skip>=0);
    return head;
}



node* new_node__blind_ssort(void)
{
    if (bufn_num-- == 0) {
        bufn = (node*) malloc(BUFSIZE * sizeof(node));
        if (bufn==NULL) {
            fprintf(stderr,"Out of mem (new_node1)\n"); exit(1);
        }
        freearr[free_num++] = (void*) bufn;
        if (free_num>=FREESIZE) {
            fprintf(stderr,"Out of mem (new_node2)\n"); exit(1);
        }
        bufn_num = BUFSIZE-1;
    }
    return bufn++;
}


/* *****************************************************
   insert a suffix in the trie rooted at *p.
   we know that the trie already contains a string
   which share the first n chars with suf
   ***************************************************** */
void insert_suffix(node* h, int32_t suf, int n, uint8_t mmchar)
{
    node* new_node__blind_ssort(void);
    int32_t t;
    uint8_t c, *s;
    node* p, **pp;

    s = Text + suf;

#if 0

    while ((t=h->skip) < n) {
        if (t < 0) break;
        c=s[t];  p=h->down;

repeat:
        if (c==p->key) {
            h = p;
            continue;
        }
        if (c>p->key && p->right!=NULL) {
            p=p->right;
            goto repeat;
        }

        fprintf(stderr,"Error in blind_sort (insert_string)\n");
        exit(1);
    }
#else
    for (t=0; t<Stack_size; t++) {
        h=Stack[t];
        if (h->skip<0 || h->skip>=n) break;
    }
#endif

    assert(s[n]!=mmchar || h->skip==-1 || h->skip==n);


    if (h->skip!=n) {
        p = new_node__blind_ssort();
        p->key = mmchar;
        p->skip = h->skip;
        p->down = h->down;
        p->right = NULL;
        h->skip = n;
        h->down = p;
    }
    assert(h->skip==n);


    c=s[n]; pp = &(h->down);
    while ((*pp)!=NULL) {
        if ((*pp)->key>=c)
            break;
        pp = &((*pp)->right);
    }

    p = new_node__blind_ssort();
    p->skip = -1;
    p->key = c;
    p->right = *pp; *pp = p;
    p->down = (void*) suf;
    return;
}

/* ************************************************************
   this procedures traverse the trie in depth first order
   so that the suffixes (stored in the leaf) are recovered
   in lexicographic order
   ************************************************************ */
void traverse_trie(node* h)
{
    node* p, *nextp;

    if (h->skip<0)
        Aux[Aux_written++] = (int32_t) h->down;
    else {
        p = h->down;
        assert(p!=NULL);
        do {
            nextp = p->right;
            if (nextp!=NULL) {
                assert(nextp->key>=p->key);


                if (nextp->key==p->key) {
                    traverse_trie(nextp);
                    traverse_trie(p);
                    p = nextp->right;
                    continue;
                }
            }
            traverse_trie(p);
            p=nextp;
        } while (p!=NULL);
    }
}



/* ***********************************************************************
   Function to compute the lcp of two strings originating from the *b1 and *b2
   the parameter is the length of s1 (which is shortest than s2)
   if s1 is a prefix of s2 we return the length of s1 -1
   The size of the unrolled loop must be at most equal to the costant
   Cmp_overshoot defined in common.h
   the function return the result of the comparison (+ or -) and writes
   in Cmp_done the number of comparisons done
   *********************************************************************** */
static
int32_t get_lcp_unrolled(uint8_t* b1, uint8_t* b2, int32_t cmp_limit)
{
    int32_t cmp2do;
    uint8_t c1, c2;
    assert(b1 != b2);



    cmp2do = cmp_limit;
    do {

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -=  1; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -=  2; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -=  3; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -=  4; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -=  5; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -=  6; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -=  7; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -=  8; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -=  9; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -= 10; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -= 11; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -= 12; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -= 13; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -= 14; break;
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            cmp2do -= 15; break;
        }
        b1++; b2++;

        cmp2do -= 16;
    } while (cmp2do>0);


    if (cmp_limit - cmp2do < cmp_limit)
        return cmp_limit-cmp2do;

    return cmp_limit-1;
}



/* ************************************************************************
   this function returns the lcp between suf1 and suf2 (that is returns n
   such that suf1[n]!=suf2[n] but suf1[i]==suf2[i] for i=0..n-1
   However, it is possible that suf1 is a prefix of suf2 (not vice-versa
   because of the initial sorting of suffixes in order of descreasing length)
   in this case the function returns n=length(suf1)-1. So in this case
   suf1[n]==suf2[n] (and suf1[n+1] does not exists).
   ************************************************************************ */
int32_t compare_suffixes(int32_t suf1, int32_t suf2, int32_t depth)
{
    int32_t get_lcp_unrolled(uint8_t*, uint8_t*, int32_t);
    int limit;
    uint8_t* s1, *s2;

    assert(suf1>suf2);
    s1  = Text + depth +suf1;
    s2  = Text + depth +suf2;
    limit = Text_size - suf1 - depth;
    return depth + get_lcp_unrolled(s1 ,s2, limit);
}



/* ******************************************************************
   comparison function used to sort suffixes in order of
   increasing length. Since suffixes are represented by their offset
   in the array, we sort these offsets in order of decreasing length.
   ****************************************************************** */
int neg_integer_cmp(const void* a, const void* b)
{
    return *((int32_t*) b) -  *((int32_t*) a);
}



void free_node_mem()
{
    int i;

    for (i=free_num-1; i>=0; i--) {
        assert(freearr[i]!=NULL);
        free(freearr[i]);
    }

    bufn_num=free_num=0;
}











/* ***********************************************************************
   Function to compare two strings originating from the *b1 and *b2
   The size of the unrolled loop must be at most equal to the costant
   Cmp_overshoot defined in common.h
   the function return the result of the comparison (+ or -) and writes
   in Cmp_done the number of successfull comparisons done
   *********************************************************************** */
static int32_t Cmp_done;

int32_t cmp_unrolled_lcp(uint8_t* b1, uint8_t* b2)
{

    uint8_t c1, c2;
    assert(b1 != b2);
    Cmp_done=0;



    do {

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done +=  1; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done +=  2; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done +=  3; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done +=  4; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done +=  5; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done +=  6; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done +=  7; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done +=  8; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done +=  9; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done += 10; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done += 11; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done += 12; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done += 13; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done += 14; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_done += 15; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        Cmp_done += 16;

    } while (b1<Upper_text_limit && b2<Upper_text_limit);


    return b2 - b1;
}

/* **************************************************************
   ternary quicksort (seward-like) with lcp information
   ************************************************************** */
#define STACK_SIZE 100
#define Swap(i,j) {tmp=a[i]; a[i]=a[j]; a[j]=tmp;}
#define Pushd(x,y,z) {stack_lo[sp]=x; stack_hi[sp]=y; stack_d[sp]=z; sp++;}
#define Popd(x,y,z)  {sp--; x=stack_lo[sp]; y=stack_hi[sp]; z=stack_d[sp];}
void qs_unrolled_lcp(int32_t* a, int n, int depth, int blind_limit)
{
    void blind_ssort(int32_t *a, int32_t n, int32_t depth);
    int32_t cmp_unrolled_lcp(uint8_t *b1, uint8_t *b2);
    uint8_t* text_depth, *text_pos_pivot;
    int32_t stack_lo[STACK_SIZE];
    int32_t stack_hi[STACK_SIZE];
    int32_t stack_d[STACK_SIZE];
    int32_t sp,r,r3,med,tmp;
    int32_t i, j, lo, hi,ris,lcp_lo,lcp_hi;


    r=sp=0;
    Pushd(0,n-1,depth);


    while (sp > 0) {
        assert(sp < STACK_SIZE);
        Popd(lo,hi,depth);
        text_depth = Text+depth;


        if (hi-lo<blind_limit) {
            blind_ssort(a+lo,hi-lo+1,depth);
            continue;
        }

        /* Random partitioning. Guidance for the magic constants
           7621 and 32768 is taken from Sedgewick's algorithms
           book, chapter 35.
        */
        r = ((r * 7621) + 1) % 32768;
        r3 = r % 3;
        if (r3 == 0) med = lo; else if (r3 == 1) med = (lo+hi)>>1; else
            med = hi;


        Swap(med,hi);
        text_pos_pivot=text_depth+a[hi];
        i=lo-1; j=hi;
        lcp_lo=lcp_hi=INT_MAX;
        while (1) {
            while (++i<hi) {
                ris=cmp_unrolled_lcp(text_depth+a[i], text_pos_pivot);
                if (ris>0) {
                    if (Cmp_done < lcp_hi) lcp_hi=Cmp_done; break;
                } else if (Cmp_done < lcp_lo) lcp_lo=Cmp_done;
            }
            while (--j>lo) {
                ris=cmp_unrolled_lcp(text_depth+a[j], text_pos_pivot);
                if (ris<0) {
                    if (Cmp_done < lcp_lo) lcp_lo=Cmp_done;
                    break;
                } else if (Cmp_done < lcp_hi) lcp_hi=Cmp_done;
            }
            if (i >= j) break;
            Swap(i,j);
        }
        Swap(i,hi);


        assert(lcp_lo<INT_MAX || i==lo);
        assert(lcp_hi<INT_MAX || i==hi);


        if (i-lo < hi-i) {
            Pushd(i+1,hi,depth+lcp_hi);
            if (i-lo>1) Pushd(lo,i-1,depth+lcp_lo);
        } else {
            Pushd(lo,i-1,depth+lcp_lo);
            if (hi-i>1) Pushd(i+1,hi,depth+lcp_hi);
        }
    }
}



/* ****************************************************************
   routine for deep-sorting the suffixes a[0] ... a[n-1]
   knowing that they have a common prefix of length "depth"
  **************************************************************** */
void deep_sort(int32_t* a, int32_t n, int32_t depth)
{
    void blind_ssort(int32_t *a, int32_t n, int32_t depth);
    int blind_limit;

    Calls_deep_sort++;
    assert(n>1);

    blind_limit=Text_size/Blind_sort_ratio;
    if (n<=blind_limit)
        blind_ssort(a,n,depth);
    else
        qs_unrolled_lcp(a,n,depth,blind_limit);
}



/* *****************************************************************
   This procedure sort the strings a[0] ... a[n-1] with the help of an
   anchor. The real sorting is done by the procedure
   anchor_sort(). Here we choose the anchor.  The parameter depth is
   the number of chars that a[0] ... a[n-1] are known to have in
   common (thus a direct comparison among a[i] and a[j] should start
   from position depth) Note that a[] is a subsection of the sa therefore
   a[0] ... a[n-1] are starting position of suffixes
   For every a[i] we look at the anchor a[i]/Anchor_dist and the one
   after that. This justifies the definition of Anchor_num (the size of
   Anchor_ofset[] and Anchor_rank[] defined in ds_sort()) as
     Anchor_num = 2 + (n-1)/Anchor_dist
   ***************************************************************** */
void helped_sort(int32_t* a, int n, int depth)
{
    int32_t i, curr_sb, diff, toffset, aoffset;
    int32_t text_pos, anchor_pos, anchor, anchor_rank;
    int32_t min_forw_offset, min_forw_offset_buc, max_back_offset;
    int32_t best_forw_anchor, best_forw_anchor_buc, best_back_anchor;
    int32_t forw_anchor_index, forw_anchor_index_buc, back_anchor_index;

    Calls_helped_sort++;
    if (n==1) goto done_sorting;


    if (Anchor_dist==0) {
        pseudo_or_deep_sort(a, n, depth);
        return;
    }


    curr_sb = Get_small_bucket(a[0]);


    min_forw_offset = min_forw_offset_buc = INT_MAX;
    max_back_offset = INT_MIN;
    best_forw_anchor = best_forw_anchor_buc = best_back_anchor = -1;
    forw_anchor_index = forw_anchor_index_buc = back_anchor_index = -1;

    for (i=0; i<n; i++) {
        text_pos = a[i];

        anchor = text_pos/Anchor_dist;
        toffset = text_pos % Anchor_dist;
        aoffset = Anchor_offset[anchor];
        if (aoffset<Anchor_dist) {
            diff = aoffset - toffset;
            assert(diff!=0);
            if (diff>0) {
                if (curr_sb!=Get_small_bucket(text_pos+diff)) {
                    if (diff<min_forw_offset) {
                        min_forw_offset = diff;
                        best_forw_anchor = anchor;
                        forw_anchor_index = i;
                    }
                } else {
                    if (diff<min_forw_offset_buc) {
                        min_forw_offset_buc = diff;
                        best_forw_anchor_buc = anchor;
                        forw_anchor_index_buc = i;
                    }
                }
            } else {
                if (diff>max_back_offset) {
                    max_back_offset = diff;
                    best_back_anchor = anchor;
                    back_anchor_index = i;
                }

                aoffset = Anchor_offset[++anchor];
                if (aoffset<Anchor_dist) {
                    diff = Anchor_dist + aoffset - toffset;
                    assert(diff>0);
                    if (curr_sb!=Get_small_bucket(text_pos+diff)) {
                        if (diff<min_forw_offset) {
                            min_forw_offset = diff;
                            best_forw_anchor = anchor;
                            forw_anchor_index = i;
                        }
                    } else {
                        if (diff<min_forw_offset_buc) {
                            min_forw_offset_buc = diff;
                            best_forw_anchor_buc = anchor;
                            forw_anchor_index_buc = i;
                        }
                    }
                }
            }
        }
    }

    if (best_forw_anchor>=0 && min_forw_offset<depth-1) {
        Calls_anchor_sort_forw++;
        assert(min_forw_offset<2*Anchor_dist);
        anchor_pos = a[forw_anchor_index] + min_forw_offset;
        anchor_rank = Anchor_rank[best_forw_anchor];
        assert(Sa[anchor_rank]==anchor_pos);
        general_anchor_sort(a,n,anchor_pos,anchor_rank,min_forw_offset);
        goto done_sorting;
    }

    if (best_back_anchor>=0) {
        uint8_t* T0, *Ti; int j;

        assert(max_back_offset>-Anchor_dist && max_back_offset<0);

        for (i=0; i<n; i++) {
            if (a[i]+max_back_offset<0)
                goto fail;
        }

        T0 = Text + a[0];
        for (i=1; i<n; i++) {
            Ti = Text + a[i];
            for (j=max_back_offset; j<= -1; j++)
                if (T0[j]!=Ti[j]) goto fail;
        }

        Calls_anchor_sort_backw++;
        anchor_pos = a[back_anchor_index] + max_back_offset;
        anchor_rank = Anchor_rank[best_back_anchor];
        assert(Sa[anchor_rank]==anchor_pos);
        general_anchor_sort(a,n,anchor_pos,anchor_rank,max_back_offset);
        goto done_sorting;
    }
fail:

    if (best_forw_anchor_buc>=0 && min_forw_offset_buc<depth-1) {
        int equal,lower,upper;

        assert(min_forw_offset_buc<2*Anchor_dist);
        anchor_pos = a[forw_anchor_index_buc] + min_forw_offset_buc;
        anchor_rank = Anchor_rank[best_forw_anchor_buc];
        assert(Sa[anchor_rank]==anchor_pos);


        equal=split_group(a,n,depth,min_forw_offset_buc,
                          forw_anchor_index_buc,&lower);
        if (equal==n) {
            Calls_anchor_sort_forw++;
            general_anchor_sort(a,n,anchor_pos,anchor_rank,min_forw_offset_buc);
        } else {

            upper = n-equal-lower;
            assert(upper>=0);


            Calls_anchor_sort_forw++;
            if (equal>1)
                general_anchor_sort(a+lower,equal,anchor_pos,anchor_rank,
                                    min_forw_offset_buc);


            if (lower>1) pseudo_or_deep_sort(a,lower,depth);
            if (upper>1) pseudo_or_deep_sort(a+lower+equal,upper,depth);
        }
        goto done_sorting;
    }





    pseudo_or_deep_sort(a, n, depth);
done_sorting:

    if (Anchor_dist>0) update_anchors(a, n);
}



/* *******************************************************************
   try pseudo_anchor sort or deep_sort
   ******************************************************************** */
void pseudo_or_deep_sort(int32_t* a, int32_t n, int32_t depth)
{

    int32_t offset, text_pos, sb, pseudo_anchor_pos, max_offset, size;


    if (Max_pseudo_anchor_offset>0) {

        max_offset = MIN(depth-1,Max_pseudo_anchor_offset);
        text_pos = a[0];
        for (offset=1; offset<max_offset; offset++) {
            pseudo_anchor_pos = text_pos+offset;
            sb = Get_small_bucket(pseudo_anchor_pos);

            if (IS_SORTED_BUCKET(sb)) {
                size=BUCKET_SIZE(sb);
                if (size>B2g_ratio*n) continue;

                pseudo_anchor_sort(a,n,pseudo_anchor_pos,offset);
                Calls_pseudo_anchor_sort_forw++;
                return;
            }
        }
    }
    deep_sort(a,n,depth);
}

/* ********************************************************************
   this routine sorts the suffixes a[0] ... a[n-1] using the fact that
   in their common prefix, after offset characters, there is a
   suffix which is in an already sorted bucket. This suffix is called
   a pseudo anchor since it is used essentially as an anchor, but
   it is not in an anchor position (=position multiple of Anchor_dist)
   ******************************************************************** */
void pseudo_anchor_sort(int32_t* a,int32_t n,int32_t pseudo_anchor_pos, int32_t offset)
{
    int32_t get_rank(int32_t);
    int32_t get_rank_update_anchors(int32_t);
    int32_t pseudo_anchor_rank;


    if (Update_anchor_ranks!=0 && Anchor_dist>0)
        pseudo_anchor_rank = get_rank_update_anchors(pseudo_anchor_pos);
    else
        pseudo_anchor_rank = get_rank(pseudo_anchor_pos);

    assert(Sa[pseudo_anchor_rank]==pseudo_anchor_pos);

    general_anchor_sort(a,n,pseudo_anchor_pos,pseudo_anchor_rank,offset);
}


/* ********************************************************
   macros for marking integers: works assuming integers have
   at least 32 bit and that the 32nd bit is not used
   This simply means that the text size can be at most 2GB
   ********************************************************* */
#define MARKER (1<<31)
#define MARK(i) {                \
  assert(( Sa[i]&MARKER) == 0);  \
  (Sa[i] |= MARKER);             \
}
#define ISMARKED(i) (Sa[i] & MARKER)
#define UNMARK(i) (Sa[i] &= ~MARKER)

/* ********************************************************************
   This routines sorts a[0] ... a[n-1] using the fact that
   in their common prefix, after offset characters, there is a
   suffix whose rank is known. In this routine we call this suffix anchor
   (and we denote its position and rank with anchor_pos and anchor_rank
   respectively) but it is not necessarily an anchor (=does not necessarily
   starts at position multiple of Anchor_dist) since this function is
   called by pseudo_anchor_sort().
   The routine works by scanning the suffixes before and after the anchor
   in order to find (and mark) those which are suffixes of a[0] ... a[n-1].
   After that, the ordering of a[0] ... a[n-1] is derived with a sigle
   scan of the marked suffixes.
   ******************************************************************** */
void general_anchor_sort(int32_t* a, int32_t n,
                         int32_t anchor_pos, int32_t anchor_rank, int32_t offset)
{
    int integer_cmp(const void*, const void*);
    int32_t sb, lo, hi;
    int32_t curr_lo, curr_hi, to_be_found, i,j;
    int32_t item;
    void* ris;

    assert(Sa[anchor_rank]==anchor_pos);
    /* ---------- get bucket of anchor ---------- */
    sb = Get_small_bucket(anchor_pos);
    lo = BUCKET_FIRST(sb);
    hi = BUCKET_LAST(sb);
    assert(sb==Get_small_bucket(a[0]+offset));

    qsort(a,n, sizeof(int32_t), integer_cmp);






    curr_hi = curr_lo = anchor_rank;


#if DEBUG
    item = anchor_pos-offset;
    assert(bsearch(&item,a,n,sizeof(int32_t), integer_cmp));
#endif

    MARK(curr_lo);

    for (to_be_found=n-1; to_be_found>0;) {

        assert(curr_lo > lo || curr_hi < hi);
        while (curr_lo > lo) {
            item = Sa[--curr_lo]-offset;
            ris = bsearch(&item,a,n,sizeof(int32_t), integer_cmp);
            if (ris)	{
                MARK(curr_lo);
                to_be_found--;
            } else	break;
        }
        while (curr_hi < hi) {
            item = Sa[++curr_hi]-offset;
            ris = bsearch(&item,a,n,sizeof(int32_t), integer_cmp);
            if (ris)	{
                MARK(curr_hi);
                to_be_found--;
            } else      break;
        }
    }

    for (j=0, i=curr_lo; i<=curr_hi; i++)
        if (ISMARKED(i)) {
            UNMARK(i);
            a[j++] = Sa[i] - offset;
        }
    assert(j==n);
}

/* ********************************************************************
   compute the rank of the suffix starting at pos.
   It is required that the suffix is in an already sorted bucket
   ******************************************************************** */
int32_t get_rank(int32_t pos)
{
    int32_t sb, lo, hi, j;

    sb = Get_small_bucket(pos);
    if (!IS_SORTED_BUCKET(sb)) {
        fprintf(stderr,"Illegal call to get_rank! (get_rank1)\n");
        exit(1);
    }
    lo = BUCKET_FIRST(sb);
    hi = BUCKET_LAST(sb);
    for (j=lo; j<=hi; j++)
        if (Sa[j]==pos) return j;
    fprintf(stderr,"Illegal call to get_rank! (get_rank2)\n");
    exit(1);
    return 1;
}

/* ********************************************************************
   compute the rank of the suffix starting at pos. At the same time
   check if the rank of the suffixes in the bucket containing pos
   can be used to update some entries in Anchor_offset[] and Anchor_rank[]
   It is required that the suffix is in an already sorted bucket
   ******************************************************************** */
uint8_t bucket_ranked[65536];
int32_t get_rank_update_anchors(int32_t pos)
{
    int32_t get_rank(int32_t pos);
    int32_t sb, lo, hi, j, toffset, aoffset, anchor, rank;

    assert(Anchor_dist>0);

    sb = Get_small_bucket(pos);
    if (!(IS_SORTED_BUCKET(sb))) {
        fprintf(stderr,"Illegal call to get_rank! (get_rank_update_anchors)\n");
        exit(1);
    }

    if (bucket_ranked[sb]) return get_rank(pos);

    bucket_ranked[sb]=1;
    rank = -1;
    lo = BUCKET_FIRST(sb);
    hi = BUCKET_LAST(sb);
    for (j=lo; j<=hi; j++) {

        toffset = Sa[j]%Anchor_dist;
        anchor  = Sa[j]/Anchor_dist;
        aoffset = Anchor_offset[anchor];
        if (toffset<aoffset) {
            Anchor_offset[anchor] = toffset;
            Anchor_rank[anchor] = j;
        }

        if (Sa[j]==pos) {
            assert(rank==-1); rank=j;
        }
    }
    assert(rank>=0);
    return rank;
}


/* ******************************************************************
   comparison function used to sort pointers as if they were integers
   if a pointer does not correspond to an int32_t we must change the
   function accordingly
   ****************************************************************** */
int integer_cmp(const void* a, const void* b)
{
    return *((int32_t*) a) -  *((int32_t*) b);
}

/* ****************************************************************
   given a SORTED array of suffixes a[0] .. a[n-1]
   updates Anchor_rank[] and Anchor_offset[]
   **************************************************************** */
void update_anchors(int32_t* a, int32_t n)
{
    int32_t i,anchor,toffset,aoffset,text_pos;

    assert(Anchor_dist>0);
    for (i=0; i<n; i++) {
        text_pos = a[i];

        anchor = text_pos/Anchor_dist;
        toffset = text_pos % Anchor_dist;
        aoffset = Anchor_offset[anchor];
        if (toffset<aoffset) {
            Anchor_offset[anchor] = toffset;
            Anchor_rank[anchor] = (a - Sa) + i;
            assert(Sa[Anchor_rank[anchor]]==
                   anchor*Anchor_dist+Anchor_offset[anchor]);
        }
    }
}



/* *******************************************************************
   This function takes as input an array a[0] .. a[n-1] of suffixes
   which share the first "depth" chars. "pivot" in an index in 0..n-1
   and offset and integer>0. The function splits a[0] .. a[n-1]
   into 3 groups: first the suffixes which are smaller than a[pivot],
   then those which are equal to a[pivot] and finally those which are
   greater than a[pivot]. Here, smaller, equal, larger refer to
   a lexicographic ordering limited to the first depth+offest chars
   (since the first depth chars are equal we only look at the chars
   in position depth, depth+1, ... depth+offset-1).
   The function returns the number "num" of suffixes equal to a[pivot],
   and stores in *first the first of these suffixes. So at the end
   the smaller suffixes are in a[0] ... a[first-1],
   the equal suffixes in a[first] ... a[first+num-1],
   the larger suffixes in a[first+num] ... a[n-1]
   The splitting is done using a modified mkq()
   ******************************************************************* */
#define swap2(a, b) { t = *(a); *(a) = *(b); *(b) = t; }
#define ptr2char(i) (*(*(i) + text_depth))

int32_t split_group(int32_t* a, int n, int depth,int offset,int32_t pivot,int* first)
{
    void vecswap2(int32_t *a, int32_t *b, int n);
    int r, partval;
    int32_t* pa, *pb, *pc, *pd, *pa_old, *pd_old, pivot_pos, t;
    uint8_t* text_depth,*text_limit;


    pivot_pos = a[pivot];
    text_depth = Text+depth;
    text_limit = text_depth+offset;



    pa = a; pd = a + n-1;

    for (; pa!=pd && (text_depth<text_limit); text_depth++) {
        assert(pa<pd);


        partval = text_depth[pivot_pos];

        pb = pa_old = pa;
        pc = pd_old = pd;
        for (;;) {
            while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) {
                if (r == 0) {
                    swap2(pa, pb);
                    pa++;
                }
                pb++;
            }
            while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) {
                if (r == 0) {
                    swap2(pc, pd);
                    pd--;
                }
                pc--;
            }
            if (pb > pc) break;
            swap2(pb, pc);
            pb++;
            pc--;
        }
        r = MIN(pa-pa_old, pb-pa); vecswap2(pa_old,  pb-r, r);
        r = MIN(pd-pc, pd_old-pd); vecswap2(pb, pd_old+1-r, r);

        pa = pa_old + (pb-pa);
        pd = pd_old - (pd-pc);

    }
    *first=pa-a;
    assert(pd-pa>=0);
    return pd-pa+1;
}






/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   shallow.c
   This is the multikey quicksort from bentley-sedgewick modified
   so that it stops recursion when depth reaches  Shallow_limit
   (that is when two or more suffixes have Shallow_limit chars in common).
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


extern uint8_t* Text;
extern uint8_t* Upper_text_limit;
extern int32_t _ds_Word_size;
extern int32_t Mk_qs_thresh;




int32_t Shallow_limit;
uint8_t* Shallow_text_limit;

#define UNROLL 1





void shallow_sort(int32_t* a, int n, int shallow_limit)
{

    Shallow_limit = shallow_limit;
    Shallow_text_limit = Text + shallow_limit;


    switch (_ds_Word_size) {
        case(1): shallow_mkq(a, n, Text+2); break;
        case(2): shallow_mkq16(a, n, Text+2); break;
        case(4): shallow_mkq32(a, n, Text+2); break;
        default:
            fprintf(stderr,
                    "Invalid word size for mkqs (%d) (shallow_sort)\n",_ds_Word_size);
            exit(1);
    }
}


/* =======================================================
   auxiliary procedures and macro for bentley-sedgewick's
   multikey quicksort
   ======================================================= */
void vecswap2(int32_t* a, int32_t* b, int n)
{
    while (n-- > 0) {
        int32_t t = *a;
        *a++ = *b;
        *b++ = t;
    }
}

#define swap2(a, b) { t = *(a); *(a) = *(b); *(b) = t; }
#define ptr2char(i) (*(*(i) + text_depth))

int32_t* med3func(int32_t* a, int32_t* b, int32_t* c, uint8_t* text_depth)
{
    int va, vb, vc;
    if ((va=ptr2char(a)) == (vb=ptr2char(b)))
        return a;
    if ((vc=ptr2char(c)) == va || vc == vb)
        return c;
    return va < vb ?
           (vb < vc ? b : (va < vc ? c : a))
               : (vb > vc ? b : (va < vc ? a : c));
}
#define med3(a, b, c) med3func(a, b, c, text_depth)


/* ********************************************************
   recursive multikey quicksort from Bentley-Sedgewick
   stops when text_depth reaches Shallow_depth_limit
   that is when we have found that the current set of strings
   have Shallow_limit chars in common
   ******************************************************** */
void shallow_mkq(int32_t* a, int n, uint8_t* text_depth)
{
    void vecswap2(int32_t *a, int32_t *b, int n);
    int d, r, partval;
    int32_t* pa, *pb, *pc, *pd, *pl, *pm, *pn, t;
    uint8_t* next_depth;


    if (n < Mk_qs_thresh) {
        shallow_inssort_lcp(a, n, text_depth);
        return;
    }


repeat:
    pl = a;
    pm = a + (n/2);
    pn = a + (n-1);
    if (n > 30) {
        d = (n/8);
        pl = med3(pl, pl+d, pl+2*d);
        pm = med3(pm-d, pm, pm+d);
        pn = med3(pn-2*d, pn-d, pn);
    }
    pm = med3(pl, pm, pn);
    swap2(a, pm);
    partval = ptr2char(a);
    pa = pb = a + 1;
    pc = pd = a + n-1;

    for (;;) {
        while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) {
            if (r == 0) {
                swap2(pa, pb);
                pa++;
            }
            pb++;
        }
        while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) {
            if (r == 0) {
                swap2(pc, pd);
                pd--;
            }
            pc--;
        }
        if (pb > pc) break;
        swap2(pb, pc);
        pb++;
        pc--;
    }

#if UNROLL
    if (pa>pd) {

        if ((next_depth = text_depth+1) >= Shallow_text_limit) {
            helped_sort(a, n, next_depth-Text);
            return;
        } else {
            text_depth = next_depth;
            goto repeat;
        }
    }
#endif

    pn = a + n;
    r = MIN(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
    r = MIN(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);

    if ((r = pb-pa) > 1)
        shallow_mkq(a, r, text_depth);

    if ((next_depth = text_depth+1) < Shallow_text_limit)
        shallow_mkq(a + r, pa-pd+n-1, next_depth);
    else
        helped_sort(a + r, pa-pd+n-1, next_depth-Text);
    if ((r = pd-pc) > 1)
        shallow_mkq(a + n-r, r, text_depth);
}



/* ************** 16 *************** */
#define ptr2char16(i) (getword16(*(i) + text_depth))
#define getword16(s) ((unsigned)((*(s) << 8) | *((s)+1)))

#if 0
__inline__ int32_t* med3func16(int32_t* a, int32_t* b, int32_t* c, uint8_t* text_depth)
{
    int va, vb, vc;
    if ((va=ptr2char16(a)) == (vb=ptr2char16(b)))
        return a;
    if ((vc=ptr2char16(c)) == va || vc == vb)
        return c;
    return va < vb ?
           (vb < vc ? b : (va < vc ? c : a))
               : (vb > vc ? b : (va < vc ? a : c));
}
#define med3_16(a, b, c) med3func16(a, b, c, text_depth)
#endif

void shallow_mkq16(int32_t* a, int n, uint8_t* text_depth)
{
    void vecswap2(int32_t *a, int32_t *b, int n);
    int d, r, partval;
    int32_t* pa, *pb, *pc, *pd, *pl, *pm, *pn, t;
    uint8_t* next_depth;


    if (n < Mk_qs_thresh) {
        shallow_inssort_lcp(a, n, text_depth);
        return;
    }


repeat:
    pl = a;
    pm = a + (n/2);
    pn = a + (n-1);
    if (n > 30) {
        d = (n/8);
        pl = med3(pl, pl+d, pl+2*d);
        pm = med3(pm-d, pm, pm+d);
        pn = med3(pn-2*d, pn-d, pn);
    }
    pm = med3(pl, pm, pn);
    swap2(a, pm);
    partval = ptr2char16(a);
    pa = pb = a + 1;
    pc = pd = a + n-1;

    for (;;) {
        while (pb <= pc && (r = ptr2char16(pb)-partval) <= 0) {
            if (r == 0) {
                swap2(pa, pb);
                pa++;
            }
            pb++;
        }
        while (pb <= pc && (r = ptr2char16(pc)-partval) >= 0) {
            if (r == 0) {
                swap2(pc, pd);
                pd--;
            }
            pc--;
        }
        if (pb > pc) break;
        swap2(pb, pc);
        pb++;
        pc--;
    }
#if UNROLL
    if (pa>pd) {

        if ((next_depth = text_depth+2) >= Shallow_text_limit) {
            helped_sort(a, n, next_depth-Text);
            return;
        } else {
            text_depth = next_depth;
            goto repeat;
        }
    }
#endif

    pn = a + n;
    r = MIN(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
    r = MIN(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);

    if ((r = pb-pa) > 1)
        shallow_mkq16(a, r, text_depth);

    if ((next_depth = text_depth+2) < Shallow_text_limit)
        shallow_mkq16(a + r, pa-pd+n-1, next_depth);
    else
        helped_sort(a + r, pa-pd+n-1, next_depth-Text);
    if ((r = pd-pc) > 1)
        shallow_mkq16(a + n-r, r, text_depth);
}


/* *************** 32 **************** */
#define ptr2char32(i) (getword32(*(i) + text_depth))
#define getword32(s) ((unsigned)( (*(s) << 24) | ((*((s)+1)) << 16) \
                                  | ((*((s)+2)) << 8) | (*((s)+3)) ))
void shallow_mkq32(int32_t* a, int n, uint8_t* text_depth)
{
    void vecswap2(int32_t *a, int32_t *b, int n);
    uint32_t partval, val;
    int32_t* pa, *pb, *pc, *pd, *pl, *pm, *pn, t, d, r;
    uint8_t* next_depth;


    if (n < Mk_qs_thresh) {
        shallow_inssort_lcp(a, n, text_depth);
        return;
    }


repeat:
    pl = a;
    pm = a + (n/2);
    pn = a + (n-1);
    if (n > 30) {
        d = (n/8);
        pl = med3(pl, pl+d, pl+2*d);
        pm = med3(pm-d, pm, pm+d);
        pn = med3(pn-2*d, pn-d, pn);
    }
    pm = med3(pl, pm, pn);
    swap2(a, pm);
    partval = ptr2char32(a);
    pa = pb = a + 1;
    pc = pd = a + n-1;

    for (;;) {
        while (pb <= pc && (val=ptr2char32(pb)) <=  partval) {
            if (val == partval) {
                swap2(pa, pb);
                pa++;
            }
            pb++;
        }
        while (pb <= pc && (val=ptr2char32(pc)) >= partval) {
            if (val == partval) {
                swap2(pc, pd);
                pd--;
            }
            pc--;
        }
        if (pb > pc) break;
        swap2(pb, pc);
        pb++;
        pc--;
    }
#if UNROLL
    if (pa>pd) {

        if ((next_depth = text_depth+4) >= Shallow_text_limit) {
            helped_sort(a, n, next_depth-Text);
            return;
        } else {
            text_depth = next_depth;
            goto repeat;
        }
    }
#endif

    pn = a + n;
    r = MIN(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
    r = MIN(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);

    if ((r = pb-pa) > 1)
        shallow_mkq32(a, r, text_depth);

    if ((next_depth = text_depth+4) < Shallow_text_limit)
        shallow_mkq32(a + r, pa-pd+n-1, next_depth);
    else
        helped_sort(a + r, pa-pd+n-1, next_depth-Text);
    if ((r = pd-pc) > 1)
        shallow_mkq32(a + n-r, r, text_depth);
}



/* >>>>>>>>>>>>>>>>>>>>>> insertion sort routines >>>>>>>>>>>>>>>>>>>
   This insertion sort routines sorts the suffixes a[0] .. a[n-1]
   which have a common prexif of length text_depth-Text.
   The comparisons are done going at most at depth Shallow_limit;
   suffixes which have Shallow_limit chars in common are sorted using
   helped_sort().
   This inserion_sort keeps trak of the lcp in order to speed up
   the sorting.
  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

/* ***********************************************************************
   Function to compare two strings originating from the *b1 and *b2
   The size of the unrolled loop must be at most equal to the costant
   Cmp_overshoot defined in common.h
   When the function is called Cmp_left must contain the maximum number of
   comparisons the algorithm can do before returning 0 (equal strings)
   At exit Cmp_left has been decreased by the # of comparisons done
   *********************************************************************** */
static int32_t Cmp_left;

int32_t cmp_unrolled_shallow_lcp(uint8_t* b1, uint8_t* b2)
{

    uint8_t c1, c2;
    assert(b1 != b2);



    do {

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -=  1; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -=  2; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -=  3; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -=  4; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -=  5; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -=  6; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -=  7; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -=  8; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -=  9; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -= 10; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -= 11; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -= 12; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -= 13; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -= 14; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        c1 = *b1; c2 = *b2;
        if (c1 != c2) {
            Cmp_left -= 15; return ((uint32_t)c1 - (uint32_t)c2);
        }
        b1++; b2++;

        Cmp_left -= 16;
        if (Cmp_left<=0) return 0;

    } while (1);

    return b2 - b1;
}

/* *****************************************************************
   this is the insertion sort routine called by multikey-quicksort
   for sorting small groups.
   During insertion sort the comparisons are done calling
   cmp_unrolled_shallow_lcp() and two strings are equal if the coincides
   for Shallow_limit characters.
   After this first phase we sort groups of "equal_string" using
   helped_sort().
   Usage of lcp.
   For i=1,...n-1 let lcp[i] denote the lcp between a[i] and a[i+1].
   assume a[0] ... a[j-1] are already ordered and that we want to
   insert a new element ai. If suf(ai) >= suf(a[j-1]) we are done.
   If suf(ai)<suf(a[j-1]) we notice that: if lcpi>lcp[j-2] then
   suf(ai)>suf(a[j-2]) and we can stop since
   j-2 mmmmmmg
   j-1 mmmmmmmmmmmm
   ai  mmmmmmmmmmmf ] lcpi
   so we write a[j-1] in position j and ai in position j-1.
   if lcpi==lcp[j-2] then we need to compare suf(ai) with suf(a[j-2])
   j-2 mmmmmmmmmmm?     we can have either ?<f or ?>f or ?==f
   j-1 mmmmmmmmmmmm
   j   mmmmmmmmmmmf
   so we move a[j-1] to position j and compare suf(ai) with suf(a[j-2])
   starting from lcpi.
   Finally, if lcpi<lcp[j-2] then
   j-2 mmmmmmmmmmmmmmmmmmg
   j-1 mmmmmmmmmmmmmmmmmmm
   j   mmmmmmmmmmmmmf
   hence we have suf(ai)<suf(a[j-2]) and we consider a[j-3];
   if lcpi<lcp[j-3] we go on look at a[j-4] and go on.
   if lcp[j]>lcp[j-3] we are in the following position:
   j-3 mmmmmmc
   j-2 mmmmmmmmmmmmmmmg
   j-1 mmmmmmmmmmmmmmmm
   j   mmmmmmmmmmf
   and we know that suf(ai) is larger than suf(a[j-3]). If we find that
   lcpi==lcp[j-3] then we must compare suf(ai) with suf(a[j-3])
   but starting with position lcpi
   ***************************************************************** */
static int lcp_aux[1+Max_thresh];
static int* lcp=lcp_aux+1;
void shallow_inssort_lcp(int32_t* a, int32_t n, uint8_t* text_depth)
{
    int32_t cmp_unrolled_shallow_lcp(uint8_t*, uint8_t*);
    int32_t i, j, j1, lcp_new, r, ai,lcpi;
    int32_t cmp_from_limit;
    uint8_t* text_depth_ai;


    lcp_aux[0] = -1;
    for (i=0; i<n; i++) lcp[i]=0;

    cmp_from_limit = Shallow_text_limit-text_depth;


    for (i = 1; i< n ; i++) {
        ai = a[i]; lcpi = 0;
        text_depth_ai = ai + text_depth;
        j=i; j1=j-1;
        while (1) {


            Cmp_left = cmp_from_limit-lcpi;
            r = cmp_unrolled_shallow_lcp(lcpi+a[j1]+text_depth,lcpi+text_depth_ai);
            lcp_new = cmp_from_limit - Cmp_left;
            assert(r!=0 || lcp_new>= cmp_from_limit);

            if (r<=0) {
                lcp[j1]=lcp_new;
                break;
            }



            lcpi = lcp_new;
            do {
                a[j] = a[j1];
                lcp[j] = lcp[j1];
                j=j1; j1--;
            } while (lcpi<lcp[j1]);

            if (lcpi>lcp[j1]) break;



        }
        a[j]=ai;
        lcp[j]=lcpi;
    }


    for (i=0; i<n-1; i=j+1) {
        for (j=i; j<n ; j++)
            if (lcp[j]<cmp_from_limit) break;
        if (j-i>0)
            helped_sort(a+i,j-i+1,Shallow_limit);
    }
}



















/* *******************************************************************
   procedure to be called by external program before calling ds_ssort()
   using this procedure external programs can choose
   the parameters Anchor_dist and Blind_sort_ratio.
   The procedure returns 0 if something goes wrong, otherwise
   it returns the overshhot, that is the amount of extra space
   required at the end of the array contanining the text
   ******************************************************************** */
int init_ds_ssort(int adist, int bs_ratio)
{
    set_global_variables();
    Anchor_dist = adist;
    Blind_sort_ratio=bs_ratio;
    Shallow_limit =  Anchor_dist + 50;
    if (check_global_variables())
        return 0;
    return compute_overshoot();
}



void set_global_variables(void)
{
    Blind_sort_ratio=2000;
    Anchor_dist = 500;
    Shallow_limit = 550;
    _ds_Verbose = 0;
    _ds_Word_size = 4;
    Mk_qs_thresh=20;
    Max_pseudo_anchor_offset=0;
    B2g_ratio=1000;
    Update_anchor_ranks=0;
}



int check_global_variables(void)
{
    if ((Anchor_dist<100) && (Anchor_dist!=0)) {
        fprintf(stderr,"Anchor distance must be 0 or greater than 99\n");
        return 1;
    }
    if (Anchor_dist>65535) {
        fprintf(stderr,"Anchor distance must be less than 65536\n");
        return 1;
    }
    if (Shallow_limit<2) {
        fprintf(stderr,"Illegal limit for shallow sort\n");
        return 1;
    }
    if (Mk_qs_thresh<0 || Mk_qs_thresh>Max_thresh) {
        fprintf(stderr,"Illegal Mk_qs_thresh parameter!\n");
        return 1;
    }
    if (Blind_sort_ratio<=0) {
        fprintf(stderr,"blind_sort ratio must be greater than 0!\n");
        return 1;
    }
    return 0;
}




int compute_overshoot(void)
{
    return 9+(Shallow_limit+Cmp_overshoot);
}


void pretty_putchar(int c)
{

    if (c>=32 && c<127)
        printf("  %c", c);
    else if (c=='\n')
        printf(" \\n");
    else if (c=='\t')
        printf(" \\t");
    else
        printf(" %02x", c);
}


int scmp3(unsigned char* p, unsigned char* q, int* l, int maxl)
{
    int i;
    i = 0;
    while (maxl>0 && *p==*q) {
        p++; q++; i++;
        maxl--;
    }
    *l = i;
    if (maxl>0) return *p-*q;
    return q-p;
}






extern int Anchor_dist;
extern int Shallow_limit;
extern int _ds_Verbose;


#define BIGFREQ(b) (ftab[((b)+1) << 8] - ftab[(b) << 8])


/* ------ "local" global variables ------- */
int32_t  Text_size;
uint8_t*  Text;
int32_t*  Sa;
uint8_t*  Upper_text_limit;
int32_t*  Anchor_rank;

uint16_t*  Anchor_offset;

int32_t Anchor_num;
int32_t ftab [65537];
int32_t runningOrder[256];


/* ------------------------------------------------------------------------
   The use of Anchor_rank[] and Anchor_offset is the following:
   Anchor_rank[i] is either -1 or contains the rank (in the list of
   the sorted suffixes of the suffix starting at position
     i*Anchor_dist + Anchor_offset[i].
   Initially Anchor_rank[i] = -1 and Anchor_offset[i]=Anchor_dist,
   then, if a suffix in position t (i*Anchor_dist <= t < (i+1)*Anchor_dist)
   appears to be in a large group which is sorted, the rank of
   t is stored in Anchor_rank[i], and the value t-(i*Anchor_dist)
   is written to Anchor_offset[i]. Both vaulues can be later updated,
   but the value in Anchor_offset[i] can only decrease, so no further
   changes are done when Anchor_offset[i] is = 0. The invariant is:
   if Anchor_rank[i]>=0 then
       Sa[Anchor_rank[i]]=i*Anchor_dist+Anchor_offset[i]
   -------------------------------------------------------------------------*/



void check_ordering(int, int);
void calc_running_order(void);


/* ************************************************************
   This is the main deep/shallow suffix sorting routines
   It divides the suffixes in buckets according to the
   first two characters. Some of the buckets are then sorted
   calling shallow_sort(). When all buckets of kind ax (x!=a) are
   sorted, we use this ordering to sort all suffixes in the
   buckets ya (for any y including y=a).
   ************************************************************* */
void ds_ssort(uint8_t* x, int32_t* p, int32_t n)
{
    void shallow_sort(int32_t*, int, int);
    int compute_overshoot(void), overshoot;
    int32_t  i, j, ss, sb, k;
    uint8_t  c1, c2;
    uint8_t   bigDone[256];
    int32_t  copyStart[256];
    int32_t  copyEnd  [256];
    int32_t  numQSorted = 0;


    Text=x;
    Text_size=n;
    Sa = p;
    Upper_text_limit = Text + Text_size;

    overshoot = compute_overshoot();
    for (i=n; i<n+overshoot; i++) Text[i]=0;


    if (Anchor_dist==0) {
        Anchor_num=0; Anchor_rank=NULL; Anchor_offset=NULL;
    } else {
        Anchor_num = 2 + (n-1)/Anchor_dist;
        Anchor_rank = (int32_t*) malloc(Anchor_num*sizeof(int32_t));
        Anchor_offset = (uint16_t*) malloc(Anchor_num*sizeof(uint16_t));
        if (!Anchor_rank || !Anchor_offset) {
            fprintf(stderr, "malloc failed (ds_sort)\n");
            exit(1);
        }
        for (i=0; i<Anchor_num; i++) {
            Anchor_rank[i]= -1;
            Anchor_offset[i] = Anchor_dist;
        }
    }


    for (i = 0; i <= 65536; i++) ftab[i] = 0;
    c1 = Text[0];
    for (i = 1; i <= Text_size; i++) {
        c2 = Text[i];
        ftab[(c1 << 8) + c2]++;
        c1 = c2;
    }
    for (i = 1; i <= 65536; i++) ftab[i] += ftab[i-1];


    c1 = Text[0];
    for (i = 0; i < Text_size; i++) {
        c2 = Text[i+1];
        j = (c1 << 8) + c2;
        c1 = c2;
        ftab[j]--;
        Sa[ftab[j]] = i;
    }

    /* decide on the running order */
    calc_running_order();
    for (i = 0; i < 256; i++) bigDone[i] = FALSE;

    /* Really do the suffix sorting */
    for (i = 0; i <= 255; i++) {

        /*--
          Process big buckets, starting with the least full.
          --*/
        ss = runningOrder[i];
        if (_ds_Verbose>2)
            fprintf(stderr,"group %3d;  size %d\n",ss,BIGFREQ(ss)&CLEARMASK);

        /*--
          Complete the big bucket [ss] by sorting
          any unsorted small buckets [ss, j].  Hopefully
          previous pointer-scanning phases have already
          completed many of the small buckets [ss, j], so
          we don't have to sort them at all.
          --*/
        for (j = 0; j <= 255; j++) {
            if (j != ss) {
                sb = (ss << 8) + j;
                if (!(ftab[sb] & SETMASK)) {
                    int32_t lo = ftab[sb]   & CLEARMASK;
                    int32_t hi = (ftab[sb+1] & CLEARMASK) - 1;
                    if (hi > lo) {
                        if (_ds_Verbose>2)
                            fprintf(stderr,"sorting [%02x, %02x], done %d "
                                    "this %d\n", ss, j, numQSorted, hi - lo + 1);
                        shallow_sort(Sa+lo, hi-lo+1,Shallow_limit);
#if 0
                        check_ordering(lo, hi);
#endif
                        numQSorted += (hi - lo + 1);
                    }
                }
                ftab[sb] |= SETMASK;
            }
        }
        assert(!bigDone[ss]);

        {
            for (j = 0; j <= 255; j++) {
                copyStart[j] =  ftab[(j << 8) + ss]     & CLEARMASK;
                copyEnd  [j] = (ftab[(j << 8) + ss + 1] & CLEARMASK) - 1;
            }

            if (ss==0) {
                k=Text_size-1;
                c1 = Text[k];
                if (!bigDone[c1])
                    Sa[ copyStart[c1]++ ] = k;
            }
            for (j = ftab[ss << 8] & CLEARMASK; j < copyStart[ss]; j++) {
                k = Sa[j]-1; if (k < 0) continue;
                c1 = Text[k];
                if (!bigDone[c1])
                    Sa[ copyStart[c1]++ ] = k;
            }
            for (j = (ftab[(ss+1) << 8] & CLEARMASK) - 1; j > copyEnd[ss]; j--) {
                k = Sa[j]-1; if (k < 0) continue;
                c1 = Text[k];
                if (!bigDone[c1])
                    Sa[ copyEnd[c1]-- ] = k;
            }
        }
        assert(copyStart[ss] - 1 == copyEnd[ss]);
        for (j = 0; j <= 255; j++) ftab[(j << 8) + ss] |= SETMASK;
        bigDone[ss] = TRUE;
    }
    if (_ds_Verbose) {
        fprintf(stderr, "\t %d pointers, %d sorted, %d scanned\n",
                Text_size, numQSorted, Text_size - numQSorted);
        fprintf(stderr, "\t %d calls to helped_sort\n",Calls_helped_sort);
        fprintf(stderr, "\t %d calls to anchor_sort (forward)\n",
                Calls_anchor_sort_forw);
        fprintf(stderr, "\t %d calls to anchor_sort (backward)\n",
                Calls_anchor_sort_backw);
        fprintf(stderr, "\t %d calls to pseudo_anchor_sort (forward)\n",
                Calls_pseudo_anchor_sort_forw);
        fprintf(stderr, "\t %d calls to deep_sort\n",Calls_deep_sort);
    }

    free(Anchor_offset);
    free(Anchor_rank);
}



/* ****************************************************************
   compute running =(sorting) order for big buckets: start with
   the least full and proceed to the largest one.
   The sorting is done using shellsort
   **************************************************************** */

void calc_running_order(void)
{
    int32_t i, j;
    for (i = 0; i <= 255; i++) runningOrder[i] = i;

    {
        int32_t vv;
        int32_t h = 1;
        do h = 3 * h + 1; while (h <= 256);
        do {
            h = h / 3;
            for (i = h; i <= 255; i++) {
                vv = runningOrder[i];
                j = i;
                while (BIGFREQ(runningOrder[j-h]) > BIGFREQ(vv)) {
                    runningOrder[j] = runningOrder[j-h];
                    j = j - h;
                    if (j <= (h - 1)) goto zero;
                }
zero:
                runningOrder[j] = vv;
            }
        } while (h != 1);
    }
}


uint8_t* transform_bwt(uint8_t* input,int32_t n,uint8_t* out,int32_t* I)
{
    int32_t overshoot;
    int32_t* sa;
    int32_t i,j;

    overshoot=init_ds_ssort(500,2000);

    uint8_t* txt = safe_malloc((n+overshoot)*sizeof(uint8_t));

    memcpy(txt,input,n);

    Text = txt;

    Upper_text_limit = Text + n;
    Text_size = n;

    sa = safe_malloc(n*sizeof(int32_t));

    Sa = sa;

    ds_ssort(txt,sa,n);

    j = 1;
    out[0] = txt[n-1];
    for (i=0; i<n; i++) {
        if (sa[i]!=0) {
            out[j] = txt[sa[i]-1];
            j++;
        } else *I = i;
    }

    free(sa);
    free(txt);

    return out;
}





