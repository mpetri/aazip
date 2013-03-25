
#include "libutil.h"
#include "liblist.h"
#include "liblupdate.h"

list_t* lupdate_createlist()
{
    uint32_t i;
    list_t* lst;

    lst = list_create();

    for (i=0; i<ALPHABET_SIZE; i++) {
        list_insert(lst,i);
    }

    return lst;
}

uint8_t*
lupdate_simple(uint8_t* bwt,uint32_t size,uint8_t* output,uint64_t* cost)
{
    uint32_t i;

    *cost = 0;
    for (i=0; i<size; i++) {
        output[i] = bwt[i];
        *cost = *cost + 1;
    }

    return output;
}

uint8_t*
lupdate_movetofront(uint8_t* bwt,uint32_t size,uint8_t* output,uint64_t* c)
{
    uint32_t i,chr;
    int32_t cost;
    list_t* lst;
    lnode_t* cur;

    lst = lupdate_createlist();

    *c = 0;
    for (i=0; i<size; i++) {
        chr = bwt[i];

        cur = list_find(lst,chr,&cost);
        list_movetofront(lst,cur);

        *c += cost + 1;

        output[i] = (uint8_t) cost;
    }

    list_free(lst);

    return output;
}

uint8_t*
lupdate_freqcount(uint8_t* bwt,uint32_t size,uint8_t* output,uint64_t* c)
{
    uint32_t chr,i;
    int32_t cost;
    list_t* lst;
    lnode_t* found,*tmp;

    lst = lupdate_createlist();

    *c = 0;
    for (i=0; i<size; i++) {
        chr = bwt[i];

        found = list_find(lst,chr,&cost);

        output[i] = (uint8_t) cost;
        *c += cost;

        /* reorder */
        found->freq++;
        tmp = found->prev;
        while (tmp != NULL && found->freq > tmp->freq) {
            tmp = tmp->prev;
        }

        if (tmp == NULL) {
            list_movetofront(lst,found);
        } else {
            list_moveafter(found,tmp);
        }
    }

    list_free(lst);

    return output;
}

int
wfc_cmp(const void* t1,const void* t2)
{
    lnode_t* n1 = (lnode_t*)t1;
    lnode_t* n2 = (lnode_t*)t2;

    if (n1->wfreq < n2->wfreq)
        return -1;
    else if (n1->wfreq > n2->wfreq)
        return 1;
    else
        return 0;
}


/*
 * weighted frequency count comparison function
 */
float
calc_wfc(int32_t t,int32_t p)
{
    if (t==1) return 1;
    if (t > 1 && t <= 64) return (1.0f/((float)t*(float)p));
    if (t > 64 && t <= 256) return (1.0f/(2.0f*((float)t*(float)p)));
    if (t > 256 && t <= 1024) return (1.0f/(4.0f*((float)t*(float)p)));
    if (t > 1024 && t <= 2048) return (1.0f/(8.0f*((float)t*(float)p)));
    return 0;
}

uint8_t*
lupdate_wfc(uint8_t* bwt,uint32_t size,uint8_t* output,uint64_t* c)
{
    uint32_t chr,j;
    int32_t cost,i;
    list_t* lst;
    lnode_t* found,*tmp;
    int32_t k,start;

    lst = lupdate_createlist();

    *c = 0;
    for (j=0; j<size; j++) {
        chr = bwt[j];

        found = list_find(lst,chr,&cost);

        output[j] = (uint8_t) cost;
        *c += cost;

        /* reset wfreq values */
        tmp = lst->head;
        while (tmp != NULL) {
            tmp->wfreq = 0;
            tmp = tmp->next;
        }

        /* calculate new wfreq values */
        start = MAX(j-512,0);
        for (k=start; k<j; k++) {
            found = list_find(lst,bwt[k],&i);
            found->wfreq = found->wfreq + calc_wfc(j-k,k);
        }

        /* sort list based on wfreq values */
        list_sort(lst,wfc_cmp);
    }

    list_free(lst);

    return output;
}

uint8_t*
lupdate_timestamp(uint8_t* bwt,uint32_t size,uint8_t* output,uint64_t* c)
{
    uint32_t chr,i;
    int32_t cost;
    list_t* lst;
    int ts;
    lnode_t* found;
    lnode_t* tmp;

    lst = lupdate_createlist();

    ts = 0;
    *c = 0;
    for (i=0; i<size; i++) {
        chr = bwt[i];

        found = list_find(lst,chr,&cost);

        output[i] = (uint8_t) cost;
        *c += cost;

        if (found->ts1 != -1) {
            tmp = lst->head;
            while (tmp != NULL && tmp != found) {
                if (tmp->ts1 < found->ts1 || (tmp->ts1 > found->ts1 && found->ts1 > tmp->ts2)) {
                    /* move in front of ts */
                    tmp = tmp->prev;
                    if (tmp == NULL) {
                        /* move to the front of the list */
                        list_movetofront(lst,found);
                    } else {
                        /* unlink found */
                        if (found->prev) found->prev->next = found->next;
                        if (found->next) found->next->prev = found->prev;

                        /* move after tmp */
                        found->next = tmp->next;
                        tmp->next->prev = found;
                        tmp->next = found;
                        found->prev = tmp;
                    }

                    break;
                }
                tmp = tmp->next;
            }
        }

        /* update ts */
        found->ts2 = found->ts1;
        found->ts1 = ts;
        ts++;
    }

    list_free(lst);

    return output;
}

