/*
 * File:   libpqueue.c
 * Author: Matthias Petri
 *
 * heap based priority queue
 *
 */

#include "libutil.h"
#include "libpqueue.h"

#include <math.h>

pqueue_t* pqueue_create()
{
    pqueue_t* pq = safe_malloc(sizeof(pqueue_t));

    pq->size = 50;
    pq->len = 0;
    pq->heap = safe_malloc(50*sizeof(pq_item*));

    return pq;
}

void pqueue_free(pqueue_t* pq)
{
    if (pq) {
        free(pq->heap);
        free(pq);
    }
}

void pqueue_enqueue(pqueue_t* pq,void* data,uint32_t prio)
{
    pq_item* tmp;
    int32_t parent,i;
    if (!pq) fatal("pq == NULL");



    if (pq->size<=pq->len) {
        pq->heap = safe_realloc(pq->heap,pq->size*2*sizeof(pq_item*));
        pq->size = pq->size *2;
    }

    tmp = safe_malloc(sizeof(pq_item));
    tmp->data = data;
    tmp->prio = prio;

    i = pq->len;
    pq->heap[i] = tmp;
    pq->len++;

    parent = (int32_t) floor((float)i/2);
    while ((i>0) && (pq->heap[parent]->prio > pq->heap[i]->prio)) {
        tmp = pq->heap[parent];
        pq->heap[parent] = pq->heap[i];
        pq->heap[i] = tmp;
        i = parent;
        parent = (int32_t) floor((float)i/2);
    }
}

void*
pqueue_dequeue(pqueue_t* pq)
{
    void* data;
    pq_item* root;
    if (!pq) fatal("pq == NULL");

    if (pq->len <= 0) return NULL;

    root = pq->heap[0];
    data = root->data;
    free(root);

    pq->heap[0] = pq->heap[pq->len-1];
    pq->len--;

    pqueue_heapify(pq->heap,pq->len,0);

    return data;
}

int32_t pqueue_isempty(pqueue_t* pq)
{
    if (pq->len > 0) return 0;
    return 1;
}

void
pqueue_heapify(pq_item** heap,int32_t n,int32_t i)
{
    pq_item* tmp;
    int32_t left,right,min;
    if (!heap) fatal("heap == NULL");

    if (n>1) {
        left = 2*i;
        right = left+1;
        min = i;
        if (left < n && heap[left]->prio < heap[min]->prio)
            min = left;
        if (right < n  && heap[right]->prio < heap[min]->prio)
            min = right;
        if (min != i) {
            tmp = heap[i];
            heap[i] = heap[min];
            heap[min] = tmp;
            pqueue_heapify(heap,n,min);
        }
    }
}