/*
 * File:   libpqueue.h
 * Author: matt
 *
 * Created on 29 November 2010, 8:50 PM
 */

#ifndef LIBPQUEUE_H
#define	LIBPQUEUE_H

#ifdef	__cplusplus
extern "C" {
#endif

    typedef struct {
        void* data;
        uint32_t prio;
    } pq_item;

    typedef struct {
        int32_t len;
        int32_t size;
        pq_item** heap;
    } pqueue_t;

    pqueue_t* pqueue_create();
    void pqueue_free(pqueue_t*);
    void pqueue_enqueue(pqueue_t* pq,void* data,uint32_t prio);
    void* pqueue_dequeue(pqueue_t* pq);
    int32_t pqueue_isempty(pqueue_t* pq);
    void pqueue_heapify(pq_item** heap,int32_t n,int32_t i);

#ifdef	__cplusplus
}
#endif

#endif	/* LIBPQUEUE_H */

