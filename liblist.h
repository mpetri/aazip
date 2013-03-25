/*
 * File:   liblist.h
 * Author: matt
 *
 * Created on 28 November 2010, 3:27 PM
 */

#ifndef LIBLIST_H
#define	LIBLIST_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "libutil.h"

    /* Link List data types */
    typedef struct lnode lnode_t;

    struct lnode {
        int             ts1;
        int             ts2;
        int             freq;
        float          wfreq;
        int             data;
        lnode_t*        next;
        lnode_t*        prev;
    };

    typedef struct {
        lnode_t*        head;
        int             len;
    } list_t;

    list_t*         list_create(void);
    void            list_free(list_t* list);
    list_t*         list_copy(list_t* list);
    void            list_insert(list_t* list,int value);
    lnode_t*        list_find(list_t* list,int item,int* cost);
    void            list_print(list_t*);
    void            list_sort(list_t*,int (*cmp_fptr)(const void*, const void*));
    lnode_t*        list_get(list_t*,int id);

    void list_movetofront(list_t* list,lnode_t* item);
    void list_moveafter(lnode_t* item,lnode_t* after);


#ifdef	__cplusplus
}
#endif

#endif	/* LIBLIST_H */

