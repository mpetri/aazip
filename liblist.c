/***************************************************************************
 *
 *   File        : liblist.c
 *
 *  List stuff from Assignment 1
 *
 ***************************************************************************/

#include "liblist.h"

/* create an empty list */
list_t*
list_create(void)
{
    list_t* list;

    if ((list = (list_t*) malloc(sizeof(list_t))) == NULL) {
        fprintf(stderr, "Memory allocation for list create failed!\n");
        exit(EXIT_FAILURE);
    }


    list->head = NULL;
    list->len = 0;

    return list;
}

/* free all memory associated with list */
void
list_free(list_t* list)
{
    lnode_t* tmp;
    while (list->head) {
        tmp = list->head->next;
        free(list->head);
        list->head = tmp;
    }
    free(list);
}

/* Insert items */
void
list_insert(list_t* list, int value)
{
    /* any required local vars go here */
    lnode_t*         new, *current, *previous;

    /*
     * allocate memory for a new list element and return FAILURE if
     * malloc() failed.
     */
    if ((new = malloc(sizeof(lnode_t))) == NULL) {
        fprintf(stderr, "Memory allocation for list insert failed!\n");
        exit(EXIT_FAILURE);
    }
    /*
     * Check to see how this function was used in main() to determine what
     * value to return for a FAILURE
     *
     * Here you'll need to write some code to search for the correct place
     * within the linked list to insert the new element. Finally return a
     * value indicating SUCCESS */
    current = list->head;
    previous = NULL;
    while (current != NULL) {
        previous = current;
        current = current->next;
    }

    new->freq = 0;
    new->wfreq = 0.0f;
    new->ts1 = -1;
    new->ts2 = -1;
    new->data = value;
    new->next = current;
    list->len++;

    if (previous == NULL) {
        list->head = new;
        new->prev = NULL;
    } else {
        previous->next = new;
        new->prev = previous;
    }
}

/* find items */
lnode_t*
list_find(list_t* plist, int value,int* cost)
{
    lnode_t* current;

    if (!plist || plist->len == 0) {
        *cost = -1;
        return NULL;
    }

    current = plist->head;

    *cost = 0;
    while (current != NULL && current->data != value) {
        current = current->next;
        *cost = *cost + 1;
    }

    return current;
}

void
list_movetofront(list_t* list,lnode_t* item)
{
    if (item == list->head) return;

    /* unlink item */
    item->prev->next = item->next;
    if (item->next) item->next->prev = item->prev;

    /* move to front */
    item->next = list->head;
    list->head->prev = item;
    item->prev = NULL;
    list->head = item;
}

void
list_moveafter(lnode_t* item,lnode_t* after)
{
    if (after->next == item) return;

    /* unlink item */
    if (item->prev) item->prev->next = item->next;
    if (item->next) item->next->prev = item->prev;

    /* move after */
    item->next = after->next;
    after->next->prev = item;
    after->next = item;
    item->prev = after;
}

void
list_print(list_t* list)
{
    lnode_t* tmp;
    tmp = list->head;

    fprintf(stdout,"L(%d): ",list->len);
    while (tmp != NULL) {
        fprintf(stdout,"%d ",tmp->data);
        tmp = tmp->next;
    }
    fprintf(stdout,"\n");
}

list_t*
list_copy(list_t* list)
{
    list_t* new;
    lnode_t* tmp;

    new = list_create();

    tmp = list->head;

    while (tmp != NULL) {
        list_insert(new,tmp->data);
        tmp = tmp->next;
    }

    return new;
}

void
list_sort(list_t* list, int (*cmp_fptr)(const void*, const void*))
{
    lnode_t* last,*cur,*next,*tmp;

    if (list->head != NULL) cur = list->head->next;
    else return;

    last = list->head;
    last->next = NULL;
    while (cur != NULL) {
        next = cur->next;
        tmp = last;

        while (tmp != NULL && cmp_fptr(tmp,cur) < 0) {
            tmp = tmp->prev;
        }

        if (tmp == NULL) {
            /* start of the list */
            tmp = list->head;
            list->head = cur;
            cur->next = tmp;
            tmp->prev = cur;
            cur->prev = NULL;
        } else {
            /* middle or end */
            if (tmp==last) {
                last = cur;
            }
            cur->next = tmp->next;
            tmp->next = cur;
            cur->prev = tmp;
            if (cur->next) cur->next->prev = cur;
        }
        cur = next;
    }
    last->next = NULL;
}


lnode_t*
list_get(list_t* list,int id)
{
    int32_t i;
    lnode_t* tmp;
    tmp = list->head;

    i = 0;
    while (tmp != NULL) {

        if (i==id) return tmp;

        tmp = tmp->next;
        i++;
    }
    return NULL;
}

