/*
 * File:   libhuff.h
 * Author: matt
 *
 * Created on 28 November 2010, 3:28 PM
 */

#ifndef LIBHUFF_H
#define	LIBHUFF_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "libutil.h"

    struct hnode {
        int32_t sym;
        int32_t freq;
        struct hnode* left;
        struct hnode* right;
    };

    typedef struct hnode hnode_t;

    void encode_huffman(uint8_t* input,uint32_t size,bit_file_t* of);

#ifdef	__cplusplus
}
#endif

#endif	/* LIBHUFF_H */

