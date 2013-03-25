/*
 * File:   liblupdate.h
 * Author: matt
 *
 */

#ifndef LIBLUPDATE_H
#define	LIBLUPDATE_H

#ifdef	__cplusplus
extern "C" {
#endif

    uint8_t* lupdate_simple(uint8_t* bwt,uint32_t size,uint8_t* input,uint64_t* cost);

    uint8_t* lupdate_movetofront(uint8_t* bwt,uint32_t size,uint8_t* input,uint64_t* cost);

    uint8_t* lupdate_freqcount(uint8_t* bwt,uint32_t size,uint8_t* input,uint64_t* cost);

    uint8_t* lupdate_wfc(uint8_t* bwt,uint32_t size,uint8_t* input,uint64_t* cost);

    uint8_t* lupdate_timestamp(uint8_t* bwt,uint32_t size,uint8_t* input,uint64_t* cost);


#ifdef	__cplusplus
}
#endif

#endif	/* LIBLUPDATE_H */

