/*
 * File:   libbwt.h
 * Author: Matthias Petri
 *
 */

#ifndef LIBBWT_H
#define	LIBBWT_H

#ifdef	__cplusplus
extern "C" {
#endif

#define Cmp_overshoot 16
#define Max_thresh 30

#define SETMASK (1 << 30)
#define CLEARMASK (~(SETMASK))
#define IS_SORTED_BUCKET(sb) (ftab[sb] & SETMASK)
#define BUCKET_FIRST(sb) (ftab[sb]&CLEARMASK)
#define BUCKET_LAST(sb) ((ftab[sb+1]&CLEARMASK)-1)
#define BUCKET_SIZE(sb) ((ftab[sb+1]&CLEARMASK)-(ftab[sb]&CLEARMASK))

    uint8_t* transform_bwt(uint8_t* input,int32_t n,uint8_t* out,int32_t* I);
    uint8_t* reverse_bwt(uint8_t* in,int32_t n,int32_t I,uint8_t* out);


    void helped_sort(int32_t* a, int32_t n, int32_t depth);
    void shallow_inssort_lcp(int32_t* a, int32_t n, uint8_t* text_depth);

    void shallow_mkq(int32_t* a, int n, uint8_t* text_depth);
    void shallow_mkq16(int32_t* a, int n, uint8_t* text_depth);
    void shallow_mkq32(int32_t* a, int n, uint8_t* text_depth);

    void ds_ssort(unsigned char* t, int* sa, int n);
    int init_ds_ssort(int adist, int bs_ratio);

#ifdef	__cplusplus
}
#endif

#endif	/* LIBBWT_H */

