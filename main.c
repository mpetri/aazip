/*
 * File:   main.c
 * Author: Matthias Petri
 *
 * main aazip
 */

#include "libutil.h"
#include "libbwt.h"
#include "libhuff.h"
#include "liblupdate.h"

enum mode_t {
    UNKNOWN,
    SIMPLE,
    MTF,
    FC,
    WFC,
    TS
};

static void
print_usage(const char* program)
{
    fprintf(stderr, "USAGE: %s -m [algorithm] <input>\n", program);
    fprintf(stderr, "  -m algorithm [simple, mtf, fc, wfc, timestamp]\n");
    fprintf(stderr, "  -h Display usage information\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "EXAMPLE: %s -m mtf test.dat\n",
            program);
    fprintf(stderr, "\n");
    return;
}

/*
 * aazip - compress files using a transform based compression system
 */
int main(int argc, char** argv)
{
    FILE* f;
    bit_file_t* of;
    char* infile,*outfile;
    uint8_t* input,*lupdate,*bwt,lumode;
    int32_t I,osize,opt;
    uint32_t size;
    mode_t lupdate_alg;
    float ient,oent;
    uint64_t cost,tstart,tstop,elapsed;

    /* parse command line parameter */
    opt = GETOPT_FINISHED;
    if (argc <= 1) {
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    while ((opt = getopt(argc, argv, "m:h")) != GETOPT_FINISHED) {
        switch (opt) {
            case 'm':
                if (strcmp(optarg, "simple") == 0) lupdate_alg = SIMPLE;
                else if (strcmp(optarg, "mtf") == 0) lupdate_alg = MTF;
                else if (strcmp(optarg, "fc") == 0) lupdate_alg = FC;
                else if (strcmp(optarg, "wfc") == 0) lupdate_alg = WFC;
                else if (strcmp(optarg, "timestamp") == 0) lupdate_alg = TS;
                else fatal("ERROR: mode <%s> unknown!\n", optarg);
                break;
            case 'h':
            default:
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    /* read input file name */
    if (optind < argc) infile = argv[optind];
    else {
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    /* read input file */
    f = safe_fopen(infile,"r");
    size = safe_filesize(f);
    input = (uint8_t*) safe_malloc(size+1);
    if (fread(input,1,size,f)!=(size_t)size) {
        fatal("read input file.");
    }
    safe_fclose(f);
    input[size] = 0;

    /* TODO calculate input entropy */
    ient = 0.0f;

    /* perform bwt */
    bwt = (uint8_t*) safe_malloc(size);

    tstart = gettime();

    bwt = transform_bwt(input,size,bwt,&I);

    /* peform list update */
    switch (lupdate_alg) {
        case SIMPLE:
            fprintf(stdout,"ALGORITHM: simple\n");
            lupdate = lupdate_simple(bwt,size,input,&cost);
            break;
        case MTF:
            fprintf(stdout,"ALGORITHM: move to front\n");
            lupdate = lupdate_movetofront(bwt,size,input,&cost);
            break;
        case FC:
            fprintf(stdout,"ALGORITHM: frequency count\n");
            lupdate = lupdate_freqcount(bwt,size,input,&cost);
            break;
        case WFC:
            fprintf(stdout,"ALGORITHM: weighted frequency count\n");
            lupdate = lupdate_wfc(bwt,size,input,&cost);
            break;
        case TS:
            fprintf(stdout,"ALGORITHM: timestamp\n");
            lupdate = lupdate_timestamp(bwt,size,input,&cost);
            break;
        default:
            fatal("unkown list update algorithm.");
    }

    fprintf(stdout,"INPUT: %s (%d bytes)\n",infile,size);
    fprintf(stdout,"COST: %lu\n",cost);

    /* TODO calculate entropy after list update*/
    oent = 0.0f;

    /* write output */
    outfile = safe_strcat(infile,".aazip");
    /* create bit file for writing */
    of = BitFileOpen(outfile, BF_WRITE);

    /* write aa zip header */
    BitFilePutChar('A', of);
    BitFilePutChar('A', of);

    /* write I */
    BitFilePutBitsInt(of,&I,32,sizeof(uint32_t));

    /* write lupdate mode */
    lumode = lupdate_alg;
    BitFilePutBitsInt(of,&lumode,8,sizeof(uint8_t));


    fprintf(stderr,"I %d lumode %d\n",I,lumode);

    /* perform huffman coding */
    encode_huffman(lupdate,size,of);

    tstop = gettime();

    elapsed = tstop - tstart;
    fprintf(stdout,"TIME: %.3f s\n",(float)elapsed/1000000);

    /* flush and get file stats */
    BitFileFlushOutput(of,0);
    f = BitFileToFILE(of);
    osize = ftell(f);

    fprintf(stdout,"OUTPUT: %s\n",outfile);
    fprintf(stdout,"ENTROPY: %.2f bps / %.2f bps\n",ient,oent);
    fprintf(stdout,"COMPRESSION: %.2f\n",((float)osize/(float)size)*100);

    /* clean up*/
    safe_fclose(f);
    free(input);
    free(bwt);

    return (EXIT_SUCCESS);
}

