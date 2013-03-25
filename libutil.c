/*
 * File:   libutil.c
 * Author: Matthias Petri
 *
 * misc support functions
 *
 */

#include "libutil.h"


/*  safe_malloc ()
 *
 *  Call calloc () and abort if the specified amount of
 *  memory cannot be allocated.
 */

void*
safe_malloc(size_t size)
{
    void* mem_block = NULL;

    if ((mem_block = calloc(1, size)) == NULL) {
        fprintf(stderr, "ERROR: safe_malloc(%zu) cannot allocate memory.", size);
        exit(EXIT_FAILURE);
    }
    return (mem_block);
}

/*  safe_realloc ()
 *
 *  Call realloc () and abort if the specified amount of
 *  memory cannot be allocated.
 */

void*
safe_realloc(void* old_mem, size_t new_size)
{
    if ((old_mem = realloc(old_mem, new_size)) == NULL) {
        fprintf(stderr, "ERROR: safe_realloc() cannot allocate"
                "%zu blocks of memory.\n", new_size);
        exit(EXIT_FAILURE);
    }
    return (old_mem);
}

/*
 * safe_strdup ()
 *
 * Safe version of strdup avoid buffer overflow, etc.
 *
 */

char*
safe_strdup(const char* str)
{
    char* copy = NULL;

    if (str == NULL) {
        fprintf(stderr, "ERROR safe_strdup(): str == NULL");
        exit(EXIT_FAILURE);
    }

    copy = safe_malloc((strlen(str) + 1) * sizeof(char));

    (void) strcpy(copy, str);

    return (copy);
}

char*
safe_strcat(char* str1, const char* str2)
{
    char* rv = NULL;
    size_t len = 0;

    if (str1 == NULL || str2 == NULL) {
        fprintf(stderr, "ERROR safe_strcat_new(): str == NULL");
        exit(EXIT_FAILURE);
    }
    len = strlen(str1) + strlen(str2) + 1;
    rv = safe_malloc(len * sizeof(char));
    (void) strcpy(rv, str1);
    rv = strcat(rv, str2);
    return (rv);
}


FILE*
safe_fopen(const char* filename,const char* mode)
{
    FILE* f;
    f = fopen(filename,mode);
    if (f == NULL) {
        fprintf(stderr, "ERROR: safe_fopen(%s,%s)\n", filename,mode);
        exit(EXIT_FAILURE);
    }
    return (f);
}

void
safe_fclose(FILE* f)
{
    if (fclose(f) != 0) {
        perror("Error: file close():");
        exit(EXIT_FAILURE);
    }
}

int
safe_filesize(FILE* f)
{
    int size;
    int cur = ftell(f);
    if (cur == -1) {
        perror("Error: file ftell():");
        exit(EXIT_FAILURE);
    } else {
        fseek(f,0,SEEK_END);
        size = ftell(f);
        if (size == -1) {
            perror("Error: file ftell():");
            exit(EXIT_FAILURE);
        }
        fseek(f,cur,SEEK_SET);
    }
    return size;
}

void
fatal(const char* format, ...)
{
    va_list vargs;

    va_start(vargs, format);
    vfprintf(stderr, format, vargs);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}

uint64_t
gettime()
{
    struct timeval tp;

    gettimeofday(&tp,NULL);

    return (tp.tv_sec*1000000) + tp.tv_usec;
}
