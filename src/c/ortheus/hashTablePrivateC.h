/* Copyright (C) 2002, 2004 Christopher Clark <firstname.lastname@cl.cam.ac.uk> */

#ifndef __HASHTABLE_PRIVATE_CWC22_H__
#define __HASHTABLE_PRIVATE_CWC22_H__
 
#include "hashTableC.h"

#include "fastCMaths.h"
#include "commonC.h"  

/*****************************************************************************/

struct entry
{
    void *k, *v;
    UNSIGNED_INT_32 h;
    struct entry *next;
};

struct hashtable {
    UNSIGNED_INT_32 tablelength;
    struct entry **table;
    UNSIGNED_INT_32 entrycount;
    UNSIGNED_INT_32 loadlimit;
    UNSIGNED_INT_32 primeindex;
    UNSIGNED_INT_32 (*hashfn) (void *k);
    INT_32 (*eqfn) (void *k1, void *k2);
    void (*keyFree)(void *);
    void (*valueFree)(void *);
    struct Chunks *entryChunks;
};

/*****************************************************************************/
UNSIGNED_INT_32
hash(struct hashtable *h, void *k);

/*****************************************************************************/
/* indexFor */
static inline UNSIGNED_INT_32
indexFor(UNSIGNED_INT_32 tablelength, UNSIGNED_INT_32 hashvalue) {
    return (hashvalue % tablelength);
}

/* Only works if tablelength == 2^N */
/*static inline UNSIGNED_INT_32
indexFor(UNSIGNED_INT_32 tablelength, UNSIGNED_INT_32 hashvalue)
{
    return (hashvalue & (tablelength - 1u));
}
*/

/*****************************************************************************/
/*#define freekey(X) free(X)*/
/*define freekey(X) ; */


/*****************************************************************************/

#endif /* __HASHTABLE_PRIVATE_CWC22_H__*/

/*
 * Copyright (c) 2002, Christopher Clark
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 
 * * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * * Neither the name of the original author; nor the names of any contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 * 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
