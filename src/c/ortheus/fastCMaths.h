/*This code is taken and adapted directly from the probcons library,
 * with an addition for random numbers.
 */
 
/////////////////////////////////////////////////////////////////
// FLOAT_32.h
//
// Routines for doing math operations in PROBCONS.
/////////////////////////////////////////////////////////////////

#ifndef FASTCMATHS_H
#define FASTCMATHS_H

#include <math.h>
#include <float.h>
#include <limits.h>  
#include <stdlib.h>

#define TRUE 1
#define FALSE 0 

#define LOG_ZERO -1e30000
//-2e20
#define LOG_ONE 0.0
#define ZERO 0.0
#define ONE 1.0

/*const FLOAT_32 LOG_ZERO = -2e20;
const FLOAT_32 LOG_ONE = 0.0;
const FLOAT_32 ZERO = 0.0;
const FLOAT_32 ONE = 1.0;*/

typedef long int INT_32; 
typedef unsigned long int UNSIGNED_INT_32;
typedef long long LONG_64;
typedef unsigned long long UNSIGNED_LONG_64;
typedef float FLOAT_32;
typedef double DOUBLE_64;
 
#define INT_STRING "%li"
#define LONG_INT_STRING "%lli"
#define FLOAT_STRING "%f"
#define DOUBLE_STRING "%f"

#define INT_32_MAX LONG_MAX
#define INT_32_MIN LONG_MIN
#define LONG_64_MAX 9223372036854775807LL
#define LONG_64_MIN -9223372036854775805LL
#define FLOAT_32_MAX FLT_MAX
#define FLOAT_32_MIN FLT_MIN

#define SMALL_CHUNK_SIZE 100
#define MEDIUM_CHUNK_SIZE 1000
#define LARGE_CHUNK_SIZE 1000000

/////////////////////////////////////////////////////////////////
// LOG()
//
// Compute the logarithm of x.
/////////////////////////////////////////////////////////////////

FLOAT_32 LOG (FLOAT_32 x);

/////////////////////////////////////////////////////////////////
// EXP()
//
// Computes exp(x).
/////////////////////////////////////////////////////////////////

inline FLOAT_32 EXP (FLOAT_32 x);

#define EXP_UNDERFLOW_THRESHOLD -4.6
#define LOG_UNDERFLOW_THRESHOLD 7.5

//const FLOAT_32 EXP_UNDERFLOW_THRESHOLD = -4.6;
//const FLOAT_32 LOG_UNDERFLOW_THRESHOLD = 7.5;

/////////////////////////////////////////////////////////////////
// LOOKUP()
//
// Computes log (exp (x) + 1), for 0 <= x <= 7.5.
/////////////////////////////////////////////////////////////////

inline FLOAT_32 LOOKUP (FLOAT_32 x);

/////////////////////////////////////////////////////////////////
// LOG_PLUS_EQUALS()
//
// Add two log probabilities and store in the first argument
/////////////////////////////////////////////////////////////////

inline void LOG_PLUS_EQUALS (FLOAT_32 *x, FLOAT_32 y);


/////////////////////////////////////////////////////////////////
// LOG_ADD()
//
// Add two log probabilities
/////////////////////////////////////////////////////////////////

inline FLOAT_32 LOG_ADD (FLOAT_32 x, FLOAT_32 y);


/////////////////////////////////////////////////////////////////
// LOG_ADD()
//
// Add three log probabilities
/////////////////////////////////////////////////////////////////

inline FLOAT_32 LOG_ADD_THREE (FLOAT_32 x1, FLOAT_32 x2, FLOAT_32 x3);

/////////////////////////////////////////////////////////////////
// MAX_EQUALS()
//
// Chooses maximum of two arguments and stores it in the first argument
/////////////////////////////////////////////////////////////////

inline void MAX_PLUS_EQUALS (FLOAT_32 *x, FLOAT_32 y);

/////////////////////////////////////////////////////////////////
// RANDOM()
//
// Return random FLOAT_32 in range [0 - 1.0 }
/////////////////////////////////////////////////////////////////
inline FLOAT_32 RANDOM();

/////////////////////////////////////////////////////////////////
// RANDOM()
//
// Return random FLOAT_32 in range [0 - 1.0 }
/////////////////////////////////////////////////////////////////
inline FLOAT_32 RANDOM_LOG();

#endif
