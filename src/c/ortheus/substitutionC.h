#ifndef SUBSTITUTIONC_H_
#define SUBSTITUTIONC_H_

#include <math.h>

#include "fastCMaths.h"
#include "commonC.h"

//subtitution stuff
#define GAP 1000
//just dna right now
#define ALPHABET_SIZE 4 

struct SubModel {
   FLOAT_32 *subMatrixForward;
   FLOAT_32 *subMatrixBackward;
   FLOAT_32 *insertionDistribution;
   FLOAT_32 *deletionDistribution;
};

struct SubModel *constructSubModel(FLOAT_32 *subMatrixForward,
                                   FLOAT_32 *subMatrixBackward,
                                   FLOAT_32 *insertionDistribution,
                                   FLOAT_32 *deletionDistribution);

void destructSubModel(struct SubModel *subModel);

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////      
//library functions for dealing with wV
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

void copyWV(FLOAT_32 *wVX, FLOAT_32 *wVY);
        
void transformWVByDistance(FLOAT_32 *wV, FLOAT_32 *subMatrix, FLOAT_32 *result);

void multiplyWV(FLOAT_32 *wVX, FLOAT_32 *wVY, FLOAT_32 *result);

void normaliseWV_GiveFac(FLOAT_32 *wV, FLOAT_32 *result, FLOAT_32 normFac);

void normaliseWV(FLOAT_32 *wV, FLOAT_32 *result); 

FLOAT_32 combineWV(FLOAT_32 *wVX, FLOAT_32 *wVY);

FLOAT_32 sumWV(FLOAT_32 *wV);

void addWV(FLOAT_32 *wVX, FLOAT_32 *wVY, FLOAT_32 *result);

FLOAT_32 * dNAMap_IUPACToWVFn(char i);

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//substitution matrix functions
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

FLOAT_32 *subMatrix_HalpernBruno(FLOAT_32 *freqColumn, FLOAT_32 *m);

FLOAT_32 *hKY(FLOAT_32 distance, FLOAT_32 freqA, FLOAT_32 freqC, FLOAT_32 freqG, FLOAT_32 freqT, 
                                FLOAT_32 transitionTransversionRatio);

struct SubModel *constructHKYSubModel(FLOAT_32 distance, FLOAT_32 freqA, FLOAT_32 freqC, FLOAT_32 freqG, FLOAT_32 freqT, 
                                            FLOAT_32 transitionTransversionRatio);

FLOAT_32 *jukesCantor(FLOAT_32 d);

struct SubModel *constructJukesCantorSubModel(FLOAT_32 distance);

FLOAT_32 *reverseSubMatrixInPlace(FLOAT_32 *wV);

#endif /*SUBSTITUTIONC_H_*/
