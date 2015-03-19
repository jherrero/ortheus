#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "substitutionC.h"
#include "fastCMaths.h"
#include "commonC.h"

//#include "sequenceGraphC.h"

struct SubModel *constructSubModel(FLOAT_32 *subMatrixForward,
                                   FLOAT_32 *subMatrixBackward,
                                   FLOAT_32 *insertionDistribution,
                                   FLOAT_32 *deletionDistribution) {
    struct SubModel *subModel;
  
    subModel = mallocLocal(sizeof(struct SubModel));
    subModel->subMatrixForward = subMatrixForward;
    subModel->subMatrixBackward = subMatrixBackward;
    subModel->insertionDistribution = insertionDistribution;
    subModel->deletionDistribution = deletionDistribution;
    return subModel;
} 

void destructSubModel(struct SubModel *subModel) {
    free(subModel->subMatrixForward);
    free(subModel->subMatrixBackward);
    free(subModel->insertionDistribution);
    free(subModel->deletionDistribution);
    free(subModel);
}
 
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////      
//library functions for dealing with wV
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
        
void transformWVByDistance(FLOAT_32 *wV, FLOAT_32 *subMatrix, FLOAT_32 *result) {
    //transform wV by given substitution matrix
    INT_32 i;
    FLOAT_32 j;
    FLOAT_32 *k;
    INT_32 l;
    static FLOAT_32 nc[ALPHABET_SIZE];
     
    for(i=0; i<ALPHABET_SIZE; i++) {
        nc[i] = 0.0f;
    } 
   // *nc = [fastMaths.ZERO_PROB]*ALPHABET_SIZE
    k = subMatrix;
    for(i=0; i<ALPHABET_SIZE; i++) { // in xrange(0, ALPHABET_SIZE):
        j = wV[i];
        for(l=0; l< ALPHABET_SIZE; l++) {
            nc[l] += j * k[l];
        }
        k += ALPHABET_SIZE;
    }
    memcpy(result, &nc, sizeof(FLOAT_32)*ALPHABET_SIZE);
}

void copyWV(FLOAT_32 *wVX, FLOAT_32 *wVY) {
    INT_32 i;
    
    for(i=0; i<ALPHABET_SIZE; i++) {
        wVY[i] = wVX[i];
    }
}

void multiplyWV(FLOAT_32 *wVX, FLOAT_32 *wVY, FLOAT_32 *result) {
    INT_32 i;
    
    for(i=0; i<ALPHABET_SIZE; i++) {
        result[i] = wVX[i] * wVY[i];
    }
}

void normaliseWV(FLOAT_32 *wV, FLOAT_32 *result) {
    normaliseWV_GiveFac(wV, result, 1.0f);
}

void normaliseWV_GiveFac(FLOAT_32 *wV, FLOAT_32 *result, FLOAT_32 normFac) {
     //make char probs divisible by one
    INT_32 i;
    FLOAT_32 j;
    
    j = wV[0];
    for(i=1; i<ALPHABET_SIZE; i++) {
        j += wV[i];
    }
    j /= normFac;
    for(i=0; i<ALPHABET_SIZE; i++) {
        result[i] = wV[i] / j;
    }
    //f = sum(wV)
    //return [ i/f for i in wV ]
}

FLOAT_32 combineWV(FLOAT_32 *wVX, FLOAT_32 *wVY) {
    INT_32 i;
    FLOAT_32 j;
    
    j = wVX[0] * wVY[0];;
    for(i=1; i<ALPHABET_SIZE; i++) {
        j += wVX[i] * wVY[i];
    }
    return LOG(j);
}

FLOAT_32 sumWV(FLOAT_32 *wV) {
    INT_32 i;
    FLOAT_32 j;
    
    j = wV[0];
    for(i=1; i<ALPHABET_SIZE; i++) {
        j += wV[i];
    }
    return j;
}

void addWV(FLOAT_32 *wVX, FLOAT_32 *wVY, FLOAT_32 *result) {
    INT_32 i;
    
    for(i=0; i<ALPHABET_SIZE; i++) {
        result[i] = wVX[i] + wVY[i];
    }
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//substitution matrix functions
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

INT_32 valuesAreClose(FLOAT_32 a, FLOAT_32 b, FLOAT_32 tolerance) {
    return a <= b + tolerance && a >= b - tolerance;
}

void checkMatrix(FLOAT_32 *freqVec, FLOAT_32 *matrix) {
    int i;
    int j;
    
    FLOAT_32 v1;
    FLOAT_32 v2;
    
    for(i=0; i<ALPHABET_SIZE; i++) {
        for(j=0; j<ALPHABET_SIZE; j++) {
            //printf("hello %f %f \n", matrix[i*ALPHABET_SIZE + j], matrix[j*ALPHABET_SIZE + i]);
            v1 = freqVec[i] * matrix[i*ALPHABET_SIZE + j];
            v2 = freqVec[j] * matrix[j*ALPHABET_SIZE + i];
            assert(v1 < v2 + 0.0001);
            assert(v2 < v1 + 0.0001);
        }
    }  
}

FLOAT_32 *reverseSubMatrixInPlace(FLOAT_32 *wV) {
    INT_32 i;
    INT_32 j;
    FLOAT_32 k;
    
    for(i=0; i<ALPHABET_SIZE; i++) {
        for(j=i+1; j<ALPHABET_SIZE; j++) {
            k = wV[i*ALPHABET_SIZE + j];
            wV[i*ALPHABET_SIZE + j] = wV[j*ALPHABET_SIZE + i];
            wV[j*ALPHABET_SIZE + i] = k;
        }
    }
    return wV;
}

FLOAT_32 *tamuraNei(FLOAT_32 distance, FLOAT_32 freqA, FLOAT_32 freqC, FLOAT_32 freqG, FLOAT_32 freqT, 
                    FLOAT_32 alphaPurine, FLOAT_32 alphaPyrimidine, FLOAT_32 beta) {
    INT_32 i;
    INT_32 j;
    INT_32 k;
    FLOAT_32 l;
    FLOAT_32 *matrix;
    
    FLOAT_32 freq[] = { freqA, freqC, freqG, freqT };
    FLOAT_32 alpha[] = { alphaPurine, alphaPyrimidine, alphaPurine, alphaPyrimidine };
   
    matrix = mallocLocal(sizeof(FLOAT_32)*ALPHABET_SIZE*ALPHABET_SIZE);
    //see page 203 of Felsenstein's Inferring Phylogenies for explanation of calculations
    for(i=0; i<ALPHABET_SIZE; i++) { //long winded, totally unoptimised method for calculating matrix
        for(j=0; j<ALPHABET_SIZE; j++) {
            l = 0.0f;
            for(k=0; k<ALPHABET_SIZE; k++) {
                l += ((j % 2) == (k % 2) ? 1.0f : 0.0f) * freq[k]; //wat-kro function
            }
            matrix[i*ALPHABET_SIZE + j] = 
            exp(-(alpha[i] + beta) * distance) * (i == j ? 1.0f : 0.0f) //kroenicker delta
            + exp(-beta*distance) * (1.0f - exp(-alpha[i] * distance)) * (freq[j] * ((i % 2) == (j % 2) ? 1.0f : 0.0f) / l)
            + (1.0f - exp(-beta * distance)) * freq[j];
        }
    }
    checkMatrix(freq, matrix);
    return matrix;
}

FLOAT_32 *hKY(FLOAT_32 distance, FLOAT_32 freqA, FLOAT_32 freqC, FLOAT_32 freqG, FLOAT_32 freqT, 
                                FLOAT_32 transitionTransversionRatio) {                                 
   FLOAT_32 alphaPurine;
   FLOAT_32 alphaPyrimidine;
   FLOAT_32 beta;
   
   FLOAT_32 freqPurine; //working parameters
   FLOAT_32 freqPyrimidine;
   FLOAT_32 p;
   
   freqPurine = freqA + freqG;
   freqPyrimidine = freqC + freqT;
   p = freqPurine/freqPyrimidine; //makes like HKY
   //p = 1.0f; //ratio of purine / pyrimidine transition rate parameters -- set to 1.0 for now
   beta = 1.0f / (2.0f * freqPurine * freqPyrimidine * (1.0f + transitionTransversionRatio));
   alphaPyrimidine = ((freqPurine * freqPyrimidine * transitionTransversionRatio) - (freqA * freqG) - (freqC * freqT))
                     / (2.0f * (1.0f + transitionTransversionRatio) * (freqPyrimidine * freqA * freqG * p + freqPurine * freqC * freqT));
   alphaPurine = p * alphaPyrimidine;
   return tamuraNei(distance, freqA, freqC, freqG, freqT, alphaPurine, alphaPyrimidine, beta);
}

FLOAT_32 *subMatrix_HalpernBruno(FLOAT_32 *freqColumn, FLOAT_32 *subMatrix) {
    //return hKY(freqColumn[0], freqColumn[1], freqColumn[2], freqColumn[3], 2.0f;)
    INT_32 i;
    INT_32 j;
    FLOAT_32 a;
    FLOAT_32 b;
    FLOAT_32 *matrix;
    
    matrix = mallocLocal(sizeof(FLOAT_32)*ALPHABET_SIZE*ALPHABET_SIZE);
    for(i=0; i<ALPHABET_SIZE; i++) {
        for(j=0; j<ALPHABET_SIZE; j++) {
            a = freqColumn[i] * subMatrix[i*ALPHABET_SIZE + j];
            b = freqColumn[j] * subMatrix[j*ALPHABET_SIZE + i];
            if (!valuesAreClose(a, b, 0.0001f)) {
                matrix[i*ALPHABET_SIZE + j] = subMatrix[i*ALPHABET_SIZE + j] * (log(b/a) / (1 - (a/b)));
            }
            else {
                matrix[i*ALPHABET_SIZE + j] = subMatrix[i*ALPHABET_SIZE + j];
            }
        }
    }
    return matrix;
}

struct SubModel *constructHKYSubModel(FLOAT_32 distance, 
                                            FLOAT_32 freqA, FLOAT_32 freqC, FLOAT_32 freqG, FLOAT_32 freqT, 
                                            FLOAT_32 transitionTransversionRatio) {
    FLOAT_32 *i;
    FLOAT_32 *j;
    INT_32 k;
    
    FLOAT_32 freq[] = { freqA, freqC, freqG, freqT };
    FLOAT_32 freqNeutral[] = { 0.25, 0.25, 0.25, 0.25 };
    
    k = sizeof(FLOAT_32)*ALPHABET_SIZE;
    i = memcpy(mallocLocal(k), freq, k);
    j = memcpy(mallocLocal(k), freqNeutral, k);
    return constructSubModel(hKY(distance, freqA, freqC, freqG, freqT, transitionTransversionRatio), 
                             reverseSubMatrixInPlace(hKY(distance, freqA, freqC, freqG, freqT, transitionTransversionRatio)), i, j);
}

FLOAT_32 *jukesCantor(FLOAT_32 d) {
    FLOAT_32 i;
    FLOAT_32 j;
    FLOAT_32 *k; 
    INT_32 l;
    
    assert(ALPHABET_SIZE == 4);
    i = 0.25 + 0.75*exp(-(4.0/3.0)*d);
    j = 0.25 - 0.25*exp(-(4.0/3.0)*d);
    k = mallocLocal(sizeof(FLOAT_32)*ALPHABET_SIZE*ALPHABET_SIZE);
    for(l=0; l<ALPHABET_SIZE*ALPHABET_SIZE; l++) {
        k[l] = j;
    }
    for(l=0; l<ALPHABET_SIZE*ALPHABET_SIZE; l+=ALPHABET_SIZE+1) {
        k[l] = i;
    }
    //k = { { i, j, j, j }, { j, i, j, j }, { j, j, i, j }, { j, j, j, i } };
    return k;
}

struct SubModel *constructJukesCantorSubModel(FLOAT_32 distance) {
    FLOAT_32 *i;
    FLOAT_32 *j;
    INT_32 k;
    
    i = mallocLocal(sizeof(FLOAT_32)*ALPHABET_SIZE);
    j = mallocLocal(sizeof(FLOAT_32)*ALPHABET_SIZE);
    for(k=0; k<ALPHABET_SIZE; k++) {
        i[k] = 0.25f;
        j[k] = 0.25f;
    }
    return constructSubModel(jukesCantor(distance), 
                             jukesCantor(distance), i, j);
}

FLOAT_32 * dNAMap_IUPACToWVFn(char i) {
    static FLOAT_32 j[ALPHABET_SIZE];
    INT_32 k;
    
    for(k=0; k<ALPHABET_SIZE; k++) {
        j[k] = 0.0f;
    }
    switch(i) {
        case 'A':
        case 'a':
            j[0] = 1.0f;
            break;
        case 'C':
        case 'c':
            j[1] = 1.0f;
            break;
        case 'G':
        case 'g':
            j[2] = 1.0f;
            break;
        case 'T':
        case 't':
            j[3] = 1.0f;
            break;
        default:
            for(k=0; k<ALPHABET_SIZE; k++) {
                j[k] = 0.25f;
            }
            break;
   }
   assert(sumWV(j) < 1.01);
   assert(sumWV(j) > 0.99);
   return j;
}


