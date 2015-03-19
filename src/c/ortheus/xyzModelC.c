#include <stdio.h>
#include <assert.h>

#include "fastCMaths.h"
#include "substitutionC.h"

#include "sequenceGraphC.h"
#include "xyzModelC.h"

#define STATE_NO 12

const INT_32 S1 = 4;
const INT_32 S2 = 2;
const INT_32 S3 = 9; //2;
const INT_32 S4 = 10; //3;
const INT_32 S5 = 6;
const INT_32 S6 = 5;
const INT_32 S7 = 7; //11;
const INT_32 S8 = 3; //7;
const INT_32 S9 = 11; //7;
const INT_32 S10 = 8;
const INT_32 S11 = 1;
const INT_32 S12 = 0;
  
/*
const INT_32 MATCH_COMB_STATE = 0;
const INT_32 INSERT_X_COMB_STATE = 1; //i.e. x,- -- x inserted
const INT_32 INSERT_Y_COMB_STATE = 2; //i.e. -,y -- y inserted
const INT_32 INSERT_X2_COMB_STATE = 3; //i.e. x,- -- x inserted
const INT_32 INSERT_Y2_COMB_STATE = 4; //i.e. -,y -- y inserted
const INT_32 DELETE_X_COMB_STATE = 5; //i.e. -,y -- x deleted
const INT_32 DELETE_Y_COMB_STATE = 6; //i.e. x,- -- y deleted
const INT_32 SILENT_IN_COMB_STATE = 7;
const INT_32 SILENT_OUT_COMB_STATE = 8;
const INT_32 DELETE_COMB_STATE = 9;*/
 
inline INT_32 isSilent(INT_32 state) { return state == S3 || state == S4 || state == S9 || state == S10; }
inline INT_32 isXInsert(INT_32 state) { return state == S2 || state == S11; }
inline INT_32 isXDelete(INT_32 state) { return state == S8; }
inline INT_32 isYInsert(INT_32 state) { return state == S1 || state == S12; }
inline INT_32 isYDelete(INT_32 state) { return state == S6; }
inline INT_32 isMatch(INT_32 state) { return state == S5; }
inline INT_32 isXYDelete(INT_32 state) { return state == S7; }

void calculateExpProb(DOUBLE_64 a, DOUBLE_64 A, DOUBLE_64 b, DOUBLE_64 B, DOUBLE_64 g, DOUBLE_64 G, 
					   DOUBLE_64 a1M, DOUBLE_64 A1M, DOUBLE_64 g1M, DOUBLE_64 G1M,
					   INT_32 state, FLOAT_32 *summedStates) {
	a = exp(a);
	A = exp(A);
	b = exp(b);
	B = exp(B);
	g = exp(g);
	G = exp(G);
	a1M = exp(a1M);
	A1M = exp(A1M);
	g1M = exp(g1M);
	G1M = exp(G1M);
		
	INT_32 MAX_ITERATIONS = 100;
	DOUBLE_64 *states;
	DOUBLE_64 *newStates;
	DOUBLE_64 *summedStates2;
	INT_32 i;
	INT_32 iteration;
	DOUBLE_64 *j;
	
	states = mallocLocal(sizeof(DOUBLE_64)*JUNCTION_STATE_NO);
	newStates = mallocLocal(sizeof(DOUBLE_64)*JUNCTION_STATE_NO);
	summedStates2 = mallocLocal(sizeof(DOUBLE_64)*JUNCTION_STATE_NO);
	
	for(i=0; i<JUNCTION_STATE_NO; i++) {
	    states[i] = 0.0f;
	    newStates[i] = 0.0f;
	    summedStates2[i] = 0.0f;
	}
	states[state] = 1.0f;
	summedStates2[state] = 1.0f;
	//junction-tree conversion 7-0, 10-1, 9-2, 3-3, 4-4
	for(iteration=0; iteration<MAX_ITERATIONS; iteration++) {
		newStates[0] = g*G*states[0] + B*states[1] + b*states[2] + b*B*states[4];
		newStates[1] = A1M * g * G1M * newStates[0];
		newStates[2] = a1M * G * g1M * newStates[0];
		newStates[3] = g1M * G1M * newStates[0];
		newStates[4] = A1M * a1M * newStates[3];
		for(i=0; i<JUNCTION_STATE_NO; i++) {
			summedStates2[i] += newStates[i];
		}
		j = states;
		states = newStates;
		newStates = j;
	}
	free(states);
	free(newStates);
	for(i=0; i<JUNCTION_STATE_NO; i++) {
		summedStates[i] = log(summedStates2[i]);
	}
	//uglyf(" " FLOAT_STRING " " FLOAT_STRING " " FLOAT_STRING " " FLOAT_STRING " " FLOAT_STRING " " FLOAT_STRING " \n ", a, A, b, B, g, G);
	//for(i=0; i<JUNCTION_STATE_NO; i++) {
	//	uglyf(" calc %f %f \n", summedStates[i], EXP(summedStates[i]));
	//}
}

struct CombinedTransitionModel *constructCombinedTransitionModel(FLOAT_32 DELETE_OPEN_PROB, //beta
                                                                 FLOAT_32 DELETE_CONTINIUE_PROB, //gamma
                                                                 FLOAT_32 INSERT_OPEN_PROB, //alpha
                                                                 FLOAT_32 INSERT_CONTINUE_PROB, //delta
                                                                 FLOAT_32 distanceX, FLOAT_32 distanceY, 
                                                                 //FLOAT_32 lengthX, FLOAT_32 lengthY, 
                                                                 struct SubModel *subModelX, 
                                                                 struct SubModel *subModelY) {
    struct CombinedTransitionModel *temp = mallocLocal(sizeof(struct CombinedTransitionModel));
    FLOAT_32 aX, bX, gX, dX, eX;
    FLOAT_32 aX1M, bX1M, gX1M, dX1M, eX1M;
     
    FLOAT_32 aY, bY, gY, dY, eY;
    FLOAT_32 aY1M, bY1M, gY1M, dY1M, eY1M;
  
    aX = INSERT_OPEN_PROB*distanceX; 
    if (aX >= 1.0f) { aX = 0.99f; }
    assert(aX <= 1.0f);
    bX = DELETE_OPEN_PROB*distanceX; 
    if (bX >= 1.0f) { bX = 0.99f; }
    assert(bX <= 1.0f);
    gX = DELETE_CONTINIUE_PROB; 
    dX = INSERT_CONTINUE_PROB; 
    eX = 1.0f - 1.0f/1000; 
    
    aY = INSERT_OPEN_PROB*distanceY; 
    if (aY >= 1.0f) { aY = 0.99f; }
    assert(aY <= 1.0f);
    bY = DELETE_OPEN_PROB*distanceY; 
    if (bY >= 1.0f) { bY = 0.99f; }
    assert(bY <= 1.0f);
    gY = DELETE_CONTINIUE_PROB; 
    dY = INSERT_CONTINUE_PROB;
    eY = 1.0f - 1.0f/1000; 
    
    aX1M = LOG(1.0f - aX);
    bX1M = LOG(1.0f - bX);
    gX1M = LOG(1.0f - gX);
    dX1M = LOG(1.0f - dX);
    eX1M = LOG(1.0f - eX);
    
    aX = LOG(aX);
    bX = LOG(bX);
    gX = LOG(gX);
    dX = LOG(dX);
    eX = LOG(eX);
    
    aY1M = LOG(1.0f - aY);
    bY1M = LOG(1.0f - bY);
    gY1M = LOG(1.0f - gY);
    dY1M = LOG(1.0f - dY);
    eY1M = LOG(1.0f - eY);
    
    aY = LOG(aY);
    bY = LOG(bY);
    gY = LOG(gY);
    dY = LOG(dY);
    eY = LOG(eY);
    
	//insert states: 1, 2, 11, 12
    temp->i3_1 = aY - eY;
    temp->i1_1 = dY - eY;
    
    temp->i3_2 = aY1M + aX - eX;
    temp->i2_2 = dX - eX;
    
    temp->i6_11 = aX + gY - eX;
    temp->i7_11 = aX + gY + gX1M - eX;
    temp->i11_11 = dX - eX;
    
    temp->i8_12 = gX + aY - eY;
    temp->i7_12 = aY + gX + gY1M - eY;
    temp->i12_12 = dY - eY;
    
    //match
    temp->m4_5 = bX1M + bY1M - eX - eY;
    
    //delete (not double)
    temp->d4_8 = bX + bY1M - eY;
    temp->d10_8 = bY1M - eY;
    
    temp->d4_6 = bX1M + bY - eX;
    temp->d9_6 = bX1M - eX;
    
    //junction-tree conversion 7-0, 10-1, 9-2, 3-3, 4-4
    
    //silent-states
    temp->s5_3 = LOG_ONE;
    temp->s6_3 = gY1M;
    temp->s7_3 = gX1M + gY1M;
    temp->s8_3 = gX1M;
    
    temp->s12_10 = dY1M; 
    temp->s8_10 = aY1M + gX;
    temp->s7_10 = aY1M + gX + gY1M;
    calculateExpProb(aX, aY, bX, bY, gX, gY, aX1M, aY1M, gX1M, gY1M, 1, temp->summedStates_10);
    
    temp->s11_9 = dX1M;
    temp->s6_9 = aX1M + gY;
    temp->s7_9 = aX1M + gX1M + gY;
    calculateExpProb(aX, aY, bX, bY, gX, gY, aX1M, aY1M, gX1M, gY1M, 2, temp->summedStates_9);
    
    temp->s1_4 = dY1M + aX1M;
    temp->s2_4 = dX1M;
    temp->s3_4 = aY1M + aX1M;
    calculateExpProb(aX, aY, bX, bY, gX, gY, aX1M, aY1M, gX1M, gY1M, 4, temp->summedStates_4);
    
    //delete-delete
    temp->s10_7 = bY;
    temp->s9_7 = bX;
    temp->s4_7 = bX + bY;
    temp->s7_7 = gX + gY;
    
    temp->eT = eX1M + eY1M;
    
    //uglyf(" temp->i3_1 " FLOAT_STRING "  temp->i1_1 " FLOAT_STRING "  temp->i3_2 " FLOAT_STRING "  temp->i2_2 " FLOAT_STRING " temp->i6_11 " FLOAT_STRING " temp->i7_11 " FLOAT_STRING " temp->i11_11 " FLOAT_STRING " temp->i8_12 " FLOAT_STRING " temp->i7_12 " FLOAT_STRING " temp->i12_12 " FLOAT_STRING " temp->m4_5 " FLOAT_STRING " temp->d4_8 " FLOAT_STRING " temp->d10_8 " FLOAT_STRING " temp->d4_6 " FLOAT_STRING " temp->d9_6 " FLOAT_STRING " temp->s5_3 " FLOAT_STRING " temp->s6_3 " FLOAT_STRING " temp->s7_3 " FLOAT_STRING " temp->s8_3 " FLOAT_STRING " temp->s12_10 " FLOAT_STRING " temp->s8_10 " FLOAT_STRING " temp->s7_10 " FLOAT_STRING " temp->s11_9 " FLOAT_STRING " temp->s6_9 " FLOAT_STRING " temp->s7_9 " FLOAT_STRING " temp->s1_4 " FLOAT_STRING " temp->s2_4 " FLOAT_STRING " temp->s3_4 " FLOAT_STRING " temp->s10_7 " FLOAT_STRING " temp->s9_7 " FLOAT_STRING " temp->s4_7 " FLOAT_STRING " temp->s7_7 " FLOAT_STRING " temp->eT " FLOAT_STRING "\n",  EXP(temp->i3_1),  EXP(temp->i1_1),  EXP(temp->i3_2),  EXP(temp->i2_2),  EXP(temp->i6_11),  EXP(temp->i7_11),  EXP(temp->i11_11),  EXP(temp->i8_12),  EXP(temp->i7_12),  EXP(temp->i12_12),  EXP(temp->m4_5),  EXP(temp->d4_8),  EXP(temp->d10_8),  EXP(temp->d4_6),  EXP(temp->d9_6),  EXP(temp->s5_3),  EXP(temp->s6_3),  EXP(temp->s7_3),  EXP(temp->s8_3),  EXP(temp->s12_10),  EXP(temp->s8_10),  EXP(temp->s7_10),  EXP(temp->s11_9),  EXP(temp->s6_9),  EXP(temp->s7_9),  EXP(temp->s1_4),  EXP(temp->s2_4),  EXP(temp->s3_4), EXP(temp->s10_7), EXP(temp->s9_7), EXP(temp->s4_7), EXP(temp->s7_7), EXP(temp->eT));
    		
    //exit(1);
    return temp;
}
 
inline void insertXFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge, 
                       void (*assignFn)(struct AlignmentDataStructures *aDS, INT_32, INT_32, FLOAT_32)) {	
	//FLOAT_32 i = edge->edgeScore; //edge->insertBranchCost + 
	assignFn(aDS, S3, S2, model->i3_2 + edge->edgeScore);
	assignFn(aDS, S2, S2, model->i2_2 + edge->edgeScore);
	
    assignFn(aDS, S6, S11, model->i6_11 + edge->edgeScore);
    assignFn(aDS, S7, S11, model->i7_11 + edge->edgeScore);
    assignFn(aDS, S11, S11, model->i11_11 + edge->edgeScore);
}

inline void insertYFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge, 
                       void (*assignFn)(struct AlignmentDataStructures *aDS, INT_32, INT_32, FLOAT_32)) {	
	//FLOAT_32 i = edge->edgeScore; //edge->insertBranchCost + 
	//insert states: 1, 2, 11, 12
    assignFn(aDS, S3, S1, model->i3_1 + edge->edgeScore);
    assignFn(aDS, S1, S1, model->i1_1 + edge->edgeScore);
    
    assignFn(aDS, S8, S12, model->i8_12 + edge->edgeScore);
    assignFn(aDS, S7, S12, model->i7_12 + edge->edgeScore);
    assignFn(aDS, S12, S12, model->i12_12 + edge->edgeScore);
}

inline void matchFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model,
            struct Edge *edgeX, struct Edge *edgeY, void (*assignFn)(struct AlignmentDataStructures *aDS, INT_32, INT_32, FLOAT_32)) {
	//match
	static FLOAT_32 m[ALPHABET_SIZE];
	multiplyWV(edgeX->wV, edgeY->wV, m);
    assignFn(aDS, S4, S5, model->m4_5 + edgeX->edgeScore + edgeY->edgeScore + combineWV(m, aDS->subModelX->insertionDistribution) - edgeX->insertBranchCost - edgeY->insertBranchCost);
    //assignFn(aDS, S4, S5, model->m4_5 + edgeX->edgeScore + edgeY->edgeScore + combineWV(edgeX->wV, edgeY->wV) - edgeX->insertBranchCost - edgeY->insertBranchCost);
}

inline void deleteXFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge, 
	                  void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32)) {
	//delete (not double)
    assignFn(aDS, S4, S8, model->d4_8 + edge->edgeScore);
    assignFn(aDS, S10, S8, model->d10_8 + edge->edgeScore);
}

inline void deleteYFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, struct Edge *edge, 
	                  void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32)) {
	assignFn(aDS, S4, S6, model->d4_6 + edge->edgeScore);
    assignFn(aDS, S9, S6, model->d9_6 + edge->edgeScore);
}
    
inline void silentFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, FLOAT_32 *cell,
					  void (*assignFn)(struct AlignmentDataStructures *aDS, INT_32, INT_32, FLOAT_32)) {   
    FLOAT_32 i, j, k;
    
    //LOG_PLUS_EQUALS(&aDS->toCell[toState], aDS->fromCell[fromState] + edgeScore);
	 
	//silent-states
    assignFn(aDS, S5, S3, model->s5_3);
    assignFn(aDS, S6, S3, model->s6_3);
    assignFn(aDS, S8, S3, model->s8_3);
    
    assignFn(aDS, S12, S10, model->s12_10);
    assignFn(aDS, S8, S10, model->s8_10);
    
    assignFn(aDS, S11, S9, model->s11_9);
    assignFn(aDS, S6, S9, model->s6_9);
    
    assignFn(aDS, S1, S4, model->s1_4);
    assignFn(aDS, S2, S4, model->s2_4); 
    assignFn(aDS, S3, S4, model->s3_4);
    
    //junction-tree conversion 7-0, 10-1, 9-2, 3-3, 4-4
    i = cell[S10]; 
    j = cell[S9];
    k = cell[S4];
     
    cell[S10] = i + model->summedStates_10[1];
    cell[S9] = j + model->summedStates_9[2];
    cell[S4] = k + model->summedStates_4[4];
    
    LOG_PLUS_EQUALS(&cell[S7], i + model->summedStates_10[0]);
    LOG_PLUS_EQUALS(&cell[S9], i + model->summedStates_10[2]);
    LOG_PLUS_EQUALS(&cell[S3], i + model->summedStates_10[3]);
    LOG_PLUS_EQUALS(&cell[S4], i + model->summedStates_10[4]);
    
    LOG_PLUS_EQUALS(&cell[S7], j + model->summedStates_9[0]);
    LOG_PLUS_EQUALS(&cell[S10], j + model->summedStates_9[1]);
    LOG_PLUS_EQUALS(&cell[S3], j + model->summedStates_9[3]);
    LOG_PLUS_EQUALS(&cell[S4], j + model->summedStates_9[4]);
    
    LOG_PLUS_EQUALS(&cell[S7], k + model->summedStates_4[0]);
    LOG_PLUS_EQUALS(&cell[S10], k + model->summedStates_4[1]);
    LOG_PLUS_EQUALS(&cell[S9], k + model->summedStates_4[2]);
    LOG_PLUS_EQUALS(&cell[S3], k + model->summedStates_4[3]);
}

inline void silentFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *aDS, INT_32, INT_32, FLOAT_32)) {   
    //silent-states
    assignFn(aDS, S5, S3, model->s5_3);
    assignFn(aDS, S6, S3, model->s6_3);
    assignFn(aDS, S7, S3, model->s7_3);
    assignFn(aDS, S8, S3, model->s8_3);
    
    assignFn(aDS, S12, S10, model->s12_10);
    assignFn(aDS, S8, S10, model->s8_10);
    assignFn(aDS, S7, S10, model->s7_10);
    
    assignFn(aDS, S11, S9, model->s11_9);
    assignFn(aDS, S6, S9, model->s6_9);
    assignFn(aDS, S7, S9, model->s7_9);
    
    assignFn(aDS, S1, S4, model->s1_4);
    assignFn(aDS, S2, S4, model->s2_4);
    assignFn(aDS, S3, S4, model->s3_4);
}

inline void deleteFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *aDS, INT_32, INT_32, FLOAT_32)) {   
    //junction-tree conversion 7-0, 10-1, 9-2, 3-3, 4-4
	assignFn(aDS, S10, S7, model->s10_7);
    assignFn(aDS, S9, S7, model->s9_7);
    assignFn(aDS, S4, S7, model->s4_7);
    assignFn(aDS, S7, S7, model->s7_7);
}

void destructCombinedTransitionModel(struct CombinedTransitionModel *model) {
    free(model);
}

INT_32 stateNo() {
    return STATE_NO;
}
 
FLOAT_32 *startStates(struct CombinedTransitionModel *model) {
    FLOAT_32 *i;
    INT_32 j;
     
    i = mallocLocal(sizeof(FLOAT_32)*STATE_NO);
    for (j=0;j<STATE_NO; j++)
        i[j] = LOG_ZERO;
    i[S3] = LOG_ONE;
    i[S4] = LOG_ONE;
    //i[INSERT_X_COMB_STATE] = LOG_ONE
    //i[INSERT_Y_COMB_STATE] = LOG_ONE
    //i[INSERT_X2_COMB_STATE] = LOG_ONE
    //i[INSERT_Y2_COMB_STATE] = LOG_ONE
    //i[DELETE_X_COMB_STATE] = LOG_ONE
    //i[DELETE_Y_COMB_STATE] = LOG_ONE
    return i;
}

FLOAT_32 *endStates(struct CombinedTransitionModel *model) {
    FLOAT_32 *i;
    INT_32 j;
    
    i = mallocLocal(sizeof(FLOAT_32)*STATE_NO);
    for (j=0;j<STATE_NO; j++)
        i[j] = LOG_ZERO;
    //i = [fastMaths.LOG_ZERO]*stateNo();
    i[S4] = LOG_ONE;
    i[S7] = LOG_ONE;
    i[S9] = LOG_ONE;
    i[S10] = LOG_ONE;
    /*i[S4] = model->eT; //LOG_ONE;
    i[S7] = model->eT; //LOG_ONE;
    i[S9] = model->eT; //LOG_ONE;
    i[S10] = model->eT; //LOG_ONE;*/
    return i;
} 
