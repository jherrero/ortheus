#ifndef XYZMODELC_H_
#define XYZMODELC_H_

#include "fastCMaths.h"

#define JUNCTION_STATE_NO 5

//stuff for xyzmodel
//combined model
struct CombinedTransitionModel {
    // this model is derived from a combination of two single-affine 
    // duration branch transducers.
	
	//insert states: 1, 2, 11, 12
    FLOAT_32 i3_1;
    FLOAT_32 i1_1;
    
    FLOAT_32 i3_2;
    FLOAT_32 i2_2;
    
    FLOAT_32 i6_11;
    FLOAT_32 i7_11;
    FLOAT_32 i11_11;
    
    FLOAT_32 i8_12;
    FLOAT_32 i7_12;
    FLOAT_32 i12_12;
    
    //match
    FLOAT_32 m4_5;
    
    //delete (not double)
    FLOAT_32 d4_8;
    FLOAT_32 d10_8;
    
    FLOAT_32 d4_6;
    FLOAT_32 d9_6;
    
    //junction-tree conversion 7-0, 10-1, 9-2, 3-3, 4-4
    
    //silent-states
    FLOAT_32 s5_3;
    FLOAT_32 s6_3;
    FLOAT_32 s8_3;
    FLOAT_32 s7_3;
    
    FLOAT_32 s12_10;
    FLOAT_32 s8_10;
    FLOAT_32 s7_10;
    FLOAT_32 summedStates_10[JUNCTION_STATE_NO];
    
    FLOAT_32 s11_9;
    FLOAT_32 s6_9;
    FLOAT_32 s7_9;
    FLOAT_32 summedStates_9[JUNCTION_STATE_NO];
    
    FLOAT_32 s1_4;
    FLOAT_32 s2_4;
    FLOAT_32 s3_4;
    FLOAT_32 summedStates_4[JUNCTION_STATE_NO];
    
    //double-delete
    FLOAT_32 s10_7;
    FLOAT_32 s9_7;
    FLOAT_32 s4_7;
    FLOAT_32 s7_7;
    
    FLOAT_32 eT;
}; 

#endif /*XYZMODELC_H_*/

