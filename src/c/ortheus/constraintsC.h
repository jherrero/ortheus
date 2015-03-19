#ifndef CONSTRAINTSC_H_
#define CONSTRAINTSC_H_

#include "fastCMaths.h"
#include "commonC.h"

#define CONSTRAINT_LESS_THAN_OR_EQUAL 1
#define CONSTRAINT_LESS_THAN 2

struct Constraints {
    INT_32 *xList;
    INT_32 *yList;
    INT_32 *constraintsList;
    INT_32 length;
    INT_32 maxLength;
};  

void destructConstraints(struct Constraints *constraints);

struct Constraints ***buildAllConstraints_FromAlignment(char **alignment, INT_32 alignmentLength, INT_32 seqNo, INT_32 *seqLengths, INT_32 relaxValue, INT_32 gap, struct BinaryTree *binaryTree, INT_32 totalConstraints);

void getXConstraint(struct Constraints *constraints, INT_32 x, INT_32 *xConstraint, INT_32 *yConstraint, INT_32 *constraintType);

void getYConstraint(struct Constraints *constraints, INT_32 y, INT_32 *xConstraint, INT_32 *yConstraint, INT_32 *constraintType);

#endif /*CONSTRAINTSC_H_*/
