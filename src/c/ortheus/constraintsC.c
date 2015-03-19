#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>  

#include "hashTableC.h"
#include "commonC.h"
#include "fastCMaths.h"

#include "constraintsC.h"

#define CONSTRAINT_UNDECIDED 3 //reserved value
#define CONSTRAINT_BASE_SIZE 1000 
#define CONSTRAINT_HASHTABLE_BASE_SIZE 1000
#define STACK_BASE_SIZE 1000

#define CONSTRAINT_MIN INT_32_MIN
#define CONSTRAINT_MAX INT_32_MAX


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//alignment constraINT_32 data structures 
/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// 
///////////////////////////////////////////////////////// 

struct Constraints *constructConstraintsBackward() {
    struct Constraints *constraints;
    constraints = mallocLocal(sizeof(struct Constraints));
    assert(CONSTRAINT_BASE_SIZE >= 2);
    constraints->xList = mallocLocal(sizeof(INT_32)*CONSTRAINT_BASE_SIZE);
    constraints->yList = mallocLocal(sizeof(INT_32)*CONSTRAINT_BASE_SIZE);
    constraints->constraintsList = mallocLocal(sizeof(INT_32)*CONSTRAINT_BASE_SIZE);
    constraints->maxLength = CONSTRAINT_BASE_SIZE;
    constraints->length = 2;
    
    constraints->xList[0] = CONSTRAINT_MAX;
    constraints->yList[0] = CONSTRAINT_MAX;
    
    constraints->xList[1] = CONSTRAINT_MIN;
    constraints->yList[1] = CONSTRAINT_MIN;
    
    constraints->constraintsList[0] = CONSTRAINT_LESS_THAN;
    constraints->constraintsList[1] = CONSTRAINT_LESS_THAN;
    
    return constraints;
}

void destructConstraints(struct Constraints *constraints) {
    //data-structure for storing set of pairwise constraints
    free(constraints->xList);
    free(constraints->yList);
    free(constraints->constraintsList);
    free(constraints);
}

void appendConstraintBackwards(struct Constraints *constraints, INT_32 x, INT_32 y, INT_32 type) {
    //add prime constraINT_32 to end of list of constraints
    assert(constraints->xList[constraints->length-2] > x);
    assert(constraints->yList[constraints->length-2] > y || (constraints->yList[constraints->length-2] == y &&
                                  constraints->constraintsList[constraints->length-2] == CONSTRAINT_LESS_THAN_OR_EQUAL &&
                                  type == CONSTRAINT_LESS_THAN));
    
    if(constraints->length == constraints->maxLength) {
        void *new;
        INT_32 newSize;
        
        newSize = 2*(constraints->maxLength);
        
        new = memcpy(mallocLocal(sizeof(INT_32)*newSize), constraints->xList, sizeof(INT_32)*(constraints->maxLength));
        free(constraints->xList);
        constraints->xList = new;
        
        new = memcpy(mallocLocal(sizeof(INT_32)*newSize), constraints->yList, sizeof(INT_32)*(constraints->maxLength));
        free(constraints->yList);
        constraints->yList = new;
        
        new = memcpy(mallocLocal(sizeof(INT_32)*newSize), constraints->constraintsList, sizeof(INT_32)*(constraints->maxLength));
        free(constraints->constraintsList);
        constraints->constraintsList = new;
        
        constraints->maxLength = newSize;
    }
    
    constraints->xList[constraints->length-1] = x;
    constraints->yList[constraints->length-1] = y;
    constraints->constraintsList[constraints->length-1] = type;
    
    constraints->xList[constraints->length] = CONSTRAINT_MIN;
    constraints->yList[constraints->length] = CONSTRAINT_MIN;
    constraints->constraintsList[constraints->length++] = CONSTRAINT_LESS_THAN;
    
}

static int bSearchComparatorX(INT_32 *i, INT_32 *j) {
    if (*i > *j) {
        return -1;
    }
    if (*i == *j) {
        return 0;
    }
    if (*i > *(j+1)) {
        return 0;
    }
    return 1;
} 

void getXConstraint(struct Constraints *constraints, INT_32 x, INT_32 *xConstraint, INT_32 *yConstraint, INT_32 *constraintType) {
    INT_32 *i;
    INT_32 j;
    j = constraints->length-1;
    if(x <= constraints->xList[j]) {
        *xConstraint = constraints->xList[j];
        *yConstraint = constraints->yList[j];
        *constraintType = constraints->constraintsList[j];
        return;
    }
    i = bsearch(&x, constraints->xList, j, sizeof(INT_32),
                (int (*)(const void *, const void *))bSearchComparatorX);
    //get prime constraINT_32 for poINT_32 x
    //i = bisect.bisect_left(self.xList, x)
    *xConstraint = *i;
    *yConstraint = constraints->yList[i - constraints->xList];
    *constraintType = constraints->constraintsList[i - constraints->xList];
    //return Constraint(constraints->xList[i], constraints->yList[i], constraints->constraintTypes[i])
}

static int bSearchComparatorY(INT_32 *i, INT_32 *j) {
    if (*i < *j) {
        return 1;
    }
    if (*i == *j) {
        return 0;
    }
    if (*i < *(j-1)) {
        return 0;
    }
    return -1;
}

void getYConstraint(struct Constraints *constraints, INT_32 y, INT_32 *xConstraint, INT_32 *yConstraint, INT_32 *constraintType) {
    //inverse of getXConstraINT_32 
    INT_32 *i;
    
    if(y >= constraints->yList[0]) {
        *xConstraint = constraints->xList[0];
        *yConstraint = constraints->yList[0];
        *constraintType = constraints->constraintsList[0];
        return;
    }
    i = bsearch(&y, constraints->yList+1, constraints->length-1, sizeof(INT_32),
                (int (*)(const void *, const void *))bSearchComparatorY);
    //get prime constraINT_32 for poINT_32 x
    //i = bisect.bisect_left(self.xList, x)
    *xConstraint = constraints->xList[i - constraints->yList];
    *yConstraint = *i;
    *constraintType = constraints->constraintsList[i - constraints->yList];
    //return Constraint(constraints->xList[i], constraints->yList[i], constraints->constraintTypes[i])
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//core constraINT_32 building functions
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
 
struct hashtable **getEmptyConstraints(const INT_32 seqNo, const INT_32 boundaries, struct Chunks *intChunks) {
    INT_32 i;
    struct hashtable **tables;
    
    tables = mallocLocal(sizeof(struct hashtable *)*seqNo*seqNo);
    for(i=0; i<seqNo*seqNo; i++) {
        tables[i] = create_hashtable(CONSTRAINT_HASHTABLE_BASE_SIZE, hashtable_intHashKey, hashtable_intEqualKey, (void (*)(void *))destructInt, (void (*)(void *))destructInt);
    }
    if(boundaries) {
        for(i=0; i<seqNo*seqNo; i++) {
            hashtable_insert(tables[i], constructChunkInt(0, intChunks), constructChunkInt(0, intChunks)); //constructInt(0), constructInt(1));
        }
    }
    //return [ [ { 0:1 } for j in xrange(0, seqNo) ] for i in xrange(0, seqNo) ]
    return tables;
}

void addConstraint(struct hashtable **constraintLists, INT_32 *seqLengths, INT_32 seqNo, INT_32 seqX, INT_32 seqY, INT_32 x, INT_32 y, struct Chunks *intChunks) {
    struct hashtable *constraints;
    INT_32 *i;
    constraints = constraintLists[seqX * seqNo + seqY];
    if(y <= seqLengths[seqY]+1 && x >= 0) {
        //if not constraints.has_key(x):
        if((i = hashtable_search(constraints, &x)) == NULL) {
            //printf(" hello %i %i %i %i \n", seqX, seqY, x, y);
            hashtable_insert(constraints, constructChunkInt(x, intChunks), constructChunkInt(y, intChunks)); //constructInt(x), constructInt(y));
        }
        else {
            if(*i > y) {
            //if(*i < y) {
                //printf(" changing of the constraINT_32 %i %i %i %i % \n", seqX, seqY, x, y, i);
                *i = y;
            }
            //assert(*i <= y);
        }
    }
}

void addLessThanConstraints(INT_32 *list, INT_32 listLength, INT_32 *indices, INT_32 relaxValue, 
                            struct hashtable **constraints, INT_32 *seqLengths, INT_32 seqNo, struct Chunks *intChunks) {
    INT_32 i;
    INT_32 j;
    INT_32 seqX;
    INT_32 seqY;
    INT_32 posX;
    INT_32 posY;
    for(i=0; i < listLength; i++) {
        //for j in xrange(i+1, len(list)):
        for(j=i+1; j<listLength; j++) {
            seqX = list[i];
            seqY = list[j];
            posX = indices[seqX];
            posY = indices[seqY];
            
            //posX = indices[seqX]-relaxValue;
            //posY = indices[seqY]+relaxValue;
            //addConstraint(constraints, seqLengths, seqNo, seqX, seqY, posX, posY+1, intChunks);
            //addConstraint(constraints, seqLengths, seqNo, seqX, seqY, posX-1, posY, intChunks);
            addConstraint(constraints, seqLengths, seqNo, seqX, seqY, posX, posY+1+relaxValue, intChunks);
            addConstraint(constraints, seqLengths, seqNo, seqX, seqY, posX-1-relaxValue, posY, intChunks);
            
            //posX = indices[seqX]+relaxValue;
            //posY = indices[seqY]-relaxValue;
            //addConstraint(constraints, seqLengths, seqNo, seqY, seqX, posY, posX+1, intChunks);
            //addConstraint(constraints, seqLengths, seqNo, seqY, seqX, posY-1, posX, intChunks);
            addConstraint(constraints, seqLengths, seqNo, seqY, seqX, posY, posX+1+relaxValue, intChunks);
            addConstraint(constraints, seqLengths, seqNo, seqY, seqX, posY-1-relaxValue, posX, intChunks);
        }
    }
}

void addMiddleConstraints(INT_32 *nonGaps, INT_32 nonGapsLength, 
                          INT_32 *gaps, INT_32 gapsLength,
                          INT_32 *indices, INT_32 relaxValue, 
                          struct hashtable **constraints, INT_32 *seqLengths, INT_32 seqNo, struct Chunks *intChunks) {
    INT_32 i;
    INT_32 j;
    INT_32 seqX;
    INT_32 seqY;
    INT_32 posX;
    INT_32 posY;
    
    for(i=0; i < nonGapsLength; i++) {
        //for j in xrange(i+1, len(nonGaps)):
        for(j=0; j<gapsLength; j++) {
            seqX = nonGaps[i];
            seqY = gaps[j];
            posX = indices[seqX];
            posY = indices[seqY];
            
            //posX = indices[seqX]-relaxValue;
            //posY = indices[seqY]+relaxValue;
            //addConstraint(constraints, seqLengths, seqNo, seqX, seqY, posX, posY+1, intChunks);
            addConstraint(constraints, seqLengths, seqNo, seqX, seqY, posX, posY+1+relaxValue, intChunks);
            
            //posX = indices[seqX]+relaxValue;
            //posY = indices[seqY]-relaxValue;
            //addConstraint(constraints, seqLengths, seqNo, seqY, seqX, posY, posX, intChunks);
            addConstraint(constraints, seqLengths, seqNo, seqY, seqX, posY-relaxValue, posX, intChunks);
        }
    }
}

struct hashtable **convertAlignmentToInputConstraints(char **alignment, INT_32 alignmentLength, INT_32 seqNo, INT_32 *seqLengths,
                                                     INT_32 relaxValue, char GAP, struct Chunks *intChunks) {
    //reads in a given alignment and converts it to a map of constraints
    struct hashtable **constraints;
    INT_32 *indices;
    INT_32 i;
    INT_32 j;
    INT_32 k;
    INT_32 l;
    INT_32 *list;
    INT_32 *list2;
  
    constraints = getEmptyConstraints(seqNo, TRUE, intChunks);
    indices = mallocLocal(sizeof(INT_32)*seqNo); //[0]*seqNo
    for(i=0; i<seqNo; i++) {
        indices[i] = 0;
    }
    list = mallocLocal(sizeof(INT_32)*seqNo);
    list2 = mallocLocal(sizeof(INT_32)*seqNo);
    //while((column = alignment()) != NULL) {
    for(i=0; i<alignmentLength; i++) {
        j=0;
        k=0;
        for(l=0; l<seqNo; l++) {
            //if column[columnIndex] != substitution.GAP:
            if(alignment[l][i] != GAP) {
                list[j++] = l;
                indices[l]++;
            }
            else {
                list2[k++] = l;
            }
        }  
        addLessThanConstraints(list, j, indices, relaxValue, constraints, seqLengths, seqNo, intChunks);
        addMiddleConstraints(list, j, list2, k, indices, relaxValue, constraints, seqLengths, seqNo, intChunks);
    }
    for(i=0; i<seqNo; i++) {
        assert(indices[i] == seqLengths[i]);
    }
    //assert [ i for i in indices ] == seqLengths
    free(indices);
    free(list);
    free(list2);
    return constraints;
}
            
INT_32 isAligned(INT_32 *nonGaps, INT_32 nonGapsLength, char **alignment, INT_32 columnIndex, char GAP) {
    INT_32 i;
    
    for(i=0; i<nonGapsLength; i++) {
        if(alignment[nonGaps[i]][columnIndex] != GAP) {
            return TRUE;
        }
    }
    return FALSE;
}

INT_32 periodIsZero_Iter(INT_32 *nonGaps, INT_32 nonGapsLength, INT_32 *gaps, INT_32 gapsLength, INT_32 from, INT_32 to, INT_32 change, char **alignment, char GAP, INT_32 *l3) {
    INT_32 i;
    INT_32 j;
    INT_32 lLength;
    INT_32 l2Length;
    INT_32 l3Length;
    
    static INT_32 *l;
    static INT_32 *l2;
    static INT_32 lSize;
    static INT_32 l2Size;
    
    l = arrayResize(l, &lSize, nonGapsLength + gapsLength, sizeof(INT_32));
    l2 = arrayResize(l2, &l2Size, nonGapsLength + gapsLength, sizeof(INT_32));
    lLength = nonGapsLength;
    l2Length = gapsLength;
    l3Length = 0; 
    memcpy(l, nonGaps, nonGapsLength*sizeof(INT_32));
    memcpy(l2, gaps, gapsLength*sizeof(INT_32));
    
    for(i=from; i != to && l2Length != 0; i += change) {
        for(j=0; j<l2Length;) {
            if(alignment[l2[j]][i] != GAP) {
                if(isAligned(l, lLength, alignment, i, GAP)) {
                    l[lLength++] = l2[j];
                }
                else {
                    l3[l3Length++] = l2[j];
                }
                memmove(l2 + j, l2 + j + 1, (l2Size-(j+1))*sizeof(INT_32));
                l2Length--;
            }
            else {
                j++;
            }
        }
    }
    return l3Length;
}

struct BinaryTree *inSingleClade(struct BinaryTree *binaryTree, INT_32 *gapSeqs, INT_32 gapsLength) {
    struct BinaryTree *binarySubTree;
    if(binaryTree == NULL) {
        return NULL;
    }
    if(leftMostLeafNo(binaryTree->traversalID) <= gapSeqs[0] && 
       rightMostLeafNo(binaryTree->traversalID) >= gapSeqs[gapsLength-1]) {
        binarySubTree = inSingleClade(binaryTree->left, gapSeqs, gapsLength);
        if(binarySubTree != NULL) {
            return binarySubTree;
        }
        binarySubTree = inSingleClade(binaryTree->right, gapSeqs, gapsLength);
        if(binarySubTree != NULL) {
            return binarySubTree;
        }
        return binaryTree;
    }
    else {
        return NULL;
    }       
}

INT_32 merge(INT_32 *l, INT_32 *l2, INT_32 lLength, INT_32 l2Length) {
    INT_32 i;
    INT_32 j;
    INT_32 k;
    
    for(i=0; i<l2Length;) {
        j = l2[i];
        for(k=0; k<lLength; k++) {
            if(j <= l[k]) {
                if(j < l[k]) {
                    memmove(l+k+1, l+k, (lLength-k)*sizeof(INT_32));
                    l[k] = j;
                    lLength++;
                }
                goto outer;
            }
        }
        l[k] = j;
        lLength++;
        outer:
        i++;
    }
    return lLength;
}

INT_32 periodIsZero(INT_32 *nonGaps, INT_32 nonGapsLength, INT_32 *gapSeqs, INT_32 gapsLength, INT_32 columnIndex, char **alignment, INT_32 alignmentLength, char GAP, struct BinaryTree *binaryTree) {
    INT_32 lLength;
    INT_32 l2Length;
    
    static INT_32 *l;
    static INT_32 *l2;
    static INT_32 lSize;
    static INT_32 l2Size;
    
    l = arrayResize(l, &lSize, gapsLength, sizeof(INT_32));
    l2 = arrayResize(l2, &l2Size, gapsLength, sizeof(INT_32));
    
    lLength = periodIsZero_Iter(nonGaps, nonGapsLength, gapSeqs, gapsLength, columnIndex-1, -1, -1, alignment, GAP, l);
    l2Length = periodIsZero_Iter(nonGaps, nonGapsLength, gapSeqs, gapsLength, columnIndex+1, alignmentLength, 1, alignment, GAP, l2);
    
    merge(l, l, 0, lLength); //sorts l
    lLength = merge(l, l2, lLength, l2Length);
    if(lLength == 0) {
        return TRUE;
    }
    binaryTree = inSingleClade(binaryTree, l, lLength);
    return leafNoInSubtree(binaryTree->traversalID) == lLength;
}
             
void convertAlignmentToPhylogeneticInputConstraints_Recursion(char **alignment, INT_32 alignmentLength, INT_32 columnIndex, struct BinaryTree *binaryTree,
                                                              INT_32 *indices, INT_32 relaxValue, struct hashtable **constraints, INT_32 *seqLengths, INT_32 seqNo, 
                                                              struct Chunks *intChunks, char GAP) {
    INT_32 i;
    INT_32 j;
    INT_32 k;
    INT_32 lLength;
    INT_32 l2Length;
    
    static INT_32 *l;
    static INT_32 *l2;
    static INT_32 lSize;
    static INT_32 l2Size;
    
    if(binaryTree->internal) {
        i = leftMostLeafNo(binaryTree->traversalID);
        j = rightMostLeafNo(binaryTree->traversalID)+1;
        l = arrayResize(l, &lSize, j-i, sizeof(INT_32));
        l2 = arrayResize(l2, &l2Size, j-i, sizeof(INT_32));
        lLength = 0;
        l2Length = 0; 
        
        for(k=i; k<j; k++) { 
            if (alignment[k][columnIndex] == GAP) { 
                l[lLength++] = k; 
            }
            else { 
                l2[l2Length++] = k; 
            }
        } 
        if(!periodIsZero(l2, l2Length, l, lLength, columnIndex, alignment, alignmentLength, GAP, binaryTree)) {
            convertAlignmentToPhylogeneticInputConstraints_Recursion(alignment, alignmentLength, columnIndex, binaryTree->left, 
                                                             indices, relaxValue, constraints, seqLengths, seqNo, intChunks, GAP);
            convertAlignmentToPhylogeneticInputConstraints_Recursion(alignment, alignmentLength, columnIndex, binaryTree->right, 
                                                             indices, relaxValue, constraints, seqLengths, seqNo, intChunks, GAP);
            return;                                                 
        }
        //success, can add constraints and return
        addLessThanConstraints(l2, l2Length, indices, relaxValue, constraints, seqLengths, seqNo, intChunks);
    }
}

struct hashtable **convertAlignmentToPhylogeneticInputConstraints(char **alignment, INT_32 alignmentLength, INT_32 seqNo, INT_32 *seqLengths,
                                                   INT_32 relaxValue, char GAP, struct Chunks *intChunks, struct BinaryTree *binaryTree) {
    struct hashtable **constraints;
    INT_32 *indices;
    INT_32 i;
    INT_32 j;
    
    //reads in a given alignment and converts it to a map of constraints
    constraints =  getEmptyConstraints(seqNo, TRUE, intChunks);
    indices = callocLocal(seqNo, sizeof(INT_32)); //[0]*seqNo
   
    for(i=0; i<alignmentLength; i++) {
        for(j=0; j<seqNo; j++) {
            if(alignment[j][i] != GAP) {
                indices[j]++;
            }
        } 
        convertAlignmentToPhylogeneticInputConstraints_Recursion(alignment, alignmentLength, i, binaryTree,
                                                                 indices, relaxValue, constraints, seqLengths, seqNo, intChunks, GAP);
    }
    for(i=0; i<seqNo; i++) {
        assert(indices[i] == seqLengths[i]);
    }
    free(indices);
    return constraints;
}

static INT_32 selectedSeq;
static struct hashtable **lessThanConstraints;
static struct hashtable **lessThanOrEqualConstraints;
static INT_32 seqNo;
static INT_32 *seqLengths;

static struct Constraints **primeConstraints;
static INT_32 *prime;
static INT_32 *constraintType;
static INT_32 *list;
static INT_32 listIndex;
static INT_32 **vertices;

static INT_32 *stack;
static INT_32 stackLength;
static INT_32 stackMaxLength;

void searchFrom(INT_32 vertexSeq, INT_32 vertexPos, INT_32 type, struct hashtable **constraintsLists) {
    INT_32 seq;
    struct hashtable *constraints;
    INT_32 *i;
    
    //for seq in xrange(seqNo-1, -1, -1):
    for(seq=seqNo-1; seq>=0; seq--) {
        if(seq != vertexSeq) {
            constraints = constraintsLists[vertexSeq*seqNo + seq];
            //if constraints.has_key(vertexPos):
            if ((i = hashtable_search(constraints, &vertexPos)) != NULL) {
                stack = arrayPrepareAppend(stack, &stackMaxLength, stackLength+3, sizeof(INT_32));
                stack[stackLength++] = type;
                stack[stackLength++] = *i;
                stack[stackLength++] = seq;
                //stack.append((seq, i, type));
            }
        }
    }
}

void search(INT_32 vertexSeq, INT_32 vertexPos, INT_32 type) {
    INT_32 vertexMark;
    
    while(TRUE) {
        vertexMark = vertices[vertexSeq][vertexPos];
        if (vertexMark == CONSTRAINT_UNDECIDED ||
        (vertexMark == CONSTRAINT_LESS_THAN_OR_EQUAL &&
         type == CONSTRAINT_LESS_THAN)) {
            vertices[vertexSeq][vertexPos] = type;
            if (prime[vertexSeq] >= vertexPos) {
                if (constraintType[vertexSeq] == CONSTRAINT_UNDECIDED && vertexSeq != selectedSeq) {
                    list[listIndex++] = vertexSeq;
                }
                prime[vertexSeq] = vertexPos;
                constraintType[vertexSeq] = type;
            }
            searchFrom(vertexSeq, vertexPos, type, lessThanOrEqualConstraints);
            searchFrom(vertexSeq, vertexPos, CONSTRAINT_LESS_THAN, lessThanConstraints);
            if(vertexPos+1 < seqLengths[vertexSeq]) {
                vertexPos++;
                type = CONSTRAINT_LESS_THAN;
                continue;
            }
        }
        if(stackLength == 0) {
            break;
        }
        vertexSeq = stack[--stackLength];
        vertexPos = stack[--stackLength];
        type = stack[--stackLength];
    }
}

void markup(INT_32 vertexSeq, INT_32 vertexPos) {
    struct hashtable *constraints;
    INT_32 vertexMark;
    INT_32 seq;
    INT_32 *i;
    
    while(TRUE) {
        vertexMark = vertices[vertexSeq][vertexPos];
        if(vertexMark != CONSTRAINT_LESS_THAN) {
            vertices[vertexSeq][vertexPos] = CONSTRAINT_LESS_THAN;
            for(seq=seqNo-1; seq >= 0; seq--) {
                constraints = lessThanOrEqualConstraints[vertexSeq * seqNo + seq];
                if((i = hashtable_search(constraints, &vertexPos)) != NULL) { //constraints.has_key(vertexPos)) {
                    //i = constraints[vertexPos];
                    stack = arrayPrepareAppend(stack, &stackMaxLength, stackLength+2, sizeof(INT_32));
                    stack[stackLength++] = *i;
                    stack[stackLength++] = seq;
                    //stack.append((seq, i))
                }
            }
        }
        if(stackLength == 0) {
            break;
        }
        vertexSeq = stack[--stackLength];
        vertexPos = stack[--stackLength];
        //vertexSeq, vertexPos = stack.pop()
    }
}

struct Constraints **buildConstraints(struct hashtable **lessThanConstraintsA, struct hashtable **lessThanOrEqualConstraintsA, INT_32 selectedSeqA, INT_32 seqNoA, INT_32 *seqLengthsA) {
    //computes set of prime constraints from set of input constraints
    //
    //exact copied, non-optimised implementation of Myers and Millers prime constraINT_32 algorithm
    //see page 14 of Progressive Multiple Alignment with Constraints, Myers et al.
    
    INT_32 i;
    INT_32 j;
    INT_32 k;
    INT_32 vertexPos;
    INT_32 seq;
    INT_32 *temp;
    
    selectedSeq = selectedSeqA;
    lessThanConstraints = lessThanConstraintsA;
    lessThanOrEqualConstraints = lessThanOrEqualConstraintsA;
    seqNo = seqNoA;
    seqLengths = seqLengthsA;
   
    primeConstraints = mallocLocal(sizeof(void *)*seqNo); ///[ [] for i in xrange(0, seqNo) ]
    for(i=0; i<seqNo; i++) {
        primeConstraints[i] = constructConstraintsBackward();
    }
   
    prime = mallocLocal(sizeof(INT_32)*seqNo); //[sys.maxint]*seqNo
    for(i=0; i<seqNo; i++) {
        prime[i] = INT_32_MAX;
    }
   
    constraintType = mallocLocal(sizeof(INT_32)*seqNo); //[CONSTRAINT_UNDECIDED]*seqNo
    for(i=0; i<seqNo; i++) {
        constraintType[i] = CONSTRAINT_UNDECIDED;
    }
   
    list = mallocLocal(sizeof(INT_32)*seqNo); //[]
    
    vertices = mallocLocal(sizeof(INT_32 *)*seqNo);
    for(i=0; i<seqNo; i++) {
        k = seqLengths[i];
        temp = mallocLocal(sizeof(INT_32)*k);
        vertices[i] = temp;
        for(j=0; j<k; j++) {
            temp[j] = CONSTRAINT_UNDECIDED;
        }
    }
        
    stackMaxLength = STACK_BASE_SIZE;
    stack = mallocLocal(sizeof(INT_32)*stackMaxLength);
    stackLength = 0;
    //vertices = [ [ CONSTRAINT_UNDECIDED for pos in xrange(0, seqLengths[seq])]
    //            for seq in xrange(0, seqNo) ]
    //for vertexPos in xrange(seqLengths[selectedSeq]-1, -1, -1):
    for(vertexPos=seqLengths[selectedSeq]-1; vertexPos>=0; vertexPos--) {
        listIndex = 0;
        search(selectedSeq, vertexPos, CONSTRAINT_LESS_THAN_OR_EQUAL);
        assert(stackLength == 0);
        //for seq in list:
        for(i=0; i<listIndex; i++) {
            seq = list[i];
            appendConstraintBackwards(primeConstraints[seq], vertexPos, prime[seq], constraintType[seq]);
            constraintType[seq] = CONSTRAINT_UNDECIDED;
        }
        markup(selectedSeq, vertexPos);
        assert(stackLength == 0);
    }
    
    //memory clean up
    free(prime);
    
    free(constraintType);
    
    free(list);
    listIndex = INT_32_MAX;
    
    for(i=0; i<seqNo; i++) {
        free(vertices[i]);
    }
    free(vertices);
    
    free(stack);
    stackLength = INT_32_MAX;
    stackMaxLength = INT_32_MAX;
    
    selectedSeq = INT_32_MAX;
    lessThanConstraints = NULL;
    lessThanOrEqualConstraints = NULL;
    seqNo = INT_32_MAX;
    seqLengths = NULL;
    //end clean up
    
    return primeConstraints;
}

struct Constraints ***buildAllConstraints(struct hashtable **lessThanConstraints, 
                                          struct hashtable **lessThanOrEqualConstraints,
                                          INT_32 seqNo, INT_32 *seqLengths) {
    INT_32 *temp;
    INT_32 i;
    INT_32 j;
    INT_32 k;
    struct Constraints ***primeConstraintsMatrix;
    struct Constraints *primeConstraints;
    
    temp = mallocLocal(sizeof(INT_32)*seqNo);
    for(i=0; i<seqNo; i++) {
        temp[i] = seqLengths[i]+2;
    }
    primeConstraintsMatrix = mallocLocal(sizeof(void *)*seqNo);
    for(i=0; i<seqNo; i++) { 
        primeConstraintsMatrix[i] = buildConstraints(lessThanConstraints, lessThanOrEqualConstraints, i, seqNo, temp);
    }
    //primeConstraints =  [ buildConstraints(lessThanConstraints, lessThanOrEqualConstraints, 
    //                                       selectedSeq, seqNo, seqLengths) for selectedSeq in xrange(0, seqNo) ]
    //memory clean up 
    free(temp);
    //end clean up
    return primeConstraintsMatrix;
}

struct Constraints ***buildAllConstraints_FromAlignment(char **alignment, INT_32 alignmentLength, INT_32 seqNo, INT_32 *seqLengths,
                                                        INT_32 relaxValue, INT_32 GAP, struct BinaryTree *binaryTree, INT_32 totalConstraints) {
    INT_32 i;
    struct hashtable **lessThanConstraints;
    struct hashtable **lessThanOrEqualConstraints;
    struct Constraints ***primeConstraintsMatrix;
    struct Chunks *intChunks;
    
    intChunks = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(INT_32));
    if(alignment != NULL) {
        if(totalConstraints) {
            lessThanConstraints = convertAlignmentToInputConstraints(alignment, alignmentLength, seqNo, seqLengths, relaxValue, GAP, intChunks);
        }
        else {
            lessThanConstraints = convertAlignmentToPhylogeneticInputConstraints(alignment, alignmentLength, seqNo, seqLengths, relaxValue, GAP, intChunks, binaryTree);
        } 
    }
    else {
        lessThanConstraints = getEmptyConstraints(seqNo, TRUE, intChunks);
    }
    lessThanOrEqualConstraints = getEmptyConstraints(seqNo, FALSE, intChunks);
  
    primeConstraintsMatrix = buildAllConstraints(lessThanConstraints, lessThanOrEqualConstraints, seqNo, seqLengths);
 
    //memory clean up
    for(i=0; i<seqNo*seqNo; i++) {
        hashtable_destroy(lessThanConstraints[i], FALSE, FALSE);
        hashtable_destroy(lessThanOrEqualConstraints[i], FALSE, FALSE);
    }
    free(lessThanConstraints);
    free(lessThanOrEqualConstraints);
    destructChunks(intChunks);
    //end memory clean up
    return primeConstraintsMatrix;                                               
}
