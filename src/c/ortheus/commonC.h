#ifndef COMMONC_H_
#define COMMONC_H_

#include <stdio.h>
#include "fastCMaths.h"
#include "hashTableC.h"

//utils functions

//logging and debugging
#define DEBUG TRUE

#define LOGGING_OFF 0
#define LOGGING_INFO 1
#define LOGGING_DEBUG 2

extern INT_32 LOG_LEVEL;

inline void setLogLevel(INT_32);

inline void logInfo(const char *string, ...);

inline void logDebug(const char *string, ...);

inline void uglyf(const char *string, ...);

//memory
struct Chunks {
    struct List *chunkList;
    void * chunk;
    INT_32 remaining;
    INT_32 chunkSize;
    INT_32 elementSize;
};

inline struct Chunks *constructChunks(INT_32 chunkSize, INT_32 elementSize);

inline void destructChunks(struct Chunks *);

inline void *mallocChunk(struct Chunks *chunk);

inline void *mallocLocal(INT_32 i);

inline void *callocLocal(INT_32 i, INT_32 j);

//general data structures you always need
//lists 
struct List { 
    INT_32 length;
    INT_32 maxLength;
    void **list;
    void (*destructElement)(void *);
};
 
inline void listAppend(struct List *list, void *i);

inline void *listRemoveFirst(struct List *list);

inline void *arrayResize(void *current, INT_32 *currentSize, INT_32 newSize, INT_32 base);

inline void listIntersection(struct List *list, struct List *list2, struct List *list3);

inline void listResize(struct List *list, INT_32 newMaxSize);

inline INT_32 listGetInt(struct List *list, INT_32 index);

inline FLOAT_32 listGetFloat(struct List *list, INT_32 index);

inline void listReverse(struct List *list);

inline void *arrayCopyResize(void *current, INT_32 *currentSize, INT_32 newSize, INT_32 base);

inline void *arrayPrepareAppend(void *current, INT_32 *maxLength, INT_32 currentLength, INT_32 base);

inline void listCopyResize(struct List *list, INT_32 newMaxSize);
 
inline void swapListFields(struct List *list1, struct List *list2);

inline void copyList(struct List *from, struct List *to);

inline struct hashtable *intListToHash(struct List *list, INT_32 *(*getKey)(void *));

//list functions 
inline struct List *copyConstructList(void **list, INT_32 length, void (*destructElement)(void *));

inline struct List *constructZeroLengthList(INT_32 length, void (*destructElement)(void *));

inline struct List *constructEmptyList(INT_32 length, void (*destructElement)(void *));

inline void destructList(struct List *list);

inline void listAppendArray(struct List *list, void **array, INT_32 length);

//ints
inline INT_32 *constructInt(INT_32 i);

inline void destructInt(INT_32 *i);

inline INT_32 *constructChunkInt(INT_32 intValue, struct Chunks *chunks);

inline LONG_64 *constructChunkLong(LONG_64 longValue, struct Chunks *chunks);
//ints
inline LONG_64 *constructLong(LONG_64 i);

inline void destructLong(LONG_64 *i);

inline FLOAT_32 *constructFloat(FLOAT_32 i);

inline void destructFloat(FLOAT_32 *i);

inline UNSIGNED_INT_32 hashtable_intHashKey( void *k );

inline INT_32 hashtable_intEqualKey( void *key1, void *key2 );

inline UNSIGNED_INT_32 hashtable_longHashKey( void *k );

inline INT_32 hashtable_longEqualKey( void *key1, void *key2 );

inline INT_32 intComparator(INT_32 *i, INT_32 *j);

inline INT_32 longComparator(LONG_64 *i, LONG_64 *j);

inline int intComparator_Int(INT_32 *i, INT_32 *j);

inline int longComparator_Int(LONG_64 *i, LONG_64 *j);

inline INT_32 intsComparator(INT_32 *ints1, INT_32 *ints2, INT_32 length);

struct TraversalID {
    //tree traversal numbers, used as nodeIDs for identifying
    //orders in the tree
    //pre == pre order traversal
    //preEnd == max pre index + 1 of node in subtree
    //mid == mid order (in-order) traversal number
    //def __init__(self, pre, preEnd, mid):
    INT_32 midStart;
    INT_32 mid;
    INT_32 midEnd;
    INT_32 leafNo; 
};

inline struct TraversalID *constructTraversalID(INT_32 midStart, INT_32 mid, INT_32 midEnd, INT_32 leafNo);

inline void destructTraversalID(struct TraversalID *traversalID);

struct BinaryTree { 
    FLOAT_32 distance;
    INT_32 internal;
    struct TraversalID *traversalID;
    struct BinaryTree *left;
    struct BinaryTree *right;
};

inline INT_32 leftMostLeafNo(struct TraversalID *traversalID);

inline INT_32 rightMostLeafNo(struct TraversalID *traversalID);

inline INT_32 leafNoInSubtree(struct TraversalID *traversalID);

inline struct BinaryTree *constructBinaryTree(FLOAT_32 distance, INT_32 internal,
                                              struct BinaryTree *left,
                                              struct BinaryTree *right);

inline void destructBinaryTree(struct BinaryTree *binaryTree);

inline void binaryTree_depthFirstNumbers(struct BinaryTree *binaryTree);

inline void printBinaryTree(FILE *file, struct BinaryTree *binaryTree, char **nodeNames);

inline void annotateTree(struct BinaryTree *bT, void *(*fn)(struct BinaryTree *i), struct List *list);

void getBinaryTreeNodesInMidOrder(struct BinaryTree *binaryTree, struct BinaryTree **labels);

FLOAT_32 linOriginRegression(struct List *pointsX, struct List *pointsY);

#endif /*COMMONC_H_*/
