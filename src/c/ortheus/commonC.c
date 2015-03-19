#include "fastCMaths.h"
#include "commonC.h"
#include "hashTableC.h"

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

INT_32 LOG_LEVEL;

inline void logInfo(const char *string, ...) {
   if(LOG_LEVEL >= LOGGING_INFO) {
       va_list ap;
       va_start(ap, string);
       vfprintf(stderr, string, ap);
       va_end(ap); 
   }
}

inline void logDebug(const char *string, ...) {
   if(LOG_LEVEL >= LOGGING_DEBUG) {
       va_list ap;
       va_start(ap, string);
       vfprintf(stderr, string, ap);
       va_end(ap); 
   }
}

inline void uglyf(const char *string, ...) {
    va_list ap;
    va_start(ap, string);
    vfprintf(stderr, string, ap);
    //vfprintf(stdout, string, ap);
    va_end(ap); 
}

inline void setLogLevel(INT_32 level) {
    LOG_LEVEL = level;
}

inline void *mallocLocal(INT_32 i) {
    void *j;
    j = malloc(i);
    if(j == 0) {
        fprintf(stderr, "Malloc failed\n");
        exit(1);
    }
    return j;
}

inline struct Chunks *constructChunks(INT_32 chunkSize, INT_32 elementSize) {
    struct Chunks *chunks;
    
    chunks = mallocLocal(sizeof(struct Chunks));
    chunks->chunkList = constructEmptyList(0, free);
    chunks->remaining = 0;
    chunks->chunkSize = chunkSize;
    chunks->elementSize = elementSize;
    return chunks;
}

inline void destructChunks(struct Chunks *chunk) {
    destructList(chunk->chunkList);
    free(chunk);
}

inline void *mallocChunk(struct Chunks *chunk) {
    if(chunk->remaining-- > 0) {
        return (chunk->chunk += chunk->elementSize);
    }
    else {
        chunk->chunk = mallocLocal(chunk->elementSize*chunk->chunkSize);
        listAppend(chunk->chunkList, chunk->chunk);
        chunk->remaining = chunk->chunkSize-1;
        return chunk->chunk;
    }
}  

inline void *callocLocal(INT_32 i, INT_32 j) {
    void *k;
    k = calloc(i, j);
    if(k == 0) {
        fprintf(stderr, "Calloc failed");
        exit(1);
    }
    return k;
}

inline void *arrayResize_NoCheck(void *current, INT_32 *currentSize, INT_32 newSize, INT_32 base) {
    assert(*currentSize <= newSize);
    if(current != NULL) { 
        free(current); 
    }
    *currentSize = newSize; 
    return mallocLocal(base*newSize);
}

inline void *arrayResize(void *current, INT_32 *currentSize, INT_32 newSize, INT_32 base) {
    if(*currentSize < newSize) {
        return arrayResize_NoCheck(current, currentSize, newSize, base);
    }
    return current;
}

inline void listResize(struct List *list, INT_32 newMaxSize) {
    list->list = arrayResize(list->list, &list->maxLength, newMaxSize, sizeof(void *));
}

inline void *arrayCopyResize_NoCheck(void *current, INT_32 *currentSize, INT_32 newSize, INT_32 base) {
    assert(*currentSize <= newSize);
    void *new;
    new = memcpy(mallocLocal(base*newSize), current, base*(*currentSize));
    if(current != NULL) {
        free(current);
    }
    *currentSize = newSize;
    return new;
}

inline void *arrayCopyResize(void *current, INT_32 *currentSize, INT_32 newSize, INT_32 base) {
    if(*currentSize < newSize) {
        return arrayCopyResize_NoCheck(current, currentSize, newSize, base);
    }
    return current;
}

inline void *arrayPrepareAppend(void *current, INT_32 *maxLength, INT_32 newLength, INT_32 base) {
    if(newLength >= *maxLength) { 
        return arrayCopyResize_NoCheck(current, maxLength, (*maxLength)*2 + newLength + SMALL_CHUNK_SIZE, base);
    } 
    return current;
}

inline void listReverse(struct List *list) {
    INT_32 i;
    void *j;
    INT_32 k;
    
    k = list->length-1;
    for(i=0; i<list->length/2; i++) {
        j = list->list[i];
        list->list[i] = list->list[k - i];
        list->list[k - i] = j;
    }
}

inline void listIntersection(struct List *list, struct List *list2, struct List *list3) {
    //currently quadratic time, watch cost closely
    //output list can be the same as the input list
    INT_32 i;
    INT_32 j;
    INT_32 k;
    static void **scratch;
    static INT_32 scratchSize;
    scratch = arrayResize(scratch, &scratchSize, list->length + 1, sizeof(void *));
    k = 0;
    for(i=0; i<list->length; i++) {
        for(j=0; j<list2->length; j++) {
            if(list->list[i] == list2->list[j]) {
                scratch[k++] = list->list[i];
                break;
            }
        }
    }
    list3->length = 0;
    for(i=0; i<k; i++) {
        listAppend(list3, scratch[i]);
    }
}

inline void listCopyResize(struct List *list, INT_32 newMaxSize) {
    list->list = arrayCopyResize(list->list, &list->maxLength, newMaxSize, sizeof(void *));
}

inline struct hashtable *intListToHash(struct List *list, INT_32 *(*getKey)(void *)) {
    INT_32 i;
    void *key;
    struct hashtable *hT;
    
    hT = create_hashtable(list->length, hashtable_intHashKey, hashtable_intEqualKey, NULL, NULL);
    for(i=0; i<list->length; i++) { 
        key = getKey(list->list[i]); 
        hashtable_insert(hT, key, list->list[i]);
    } 
    return hT;
}

inline void swapListFields(struct List *list1, struct List *list2) {
    assert(list1 != list2);
    assert(list1->list == NULL || list1->list != list2->list);
    void **list;
    INT_32 length;
    INT_32 maxLength;
    
    list = list1->list;
    length = list1->length;
    maxLength = list1->maxLength;
    
    list1->list = list2->list;
    list1->length = list2->length;
    list1->maxLength = list2->maxLength;
    
    list2->list = list;
    list2->length = length;
    list2->maxLength = maxLength;
}

inline void copyList(struct List *from, struct List *to) {
    assert(from != to);
    assert(from->list == NULL || from->list != to->list);
    if(from->length > to->maxLength) {
        listResize(to, from->length);
    }
    to->length = from->length;
    to->list = memcpy(to->list, from->list, sizeof(void *)*from->length);
}

inline void listAppend(struct List *list, void *item) {
    if(list->length >= list->maxLength) {
        list->list = arrayCopyResize_NoCheck(list->list, &list->maxLength, list->maxLength*2 + SMALL_CHUNK_SIZE, sizeof(void *));
    }
    list->list[list->length++] = item;
}

inline void *listRemoveFirst(struct List *list) {
    INT_32 i;
    void *j;
    
    j = list->list[0];
    for(i=1; i<list->length; i++) {
        list->list[i-1] = list->list[i];
    }
    list->length--;
    return j;
}

inline void listAppendArray(struct List *list, void **array, INT_32 length) {
    INT_32 i;
    
    if(list->length + length > list->maxLength) {
        list->list = arrayCopyResize_NoCheck(list->list, &list->maxLength, list->maxLength*2 + length + SMALL_CHUNK_SIZE, sizeof(void *));
    }
    for(i=0; i<length; i++) {
        list->list[list->length++] = array[i];
    }
}

inline INT_32 listGetInt(struct List *list, INT_32 index) {
    assert(list != NULL);
    assert(index >= 0);
    assert(index < list->length);
    return *((INT_32 *)list->list[index]);
}

inline FLOAT_32 listGetFloat(struct List *list, INT_32 index) {
    assert(list != NULL);
    assert(index >= 0);
    assert(index < list->length);
    return *((FLOAT_32 *)list->list[index]);
}

//list functions 
inline struct List *copyConstructList(void **list, INT_32 length, void (*destructElement)(void *)) {
    struct List *i;
    INT_32 j;
    
    i = mallocLocal(sizeof(struct List));
    i->length = length; 
    i->maxLength = length;
    j = sizeof(void *)*length;
    i->list = mallocLocal(j);
    memcpy(i->list, list, j);
    i->destructElement = destructElement;
    return i;
}

inline struct List *constructZeroLengthList(INT_32 length, void (*destructElement)(void *)) {
    struct List *l;
    l = constructEmptyList(length, destructElement);
    l->length = 0;
    return l;
}

inline struct List *constructEmptyList(INT_32 length, void (*destructElement)(void *)) {
    struct List *i;

    i = mallocLocal(sizeof(struct List));
    i->length = length;
    i->maxLength = length;
    i->list = mallocLocal(sizeof(void *)*length);
    i->destructElement = destructElement;
    return i;
}

inline void destructList(struct List *list) {
    INT_32 i;
   
    if (list->destructElement != NULL) { 
        for(i=0; i<list->length; i++) { //only free up to known area of list
            list->destructElement(list->list[i]);
        }
    }
    free(list->list);
    free(list);
}

inline INT_32 *constructChunkInt(INT_32 intValue, struct Chunks *chunks) {
    INT_32 *i;
    
    i = mallocChunk(chunks);
    *i = intValue;
    return i;
}

inline LONG_64 *constructChunkLong(LONG_64 longValue, struct Chunks *chunks) {
    LONG_64 *i;
    
    i = mallocChunk(chunks);
    *i = longValue;
    return i;
}

//ints
inline FLOAT_32 *constructFloat(FLOAT_32 i) { 
    FLOAT_32 *j; 
    
    j = mallocLocal(sizeof(FLOAT_32)); 
    *j = i; 
    return j;
}

inline void destructFloat(FLOAT_32 *i) {
    free(i);
}

inline INT_32 *constructInt(INT_32 i) { 
    INT_32 *j; 
    
    j = mallocLocal(sizeof(INT_32)); 
    *j = i; 
    return j;
}

inline void destructInt(INT_32 *i) {
    free(i);
}

//ints
inline LONG_64 *constructLong(LONG_64 i) {
    LONG_64 *j; 
    
    j = mallocLocal(sizeof(LONG_64)); 
    *j = i; 
    return j;
}

inline void destructLong(LONG_64 *i) {
    free(i);
} 

inline UNSIGNED_INT_32 hashtable_intHashKey( void *k ) { 
    return *((INT_32 *)k); 
}

inline INT_32 hashtable_intEqualKey( void *key1, void *key2 ) {
     return *((INT_32 *)key1) == *((INT_32 *)key2);
}

inline UNSIGNED_INT_32 hashtable_longHashKey( void *k ) {
    return *((LONG_64 *)k); 
}

inline INT_32 hashtable_longEqualKey( void *key1, void *key2 ) {
     return *((LONG_64 *)key1) == *((LONG_64 *)key2);
}

inline INT_32 intComparator(INT_32 *i, INT_32 *j) {
    return *i < *j ? -1 : *i > *j ? 1 : 0;
}

inline INT_32 longComparator(LONG_64 *i, LONG_64 *j) {
    return *i < *j ? -1 : *i > *j ? 1 : 0;
}

inline int intComparator_Int(INT_32 *i, INT_32 *j) {
    return *i < *j ? -1 : *i > *j ? 1 : 0;
}

inline int longComparator_Int(LONG_64 *i, LONG_64 *j) {
    return *i < *j ? -1 : *i > *j ? 1 : 0;
}

inline INT_32 intsComparator(INT_32 *ints1, INT_32 *ints2, INT_32 length) {
    INT_32 i;
    INT_32 j;
    INT_32 k;
    
    //assertLocal INT_32s1->length == INT_32s2->length;
    if (ints1 == ints2) {
        return 0;
    }
        
    for (i = 0; i < length; ++i) {
        j = ints1[i];
        k = ints2[i];
        if (j < k) {
            return -1;
        }
        if (j > k) {
            return 1;
        }
    }
    return 0;
}

inline struct TraversalID *constructTraversalID(INT_32 midStart, INT_32 mid, INT_32 midEnd, INT_32 leafNo) {
    struct TraversalID *traversalID;
    
    traversalID = mallocLocal(sizeof(struct TraversalID));
    traversalID->midStart = midStart;
    traversalID->mid = mid;
    traversalID->midEnd = midEnd;
    traversalID->leafNo = leafNo;
    return traversalID;
}

inline void destructTraversalID(struct TraversalID *traversalID) {
    free(traversalID);
}

inline struct BinaryTree *constructBinaryTree(FLOAT_32 distance, INT_32 internal,
                                       struct BinaryTree *left,
                                       struct BinaryTree *right) {
    struct BinaryTree *binaryTree;
  
    binaryTree = mallocLocal(sizeof(struct BinaryTree));
    binaryTree->distance = distance;
    binaryTree->internal = internal;
    binaryTree->left = left;
    binaryTree->right = right;
    return binaryTree;
}

inline void destructBinaryTree(struct BinaryTree *binaryTree) {
    destructTraversalID(binaryTree->traversalID);
    if(binaryTree->left) {
        destructBinaryTree(binaryTree->left);
    }
    if(binaryTree->right) {
        destructBinaryTree(binaryTree->right);
    }
    free(binaryTree);
}

inline INT_32 leftMostLeafNo(struct TraversalID *traversalID) {
    return traversalID->midStart/2;
}

inline INT_32 rightMostLeafNo(struct TraversalID *traversalID) {
    return traversalID->midEnd/2;
}

inline INT_32 leafNoInSubtree(struct TraversalID *traversalID) {
    return rightMostLeafNo(traversalID) - leftMostLeafNo(traversalID) + 1;
}

inline static void binaryTree_depthFirstNumbers_Traverse(struct BinaryTree *binaryTree,
                                          INT_32 *mid, INT_32 *leafNo) {
    INT_32 i;
    INT_32 j;
    
    if(binaryTree->internal) { 
        //_isInternal(binaryTree):
        i = *mid;
        binaryTree_depthFirstNumbers_Traverse(binaryTree->left, mid, leafNo);
        j = (*mid)++;
        binaryTree_depthFirstNumbers_Traverse(binaryTree->right, mid, leafNo);
        binaryTree->traversalID = constructTraversalID(i, j, *mid, INT_32_MAX);
    }
    else {
        i = (*mid)++;
        binaryTree->traversalID = constructTraversalID(i, i, *mid, (*leafNo)++);
    } 
}

inline void binaryTree_depthFirstNumbers(struct BinaryTree *binaryTree) {
    //get pre-order, post-order and mid-order depth first tree numbers
    INT_32 mid = 0;
    INT_32 leafNo = 0;
    
    binaryTree_depthFirstNumbers_Traverse(binaryTree, &mid, &leafNo);
}

inline void printBinaryTree(FILE *file, struct BinaryTree *binaryTree, char **nodeNames) {
    if(binaryTree->internal) {
        printBinaryTree(file, binaryTree->left, nodeNames);
        fprintf(file, "Internal node\t" INT_STRING "\t" INT_STRING "\t" INT_STRING "\t" INT_STRING "\t\t%f\t%s\n", binaryTree->traversalID->midStart, binaryTree->traversalID->mid, binaryTree->traversalID->midEnd, binaryTree->traversalID->leafNo, binaryTree->distance, nodeNames[binaryTree->traversalID->mid]);
        printBinaryTree(file, binaryTree->right, nodeNames);
    }
    else {
        fprintf(file, "Leaf node\t" INT_STRING "\t" INT_STRING "\t" INT_STRING "\t" INT_STRING "\t\t%f\t%s\n", binaryTree->traversalID->midStart, binaryTree->traversalID->mid, binaryTree->traversalID->midEnd, binaryTree->traversalID->leafNo, binaryTree->distance, nodeNames[binaryTree->traversalID->mid]);
    }
}

inline void annotateTree_Fn(struct BinaryTree *bT, void *(*fn)(struct BinaryTree *i), struct List *list) {
    list->list[bT->traversalID->mid] = fn(bT);
    if(bT->internal) {
        annotateTree_Fn(bT->left, fn, list);
        annotateTree_Fn(bT->right, fn, list);
    }
}

inline void annotateTree(struct BinaryTree *bT, void *(*fn)(struct BinaryTree *i), struct List *list) {
    INT_32 i;
    
    list->length = 0;
    for(i=0; i<bT->traversalID->midEnd; i++) {
        listAppend(list, NULL);
    }
    annotateTree_Fn(bT, fn, list);
}

void getBinaryTreeNodesInMidOrder(struct BinaryTree *binaryTree, struct BinaryTree **labels) {
    labels[binaryTree->traversalID->mid] = binaryTree;
    if(binaryTree->internal) {
        getBinaryTreeNodesInMidOrder(binaryTree->left, labels);
        getBinaryTreeNodesInMidOrder(binaryTree->right, labels);
    }
}

FLOAT_32 linOriginRegression(struct List *pointsX, struct List *pointsY) {
    INT_32 i;
    FLOAT_32 j;
    FLOAT_32 k;
    
    j = 0.0;
    k = 0.0;
    assert(pointsX->length == pointsY->length);
    for(i=0; i<pointsX->length; i++) {
        j += *((FLOAT_32 *)pointsX->list[i]);
        k += *((FLOAT_32 *)pointsY->list[i]);
    }
    if(j != 0.0) {
        return k/j;
    }
    return 1.0;
}
