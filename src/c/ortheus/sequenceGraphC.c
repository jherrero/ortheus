#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
 
#include "fastCMaths.h"
#include "hashTableC.h"
#include "heapC.h" 
#include "commonC.h" 
#include "substitutionC.h" 
 
#include "xyzModelC.h" 
#include "sequenceGraphC.h"
#include "constraintsC.h" 

#define CELL_TABLE_INITIAL_SIZE_PER_EDGE 10

void calculateConstraints_GreaterThan(INT_32 *leafPositionsX, INT_32 leftMostLeafNoX, INT_32 leftMostLeafNoY, INT_32 leafSeqNoX, INT_32 leafSeqNoY, 
                                     struct Constraints ***allConstraints, INT_32 *constraints) {
    //Like calculateConstraints, except for 'x-is-greater-than-y' constraints.
    INT_32 seqX;
    INT_32 seqY;
    INT_32 x, x2;
    INT_32 dummy; //these values are dummies
    
    for(seqY=0; seqY<leafSeqNoY; seqY++) {
        //x = constraints[seqY];
        getYConstraint(allConstraints[leftMostLeafNoY + seqY][leftMostLeafNoX + 0], leafPositionsX[0], &x, &dummy, &dummy);
        for(seqX=1; seqX<leafSeqNoX; seqX++) { 
            getYConstraint(allConstraints[leftMostLeafNoY + seqY][leftMostLeafNoX + seqX], leafPositionsX[seqX], &x2, &dummy, &dummy);
            if(x2 > x) {
                x = x2;
            } 
        }
        constraints[seqY] = x;
    }
}

void calculateConstraints_LessThan(INT_32 *leafPositionsX, INT_32 leftMostLeafNoX, INT_32 leftMostLeafNoY, INT_32 leafSeqNoX, INT_32 leafSeqNoY, 
                                   struct Constraints ***allConstraints, INT_32 *constraints) {
    //Finds the constraints for a given set of leaf positions (the x coordinates) on the y sequences. 
    //Default function is for 'x-is-less-than-y' constraints  
    //ignores constraINT_32 type
    INT_32 seqX;
    INT_32 seqY;
    INT_32 y, y2;
    INT_32 dummy; //these values are dummies
    
    for(seqY=0; seqY<leafSeqNoY; seqY++) {
        //y = constraints[seqY];
        getXConstraint(allConstraints[leftMostLeafNoX + 0][leftMostLeafNoY + seqY], leafPositionsX[0], &dummy, &y, &dummy);
        for(seqX=1; seqX<leafSeqNoX; seqX++) {
            getXConstraint(allConstraints[leftMostLeafNoX + seqX][leftMostLeafNoY + seqY], leafPositionsX[seqX], &dummy, &y2, &dummy);
            if(y2 < y) {
                y = y2;
            }
        }
        constraints[seqY] = y; 
    }
}

void makeNonRedundantCollection(struct List *collection, INT_32 (*legalCollectionFn)(INT_32 *i, INT_32 *j, INT_32 k), INT_32 length) {
    //Removes redundancy from a set. Note legalFunction is not necessarily ymmetric.
    INT_32 i;
    INT_32 j;
    void *k;
    void *l;
    
    static struct List list;
    
    listResize(&list, collection->length);
    list.length = 0;
    
    //for i in collection:
    for(i=0; i<collection->length; i++) {
        k = collection->list[i];
        //for j in collection:
        for(j=0; j<i; j++) {
            l = collection->list[j];
            if (!legalCollectionFn(k, l, length)) {
                goto outer;
            }
        }
        for(j=i+1; j<collection->length; j++) {
            l = collection->list[j];
            if ((!legalCollectionFn(k, l, length)) && memcmp(k, l, length*sizeof(INT_32))) {
                goto outer;
            }
        }
        list.list[list.length++] = k;
        outer:
        continue;
    }
    collection->length = list.length;
    memcpy(collection->list, list.list, list.length*sizeof(void *));
}

void filterCollection(struct List *collection, struct List *collection2, struct List *collection3, INT_32 (*legalCollectionFn)(INT_32 *i, INT_32 *j, INT_32 k), INT_32 length) {
    //Removes redundancy from a set. Note legalFunction is not necessarily ymmetric.
    INT_32 i;
    INT_32 j;
    void *k;
    
    listResize(collection3, collection->length);
    collection3->length = 0;
   
    //for i in collection:
    for(i=0; i<collection->length; i++) {
        k = collection->list[i];
        //for j in collection:
        for(j=0; j<collection2->length; j++) {
            if (!legalCollectionFn(k, collection2->list[j], length)) {
                goto outer;
            }
        }
        collection3->list[collection3->length++] = k;
        outer:
        continue;
    }
}
    
void mergeNonRedundantCollections(struct List *collection1, struct List *collection2, INT_32 (*legalCollectionFn)(INT_32 *i, INT_32 *j, INT_32 k), INT_32 length) {
    //Merges two sets of sequenceCoordinates, by finding minimal number of explainatory sets.
    static struct List collection3;
    INT_32 i;
    
    filterCollection(collection1, collection2, &collection3, legalCollectionFn, length);
    filterCollection(collection2, &collection3, collection1, legalCollectionFn, length);
    listCopyResize(collection1, collection1->length + collection3.length);
    for(i=0; i<collection3.length; i++) {
        collection1->list[collection1->length++] = collection3.list[i];
    }
}

inline INT_32 legalCollection_GreaterThan(INT_32 *a, INT_32 *b, INT_32 length) {
    INT_32 i;

    for(i=0; i<length; i++) { 
        if (a[i] > b[i]) {
            return TRUE;
        }
    }
    return FALSE;
}

inline INT_32 legalCollection_LessThan(INT_32 *a, INT_32 *b, INT_32 length) {
    INT_32 i;

    for(i=0; i<length; i++) { 
        if (a[i] < b[i]) {  
            return TRUE;
        }
    }
    return FALSE;
}
    
inline INT_32 lessThan(INT_32 i, INT_32 j) {
    return i < j;
}

inline INT_32 greaterThan(INT_32 i, INT_32 j) {
    return i > j;
}

inline INT_32 edgeFrom(struct Edge *edge) {
    return edge->from;
}

inline INT_32 edgeTo(struct Edge *edge) {
    return edge->to;
}

struct List *cloneCollection(struct List *oldList, INT_32 length, struct Chunks *chunks) {
    INT_32 i;
    struct List *newList;
    
    newList = constructEmptyList(oldList->length, NULL);
    for(i=0; i<oldList->length; i++) {
        newList->list[i] = memcpy(mallocChunk(chunks), oldList->list[i], sizeof(INT_32)*length);
    }
    return newList;
}

void getActiveLeaves(INT_32 *activeLeaves, INT_32 leftMostSeqNo, struct TreeNode *treeNode) {
    if(treeNode->type == TREE_NODE_INTERNAL) {
        getActiveLeaves(activeLeaves, leftMostSeqNo, treeNode->treeNodeX);
        getActiveLeaves(activeLeaves, leftMostSeqNo, treeNode->treeNodeY);
    }
    else if(treeNode->type == TREE_NODE_LEAF) {
        activeLeaves[treeNode->traversalID->leafNo - leftMostSeqNo] = TRUE;
    }
}

void calculateMergedSequenceCoordinates(INT_32 **vertexSequenceCoordinates, struct SequenceGraph *sequenceGraph, 
                                        INT_32 leftMostSeqNo, INT_32 leafSeqNo, void ***mergedVertices, void ***mergedEdges, 
                                        INT_32 vertexStart, INT_32 vertexEnd, INT_32 vertexIncrement, INT_32 endIncrement, 
                                        struct List **connections, INT_32 (*legalCollectionFn)(INT_32 *, INT_32 *, INT_32), 
                                        INT_32 (*edgeTo)(struct Edge *), INT_32 (*lessThan)(INT_32 i, INT_32 j), struct Chunks *chunks) {
    //Calculates the set of correct sequenceCoordinates for each vertex and non-silent edge, taking into
    //account the effect of skipping out sequences in the graph contained in insert
    //edges. Multiple sequence coordinate sets for each vertex are possible, due to multiple possible 
    //paths to the right of a given vertex. Function will try to find minimum number possible.
    void **mergedSequenceCoordinates_Vertices;
    void **mergedSequenceCoordinates_Edges;
    INT_32 vertex;
    INT_32 i; 
    INT_32 j;
    INT_32 seq;
    struct List *edges;
    struct Edge *edge;
    struct List *tempList;
    INT_32 *tempInt;
    INT_32 *activeLeaves;
    struct List *sequenceCoordinatesCollection;
    INT_32 *sequenceCoordinatesTo;
    INT_32 *sequenceCoordinatesFrom;
    INT_32 *sequenceCoordinates;
                                     
    static struct List mergeList;
  
    mergedSequenceCoordinates_Vertices = mallocLocal(sizeof(void *)*sequenceGraph->vertexNo);
    mergedSequenceCoordinates_Edges = mallocLocal(sizeof(void *)*sequenceGraph->edges->length);
    activeLeaves = mallocLocal(sizeof(INT_32)*leafSeqNo);
   
    //for vertex in vertexOrder(sequenceGraph):
    for(vertex=vertexStart; vertex!=vertexEnd; vertex+=vertexIncrement) {
        edges = connections[vertex];
        if (edges->length == 0) {
            tempList = constructEmptyList(1, NULL);
            mergedSequenceCoordinates_Vertices[vertex] = tempList;
            tempInt = mallocChunk(chunks);
            tempList->list[0] = tempInt;
            for(i=0; i<leafSeqNo; i++) {
                tempInt[i] = vertexSequenceCoordinates[vertex][i]+endIncrement;
            }
        }
        else {
            //for edge in edges: 
            mergeList.length = 0;
            for(i=0; i<edges->length; i++) {
                edge = edges->list[i];
                if (!edge->silent) {
                    //modifies sequenceCoordinates from to vertices by transitions along edge
                    sequenceCoordinatesCollection = cloneCollection(mergedSequenceCoordinates_Vertices[edgeTo(edge)], leafSeqNo, chunks);
                    sequenceCoordinatesTo = vertexSequenceCoordinates[edge->to];
                    sequenceCoordinatesFrom = vertexSequenceCoordinates[edge->from];
                    //for seq in xrange(0, leafSeqNo):
                    memset(activeLeaves, FALSE, sizeof(INT_32)*leafSeqNo);
                    getActiveLeaves(activeLeaves, leftMostSeqNo, edge->treeNode);
                    for(seq=0; seq<leafSeqNo; seq++) {
                        if (sequenceCoordinatesTo[seq] - sequenceCoordinatesFrom[seq] > 0) {
                            //for sequenceCoordinates in sequenceCoordinatesCollection:
                            for(j=0; j<sequenceCoordinatesCollection->length; j++) {
                                sequenceCoordinates = sequenceCoordinatesCollection->list[j];
                                if (lessThan(sequenceCoordinatesTo[seq], sequenceCoordinates[seq]) && activeLeaves[seq]) { //this is correct, even though you probably think it's not
                                    sequenceCoordinates[seq] = sequenceCoordinatesTo[seq]; //sequenceCoordinatesTo[seq]; 
                                }
                            }
                        }
                    }  
                    makeNonRedundantCollection(sequenceCoordinatesCollection, legalCollectionFn, leafSeqNo);
                    mergedSequenceCoordinates_Edges[edge->iD] = sequenceCoordinatesCollection;
                    mergeNonRedundantCollections(&mergeList, sequenceCoordinatesCollection, legalCollectionFn, leafSeqNo);
                }
                else {
                    mergedSequenceCoordinates_Edges[edge->iD] = mergedSequenceCoordinates_Vertices[edgeTo(edge)];
                    mergeNonRedundantCollections(&mergeList, mergedSequenceCoordinates_Vertices[edgeTo(edge)], legalCollectionFn, leafSeqNo);
                }
            }
            mergedSequenceCoordinates_Vertices[vertex] = copyConstructList(mergeList.list, mergeList.length, NULL);
        }
    }
    free(activeLeaves);
    *mergedVertices =  mergedSequenceCoordinates_Vertices;
    *mergedEdges =  mergedSequenceCoordinates_Edges;
}

void calculateMergedSequenceCoordinates_LeftToRight(INT_32 **vertexSequenceCoordinates, struct SequenceGraph *sequenceGraph, INT_32 leftMostSeqNo, INT_32 leafSeqNo,
                                                   void ***mergedVertices, void ***mergedEdges, struct Chunks *chunks) {
    //As with calculateMergedSequenceCoordinates, but going from left to right 
    //instead of right to left
    calculateMergedSequenceCoordinates(vertexSequenceCoordinates, sequenceGraph, leftMostSeqNo, leafSeqNo, mergedVertices, mergedEdges,
                                       0, sequenceGraph->vertexNo, 1, -1, sequenceGraph->edgesArrangedByToVertex,
                                       legalCollection_LessThan, edgeFrom, greaterThan, chunks);
}


void calculateMergedSequenceCoordinates_RightToLeft(INT_32 **vertexSequenceCoordinates, struct SequenceGraph *sequenceGraph, INT_32 leftMostSeqNo, INT_32 leafSeqNo,
                                                   void ***mergedVertices, void ***mergedEdges, struct Chunks *chunks) {
    //As with calculateMergedSequenceCoordinates, but going from left to right 
    //instead of right to left
    calculateMergedSequenceCoordinates(vertexSequenceCoordinates, sequenceGraph, leftMostSeqNo, leafSeqNo, mergedVertices, mergedEdges,
                                       sequenceGraph->vertexNo-1, -1, -1, 1, sequenceGraph->edgesArrangedByFromVertex,
                                       legalCollection_GreaterThan, edgeTo, lessThan, chunks);
}


void **convertSequenceCoordinatesToConstraints(void **mergedSequenceCoordinates, INT_32 mergedSequenceCoordinateLength,
                                             INT_32 leftMostLeafNoX, INT_32 leftMostLeafNoY, INT_32 leafSeqNoX, INT_32 leafSeqNoY,
                                             struct Constraints ***allConstraints, 
                                             void (*constraintsFn)(INT_32 *, INT_32 , INT_32 , INT_32 , INT_32 , struct Constraints ***, INT_32 *),
                                             INT_32 (*legalCollectionFn)(INT_32 *i, INT_32 *j, INT_32 length), struct Chunks *chunks) {
    //Converts sequence coordinates to constraints on other sequence. Used to convert merged
    //sequence coordinates into constraints
    INT_32 i;
    INT_32 j;
    struct List *collection;
    struct List *collection2;
    INT_32 *temp;
    void **mergedSequenceConstraints;
    
    mergedSequenceConstraints = mallocLocal(sizeof(void *)*mergedSequenceCoordinateLength);
    
    //for key in mergedSequenceCoordinates.keys():
    for(i=0; i<mergedSequenceCoordinateLength; i++) {
       
        collection = mergedSequenceCoordinates[i];
        collection2 = constructEmptyList(collection->length, NULL);
        mergedSequenceConstraints[i] = collection2;
        
        //for sequenceCoordinates in sequenceCoordinatesCollection:
        for(j=0; j<collection->length; j++) {
            temp = mallocChunk(chunks);//mallocLocal(sizeof(INT_32)*leafSeqNoY);
            constraintsFn(collection->list[j], leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY, allConstraints, temp);
            collection2->list[j] = temp;
        }
        makeNonRedundantCollection(collection2, legalCollectionFn, leafSeqNoY);
        //(struct List *collection, INT_32 (*legalCollectionFn)(INT_32 *i, INT_32 *j, INT_32 k), INT_32 length)
    }
    return mergedSequenceConstraints;
}

void **convertSequenceCoordinatesToConstraints_LeftToRight(void **mergedSequenceCoordinates, INT_32 mergedSequenceCoordinateLength,
                                                         INT_32 leftMostLeafNoX, INT_32 leftMostLeafNoY, INT_32 leafSeqNoX, INT_32 leafSeqNoY, 
                                                         struct Constraints ***allConstraints, struct Chunks *chunks) {
    return convertSequenceCoordinatesToConstraints(mergedSequenceCoordinates, mergedSequenceCoordinateLength,
                                                   leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY, allConstraints, 
                                                   calculateConstraints_GreaterThan, legalCollection_LessThan, chunks);                                            
}


void **convertSequenceCoordinatesToConstraints_RightToLeft(void **mergedSequenceCoordinates, INT_32 mergedSequenceCoordinateLength,
                                                         INT_32 leftMostLeafNoX, INT_32 leftMostLeafNoY, INT_32 leafSeqNoX, INT_32 leafSeqNoY,
                                                         struct Constraints ***allConstraints, struct Chunks *chunks) {
    return convertSequenceCoordinatesToConstraints(mergedSequenceCoordinates, mergedSequenceCoordinateLength,
                                                   leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY, allConstraints, 
                                                   calculateConstraints_LessThan, legalCollection_GreaterThan, chunks);                                                               
}

void cleanUpSequenceCoordinates(struct SequenceGraph *sequenceGraphX,
                                void **mergedSequenceCoordinates_Vertices,
                                void **mergedSequenceCoordinates_Edges) {
    INT_32 i;
    //memory clean up
    for(i=0; i<sequenceGraphX->vertexNo; i++) {
        destructList(mergedSequenceCoordinates_Vertices[i]);
    }
    for(i=0; i<sequenceGraphX->edges->length; i++) {
        if(!((struct Edge *)sequenceGraphX->edges->list[i])->silent) {
            destructList(mergedSequenceCoordinates_Edges[i]);
        }
    } 
    free(mergedSequenceCoordinates_Vertices);
    free(mergedSequenceCoordinates_Edges);
}

void calculateMergedConstraints_RightToLeft(INT_32 **vertexXSequenceCoordinates,
                                            INT_32 leftMostLeafNoX, INT_32 leftMostLeafNoY, INT_32 leafSeqNoX, INT_32 leafSeqNoY,
                                            struct Constraints ***allConstraints, struct SequenceGraph *sequenceGraphX,
                                            void ***mergedVertices, void ***mergedEdges, struct Chunks *constraintChunks) {
    void **mergedSequenceCoordinates_Vertices;
    void **mergedSequenceCoordinates_Edges;  
    struct Chunks *chunks;
    
    chunks = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(INT_32)*leafSeqNoX);                                                                           
    //Convenience function to wrap together two functions for calculating the constraints for each vertex
    //and edge in a graph in terms of a given set of sequences present on an opposite graph
    calculateMergedSequenceCoordinates_RightToLeft(vertexXSequenceCoordinates, sequenceGraphX, leftMostLeafNoX, leafSeqNoX,
                                                   &mergedSequenceCoordinates_Vertices, &mergedSequenceCoordinates_Edges, chunks);                                            
    *mergedVertices = convertSequenceCoordinatesToConstraints_RightToLeft(mergedSequenceCoordinates_Vertices, sequenceGraphX->vertexNo, 
                                                                          leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY,
                                                                          allConstraints, constraintChunks);                                              
    *mergedEdges = convertSequenceCoordinatesToConstraints_RightToLeft(mergedSequenceCoordinates_Edges, sequenceGraphX->edges->length, 
                                                                       leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY,
                                                                       allConstraints, constraintChunks);                                                                  
    cleanUpSequenceCoordinates(sequenceGraphX, mergedSequenceCoordinates_Vertices, mergedSequenceCoordinates_Edges); 
    destructChunks(chunks);                                                                                                                                                        
}

void calculateMergedConstraints_LeftToRight(INT_32 **vertexXSequenceCoordinates,
                                            INT_32 leftMostLeafNoX, INT_32 leftMostLeafNoY, INT_32 leafSeqNoX, INT_32 leafSeqNoY,
                                            struct Constraints ***allConstraints, struct SequenceGraph *sequenceGraphX,
                                            void ***mergedVertices, void ***mergedEdges, struct Chunks *constraintChunks) {
    void **mergedSequenceCoordinates_Vertices;
    void **mergedSequenceCoordinates_Edges; 
    struct Chunks *chunks;  
    
    chunks = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(INT_32)*leafSeqNoX);                                                     
    //Convenience function to wrap together two functions for calculating the constraints for each vertex
    //and edge in a graph in terms of a given set of sequences present on an opposite graph
    calculateMergedSequenceCoordinates_LeftToRight(vertexXSequenceCoordinates, sequenceGraphX, leftMostLeafNoX, leafSeqNoX,
                                                   &mergedSequenceCoordinates_Vertices, &mergedSequenceCoordinates_Edges, chunks);
    *mergedVertices = convertSequenceCoordinatesToConstraints_LeftToRight(mergedSequenceCoordinates_Vertices, sequenceGraphX->vertexNo, leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY,
                                                                          allConstraints, constraintChunks);                                                     
    *mergedEdges = convertSequenceCoordinatesToConstraints_LeftToRight(mergedSequenceCoordinates_Edges, sequenceGraphX->edges->length, leftMostLeafNoX, leftMostLeafNoY, leafSeqNoX, leafSeqNoY,
                                                                       allConstraints, constraintChunks);  
    cleanUpSequenceCoordinates(sequenceGraphX, mergedSequenceCoordinates_Vertices, mergedSequenceCoordinates_Edges);      
    destructChunks(chunks);                                                                                                                                                                         
}

INT_32 *calculateSilentVertices(struct SequenceGraph *sequenceGraph) {
    //labels vertices silent if has (only) silent incoming edges
    INT_32 vertex;
    INT_32 *silentVertices;
    INT_32 i;
    struct Edge *edge;
    
    silentVertices = callocLocal(sequenceGraph->vertexNo, sizeof(INT_32));
    //for vertex in xrange(0, sequenceGraph.vertexNo):
    for(vertex=0; vertex<sequenceGraph->vertexNo; vertex++) {
        //for edge in sequenceGraph.edgesArrangedByToVertex[vertex]:
        for(i=0; i<sequenceGraph->edgesArrangedByToVertex[vertex]->length; i++) {
            edge = sequenceGraph->edgesArrangedByToVertex[vertex]->list[i];
            if (edge->silent) {
                silentVertices[vertex] = TRUE;
                break;
            }
        }
    }
    return silentVertices;
}  

//treeNodes
struct TreeNode *copyConstructTreeNode(INT_32 type, INT_32 transitionID, struct TraversalID *traversalID,
                                   struct TreeNode *treeNodeX, struct TreeNode *treeNodeY, FLOAT_32 *wV) {                         
    struct TreeNode *treeNode = mallocLocal(sizeof(struct TreeNode));
    treeNode->left = NULL;
    treeNode->refCount = 0;
    
    treeNode->type = type;
    treeNode->transitionID = transitionID;
    treeNode->traversalID = traversalID;
    treeNode->treeNodeX = treeNodeX;
    treeNode->treeNodeY = treeNodeY;
    treeNode->wV = NULL;
    if (treeNodeX != NULL) {
        treeNodeX->refCount++;
    }
    if (treeNodeY != NULL) {
        treeNodeY->refCount++;
    }
    treeNode->wV = wV ? memcpy(mallocLocal(ALPHABET_SIZE*sizeof(FLOAT_32)), wV, ALPHABET_SIZE*sizeof(FLOAT_32)) : NULL;
    return treeNode;
}

void destructTreeNode(struct TreeNode *treeNode) {
    assert(treeNode->refCount >= 0);
    if(treeNode->refCount == 0) {
        if (treeNode->left != NULL) {
            treeNode->left->refCount--;
            destructTreeNode(treeNode->left);
        }
        if (treeNode->treeNodeX != NULL) {
            treeNode->treeNodeX->refCount--;
            destructTreeNode(treeNode->treeNodeX);
        }
        if (treeNode->treeNodeY != NULL) {
            treeNode->treeNodeY->refCount--;
            destructTreeNode(treeNode->treeNodeY);
        }
        if(treeNode->wV != NULL) {
            free(treeNode->wV);
        }
        free(treeNode);
    }
}

//edges
struct Edge *copyConstructEdge(INT_32 from, INT_32 to, FLOAT_32 edgeScore, FLOAT_32 insertBranchCost, 
                               FLOAT_32 deleteBranchCost, FLOAT_32 *wV, INT_32 silent, void *treeNode, INT_32 iD) {
    struct Edge *temp;
    //INT_32 i;
    temp = mallocLocal(sizeof(struct Edge));
    temp->from = from;
    temp->to = to;
    temp->edgeScore = edgeScore;
    temp->insertBranchCost = insertBranchCost;
    temp->deleteBranchCost = deleteBranchCost;
    //temp->wV = wV;
    temp->wV = wV ? memcpy(mallocLocal(ALPHABET_SIZE*sizeof(FLOAT_32)), wV, ALPHABET_SIZE*sizeof(FLOAT_32)) : NULL;
    temp->silent = silent;
    temp->treeNode = treeNode;
    temp->iD = iD;
    
    return temp;
}

void destructEdge(struct Edge *edge) {
    if(edge->wV)
        free(edge->wV);
    //free(edge->treeNode);
    destructTreeNode(edge->treeNode);
    free(edge);
}

struct TraceBackEdge *constructTraceBackEdge(LONG_64 from, LONG_64 to, FLOAT_32 edgeScore, struct Edge *edgeX, struct Edge *edgeY, char silent, void *getTreeNode) {
//struct TraceBackEdge *constructTraceBackEdge(LONG_64 from, LONG_64 to, FLOAT_32 edgeScore, struct Edge *edgeX, struct Edge *edgeY) {
    struct TraceBackEdge *temp;

    temp = mallocLocal(sizeof(struct TraceBackEdge));
    temp->from = from;
    temp->to = to;
    temp->edgeScore = edgeScore;
    temp->edgeX = edgeX;
    temp->edgeY = edgeY;
    temp->silent = silent;
    temp->getTreeNode = getTreeNode;
    return temp;
}

void destructTraceBackEdge(struct TraceBackEdge *edge) {
    free(edge);
}

inline INT_32 edgeComparator(struct Edge *edge1, struct Edge *edge2) {
    INT_32 i;
    INT_32 j;
    FLOAT_32 k;
    FLOAT_32 l;
    
    i = edge1->to;
    j = edge2->to;
    if (i < j) {
        return -1;
    }
    if (i > j) {
        return 1;
    }
    i = edge1->from;
    j = edge2->from;
    if (i < j) {
        return -1;
    }
    if (i > j) {
        return 1;
    }
    if (edge1 == edge2) {
        return 0;
    }
    k = edge1->edgeScore;
    l = edge2->edgeScore;
    if (k < l) {
        return -1;
    }
    if (k > l) {
        return 1;
    }
    i = edge1->silent;
    j = edge2->silent;
    if (i < j) {
        return -1;
    }
    if (i > j) {
        return 1;
    }
    return edge1->treeNode < edge2->treeNode ? -1 : edge1->treeNode > edge2->treeNode ? 1 : 0;
    //return 0; //a possible but bizarre occurrence!
}

INT_32 edgeComparatorStub(struct AlignmentDataStructures *aDS, struct Edge *edge1, struct Edge *edge2) {
    return edgeComparator(edge1, edge2);
}

int edgePointerComparator(struct Edge **edge1, struct Edge **edge2) {
    return edgeComparator(*edge1, *edge2);
} 

void sortEdges(void **list, INT_32 length) {
    //edgeComparator(struct Edge *edge1, struct Edge *edge2)
    qsort(list, length, sizeof(void *), (int (*)(const void *, const void*))edgePointerComparator);
}

//sequence graphs
struct SequenceGraph *constructSequenceGraph(struct List *edges, INT_32 vertexNo) {
    INT_32 i;
    struct Edge *edge;
    struct SequenceGraph *sequenceGraph;
    INT_32 *counts;
    
    sequenceGraph = mallocLocal(sizeof(struct SequenceGraph));
    sequenceGraph->edges = edges;
    sequenceGraph->vertexNo = vertexNo;
    sequenceGraph->edgesArrangedByToVertex = mallocLocal(sizeof(struct List *)*vertexNo);
    sequenceGraph->edgesArrangedByFromVertex = mallocLocal(sizeof(struct List *)*vertexNo);
    counts = (INT_32 *)mallocLocal(sizeof(INT_32)*vertexNo);
    //tos
    for(i=0; i<vertexNo; i++) {
        counts[i] = 0;
    }
    //for edge in edges:
    for(i=0; i<edges->length; i++) {
        edge = edges->list[i];
        counts[edge->to]++;
    }
    //for edge in edges:
    for(i=0; i<vertexNo; i++) {
        sequenceGraph->edgesArrangedByToVertex[i] = constructEmptyList(counts[i], NULL); //assumes that will be destructed in edges list
        counts[i] = 0;
    }
    for(i=0; i<edges->length; i++) {
        edge = edges->list[i];
        sequenceGraph->edgesArrangedByToVertex[edge->to]->list[counts[edge->to]++] = edge;
    }
    //froms
    for(i=0; i<vertexNo; i++) {
        counts[i] = 0;
    }
    //for edge in edges:
    for(i=0; i<edges->length; i++) {
        edge = edges->list[i];
        counts[edge->from]++;
    }
    //for edge in edges:
    for(i=0; i<vertexNo; i++) {
        sequenceGraph->edgesArrangedByFromVertex[i] = constructEmptyList(counts[i], NULL);
        counts[i] = 0;
    }
    for(i=0; i<edges->length; i++) { 
        edge = edges->list[i];
        sequenceGraph->edgesArrangedByFromVertex[edge->from]->list[counts[edge->from]++] = edge;
    }
    //for(i=0; i<vertexNo; i++) { not necessary
    //    sortEdges(sequenceGraph->edgesArrangedByFromVertex[i]->list, sequenceGraph->edgesArrangedByFromVertex[i]->length);
    //    sortEdges(sequenceGraph->edgesArrangedByToVertex[i]->list, sequenceGraph->edgesArrangedByToVertex[i]->length);
    //}
    free(counts);
    return sequenceGraph;
}

void destructSequenceGraph(struct SequenceGraph *sequenceGraph, INT_32 freeIncomingEdges) {
    INT_32 i;
   
    if(freeIncomingEdges) {
        destructList(sequenceGraph->edges);
    }
    for (i = 0; i < sequenceGraph->vertexNo; ++i) {
        destructList(sequenceGraph->edgesArrangedByFromVertex[i]);
        destructList(sequenceGraph->edgesArrangedByToVertex[i]);
    }
    free(sequenceGraph->edgesArrangedByFromVertex);
    free(sequenceGraph->edgesArrangedByToVertex);
    free(sequenceGraph);
}

/*INT_32 *randomChoices(FLOAT_32 *probs, INT_32 sizeA, INT_32 pathWeight) {
    //method for choosing a number of items according to their given probs
    static INT_32 *choices;
    static INT_32 choicesSize;
    FLOAT_32 maxProb; 
    FLOAT_32 totalProb;
    FLOAT_32 cumulativeProbs[sizeA];
    INT_32 i, j, k;
    
    choices = arrayResize(choices, &choicesSize, sizeA, sizeof(INT_32)); 
    totalProb = probs[0]; 
    cumulativeProbs[0] = totalProb;
    choices[0] = 0;
    for(i=1; i<sizeA; i++) {
    	choices[i] = 0;
        LOG_PLUS_EQUALS(&totalProb, probs[i]);
        //cumulativeProbs[i] = totalProb;
    }
        
    assert(totalProb != LOG_ZERO);
    
    for(i=0; i<sizeA; i++) {
    	probs[i] -= totalProb;
    	probs[i] += LOG(pathWeight);
        //LOG_PLUS_EQUALS(&totalProb, probs[i]);
        //cumulativeProbs[i] = totalProb;
    	//printf(" scores %f %i \n", EXP(probs[i]), pathWeight);
    }
    
    //choiceProbs = [ LOG(RANDOM()) + totalProb for i in xrange(0, pathWeight) ]
    for(i=0; i< pathWeight; ++i) {
    	k = 0;
    	maxProb = LOG_ZERO;
    	for (j = 0; j < sizeA; ++j) {
			if (probs[j] >= maxProb) {
				k = j;
				maxProb = probs[j];
			}
		}
    	choices[k] += 1;
    	probs[k] = LOG(EXP(probs[k]) -1.0);
    }
    for(i=0; i<sizeA; i++) {
    	//printf(" hey %i %i \n",  choices[i], pathWeight);
    }
    return choices;
}*/

INT_32 *randomChoices(FLOAT_32 *probs, INT_32 sizeA, INT_32 pathWeight) {
    //method for choosing a number of items according to their given probs
    static INT_32 *choices;
    static INT_32 choicesSize;
    FLOAT_32 randomChoice; 
    FLOAT_32 totalProb;
    FLOAT_32 cumulativeProbs[sizeA];
    INT_32 i, j;
    
    choices = arrayResize(choices, &choicesSize, sizeA, sizeof(INT_32)); 
    totalProb = probs[0]; 
    cumulativeProbs[0] = totalProb;
    choices[0] = 0;
    for(i=1; i< sizeA; i++) {
    	choices[i] = 0;
        LOG_PLUS_EQUALS(&totalProb, probs[i]);
        cumulativeProbs[i] = totalProb;
    }
        
    assert(totalProb != LOG_ZERO);
    
    //choiceProbs = [ LOG(RANDOM()) + totalProb for i in xrange(0, pathWeight) ]
    for(i=0; i< pathWeight; ++i) {
    	randomChoice = RANDOM_LOG() + totalProb;
    	for (j = 0; j < sizeA; ++j) {
			if (cumulativeProbs[j] >= randomChoice) {
				choices[j]++;
				break;
			} 
		}       
    }
    
    return choices;
}

//graph member holders

struct GraphMemberHolder *constructGraphMember(void *graphMember, INT_32 *sequenceConstraints, void (*destructGraphMember)(void *)) {
	struct GraphMemberHolder *graphMemberHolder;
    graphMemberHolder = mallocLocal(sizeof(struct GraphMemberHolder));
	graphMemberHolder->graphMember = graphMember;
	graphMemberHolder->sequenceConstraints = sequenceConstraints;
    graphMemberHolder->destructGraphMember = destructGraphMember;
	return graphMemberHolder;
}

void destructGraphMember(struct GraphMemberHolder *graphMemberHolder) {
    //don't destruct sequence constraints, as already handled
    if(graphMemberHolder->destructGraphMember != NULL) {
        graphMemberHolder->destructGraphMember(graphMemberHolder->graphMember);
    }
    free(graphMemberHolder);
}

INT_32 getIDFromGraphMemberEdge(void *graphMember) {
    return ((struct Edge *)graphMember)->iD;
}

INT_32 getIDfromGraphMemberVertex(void *graphMember) {
    return *((INT_32 *)graphMember);
}

void computeUnionScratch(struct AlignmentDataStructures *aDS, struct List *list1, struct List *list2, 
                        struct List *scratchList,
                        INT_32 (*compare)(struct AlignmentDataStructures *, void *, void *)) {
    INT_32 i = 0;
    INT_32 j = 0;
    INT_32 k = 0;
    INT_32 l = 0;
    void *graphMember1;
    void *graphMember2; 
    
    INT_32 length1 = list1->length;
    INT_32 length2 = list2->length;
    
    void **array1 = list1->list;
    void **array2 = list2->list;
    void **list3;
    
    listResize(scratchList, length1 + length2);
    list3 = scratchList->list;
    
    assert(list1->destructElement == list2->destructElement);
    while (i < length1) {
        if (j < length2) {
            graphMember1 = array1[i];
            graphMember2 = array2[j];
            k = (*compare)(aDS, graphMember1, graphMember2);
            if (k < 0) {
                list3[l] = graphMember1;
                i += 1;
            }
            else if (k > 0) {
                list3[l] = graphMember2;
                j += 1;
            }
            else {
                list3[l] = graphMember1;
                i += 1;
                j += 1;
            }
            l += 1;
        }
        else {
            k = (length1-i);
            memcpy(list3 + l, array1 + i, k*sizeof(void *));
            scratchList->length = l + k;
            return;
        }
    }
    k = (length2-j);
    memcpy(list3 + l, array2 + j, k*sizeof(void *));
    scratchList->length = l + k;       
}

struct List *computeUnion(struct AlignmentDataStructures *aDS, struct List *list1, struct List *list2, 
                   INT_32 (*compare)(struct AlignmentDataStructures *, void *, void *)) {
    static struct List scratchList;
    
    computeUnionScratch(aDS, list1, list2, &scratchList, compare);
    return copyConstructList(scratchList.list, scratchList.length, list1->destructElement);
}

void diffScratch(struct AlignmentDataStructures *aDS, struct List *list1, struct List *list2, 
                struct List *scratchList,
                INT_32 (*compare)(struct AlignmentDataStructures *, void *, void *)) {
    INT_32 i = 0;
    INT_32 j = 0;
    INT_32 k = 0;
    INT_32 l = 0;
    void *graphMember1;
    void *graphMember2; 
    
    static void **list3;
    
    listResize(scratchList, list1->length);
    list3 = scratchList->list;

    assert(list1->destructElement == list2->destructElement);
    while (i < list1->length) {
        graphMember1 = list1->list[i];
        if (j < list2->length) {
            graphMember2 = list2->list[j];
            k = (*compare)(aDS, graphMember1, graphMember2);
            if (k < 0) { 
                list3[l] = graphMember1;
                l += 1;
                i += 1;
            }
            else if (k > 0) {
                j += 1; 
            }
            else {
                i += 1;
                while (i < list1->length) {
                    graphMember1 = list1->list[i];
                    k = (*compare)(aDS, graphMember1, graphMember2);
                    if (k > 0) {
                        break;
                    }
                    assert(k == 0);
                    i += 1;
                }
                j += 1;
            }
        }
        else {
            k = (list1->length - i);
            memcpy(list3 + l, list1->list + i, k*sizeof(void *));
            scratchList->length = l + k;
            return;
        }
    }
    scratchList->length = l;   
}

void filterScratch(struct AlignmentDataStructures *aDS, struct List *list, 
                  struct List *scratchList,
                  INT_32 (*filterFn)(struct AlignmentDataStructures *, void *)) {
    INT_32 i = 0; 
    INT_32 j = 0;
    
    static void **list2;
    
    void **array = list->list;
    INT_32 length = list->length;
    
    listResize(scratchList, length);
    list2 = scratchList->list;
  
    while (i < length) {
        if ((*filterFn)(aDS, array[i])) {
            list2[j] = array[i]; 
            j += 1;
        }
        i += 1;
    } 
    scratchList->length = j;
}                  

void divideScratch(struct AlignmentDataStructures *aDS, struct List *list, 
                   struct List *scratchList1, struct List *scratchList2,
                   INT_32 (*filterFn)(struct AlignmentDataStructures *, void *)) {
    INT_32 i = 0;
    INT_32 j = 0;
    INT_32 k = 0;
    
    static void **list2;
    static void **list3;
    
    listResize(scratchList1, list->length);
    list2 = scratchList1->list;
    
    listResize(scratchList2, list->length);
    list3 = scratchList2->list;

    while (i < list->length) { 
        if ((*filterFn)(aDS, list->list[i])) {
            list2[j++] = list->list[i];
        }
        else { 
            list3[k++] = list->list[i];
        }
        i += 1;
    }
    scratchList1->length = j;
    scratchList2->length = k;
} 

INT_32 graphMember_EdgeComparator(struct AlignmentDataStructures *aDS, struct GraphMemberHolder *graphMemberHolder1, 
                               struct GraphMemberHolder *graphMemberHolder2) {
    INT_32 i;
    
    i = edgeComparatorStub(aDS, graphMemberHolder1->graphMember, graphMemberHolder2->graphMember);
    if (i == 0) {
        return intsComparator(graphMemberHolder1->sequenceConstraints, 
                              graphMemberHolder2->sequenceConstraints, aDS->leafSeqNoX);
    }
    return i;
}

INT_32 graphMember_VertexComparator(struct AlignmentDataStructures *aDS, struct GraphMemberHolder *graphMemberHolder1, 
                                 struct GraphMemberHolder *graphMemberHolder2) {
    INT_32 i;
    
    i = intComparator(graphMemberHolder1->graphMember, graphMemberHolder2->graphMember);
    if (i == 0) {
        return intsComparator(graphMemberHolder1->sequenceConstraints, graphMemberHolder2->sequenceConstraints, aDS->leafSeqNoX);
    }
    return i;
}
                                 
INT_32 edge_graphMember_EdgeComparator(struct AlignmentDataStructures *aDS, struct Edge *edge, 
                                    struct GraphMemberHolder *graphMemberHolder) {
    return edgeComparatorStub(aDS, edge, graphMemberHolder->graphMember);
}

INT_32 edgeVertexFilter(struct AlignmentDataStructures *aDS, struct Edge *edge) {
    INT_32 *i;
     
    i = bsearch(&(edge->from), aDS->newVertices, aDS->noOfNewVertices, sizeof(INT_32),
                (int (*)(const void *, const void *))intComparator_Int);
    return i != NULL;
} 

INT_32 greaterThanConditionFilter(struct AlignmentDataStructures *aDS,
                               struct GraphMemberHolder *graphMemberHolder) {
    INT_32 i;
    INT_32 seq;
    INT_32 seqCoordinate;
    
    //for seq, seqCoordinate in aDS.changes 
    for (i=0; i<aDS->noOfChanges; i+=2) {
        seq = aDS->changes[i];
        seqCoordinate = aDS->changes[i+1];
        if (seqCoordinate >= graphMemberHolder->sequenceConstraints[seq]) {
            return FALSE;
        }
    }
    return TRUE;
}

INT_32 lessThanConditionFilter(struct AlignmentDataStructures *aDS,
                            struct GraphMemberHolder *graphMemberHolder) {                         
    INT_32 i;
    INT_32 j;
    INT_32 k;
    struct List *sequenceConstraintsCollection;
    INT_32 *sequenceConstraints;
    struct List *sequenceCoordinatesCollection;
    INT_32 *sequenceCoordinates;
    
    sequenceConstraintsCollection = aDS->mergedEndConstraints_GraphMembers[aDS->getIDFromGraphMember(graphMemberHolder->graphMember)];
    sequenceCoordinatesCollection = aDS->sequenceCoordinatesCollection;
    //for sequenceConstraints in aDS.mergedEndConstraints_GraphMembers[graphMember.graphMember]:
    for (i=0; i<sequenceConstraintsCollection->length; i++) {
        sequenceConstraints = sequenceConstraintsCollection->list[i];
        //for sequenceCoordinates in aDS.sequenceCoordinatesCollection 
        j=0;
        outer:
        while (j<sequenceCoordinatesCollection->length) {
            sequenceCoordinates = sequenceCoordinatesCollection->list[j];
            //for seq in xrange(0, aDS.leafSeqNoX):
            for (k=0; k<aDS->leafSeqNoX; k++) {
                if (sequenceCoordinates[k] <= sequenceConstraints[k]) {
                    j++;
                    goto outer;
                }
            }
            return TRUE;
        }
    }
    return FALSE;
}

void deleteEndsScratch(struct AlignmentDataStructures *aDS, struct List *previousVertices, struct List *previousEdges, 
                struct List *previousIllegalEdges, struct List *sequenceCoordinatesCollection,
                struct List *previousSequenceCoordinatesCollection,
                struct List *newVertices, struct List *newEdges, struct List *newIllegalEdges) {            
    INT_32 i;
    static struct List scratchList1; 
    //Deletes edges that were previously legal, but are not longer.
    //if (previousSequenceCoordinatesCollection == sequenceCoordinatesCollection)
    //    return previousVertices, previousEdges, previousIllegalEdges
    if(sequenceCoordinatesCollection->length == previousSequenceCoordinatesCollection->length) {
        for (i=0; i<sequenceCoordinatesCollection->length; i++) {
            if(memcmp(sequenceCoordinatesCollection->list[i], previousSequenceCoordinatesCollection->list[i], aDS->leafSeqNoX*sizeof(INT_32))) {
                goto outer;
            }
        }
        copyList(previousVertices, newVertices); 
        copyList(previousEdges, newEdges); 
        copyList(previousIllegalEdges, newIllegalEdges); 
        return;
    } 
    outer:
    aDS->sequenceCoordinatesCollection = sequenceCoordinatesCollection;
    aDS->mergedEndConstraints_GraphMembers = aDS->mergedEndConstraints_Edges;
    aDS->getIDFromGraphMember = getIDFromGraphMemberEdge; 
    
    divideScratch(aDS, previousEdges, newEdges, newIllegalEdges, (INT_32 (*)(struct AlignmentDataStructures *, void *))lessThanConditionFilter);
    if(newIllegalEdges->length == 0) {
        copyList(previousVertices, newVertices);
        copyList(previousIllegalEdges, newIllegalEdges);
        return;
    } 
    aDS->mergedEndConstraints_GraphMembers = aDS->mergedEndConstraints_Vertices;
    aDS->getIDFromGraphMember = getIDfromGraphMemberVertex;
    
    filterScratch(aDS, previousVertices, newVertices, (INT_32 (*)(struct AlignmentDataStructures *, void *))lessThanConditionFilter);
    //dealing with illegal edges
    //newIllegalEdges = [ i.graphMember for i in newIllegalEdges ]
    for (i = 0; i < newIllegalEdges->length; ++i) {
        newIllegalEdges->list[i] = ((struct GraphMemberHolder *)newIllegalEdges->list[i])->graphMember;
    }
    computeUnionScratch(aDS, newIllegalEdges, previousIllegalEdges, &scratchList1,
                       (INT_32 (*)(struct AlignmentDataStructures *, void *, void *))edgeComparatorStub);
    //aDS.newVertices = [ i.graphMember for i in newVertices ] 
    aDS->newVertices = arrayResize(aDS->newVertices, &aDS->newVertices_Size, newVertices->length, sizeof(INT_32)); 
    for (i = 0; i < newVertices->length; ++i) {
        aDS->newVertices[i] = *((INT_32 *)((struct GraphMemberHolder *)newVertices->list[i])->graphMember);
    } 
    aDS->noOfNewVertices = newVertices->length;
    filterScratch(aDS, &scratchList1, newIllegalEdges, (INT_32 (*)(struct AlignmentDataStructures *, void *))edgeVertexFilter); 
}

struct Edge *getMinEdge(void** edges, INT_32 *size) { //not the best way of doing this, but still
    INT_32 i;
    INT_32 j;
    INT_32 k;
    struct Edge *edge1;
    struct Edge *edge2;
    
    edge1 = edges[0];
    k = 0;
    for(i=1; i<(*size); i++) {
        edge2 = edges[i]; 
        j = edgeComparator(edge1, edge2); 
        if(j > 0) {
            k = i;
            edge1 = edge2;
        }
    }
    edges[k] = edges[--(*size)];
    return edge1;
}

inline INT_32 extendEnds_isLegalEnds(struct AlignmentDataStructures *aDS, INT_32 *sequenceCoordinates, INT_32 *sequenceConstraints) {
    INT_32 seq;
    for (seq=0; seq<aDS->leafSeqNoX; seq++) { 
        if (sequenceCoordinates[seq] <= sequenceConstraints[seq]) {
            return FALSE;
        }
    } 
    return TRUE;
}

void extendEnds(struct AlignmentDataStructures *aDS, struct List *vertices, struct List *edges, struct List *illegalEdges, 
                struct List *sequenceCoordinatesCollection_Ends,
                struct List **newVertices, struct List **newEdges, struct List **newIllegalEdges) {
    //Extends the legal vertices along the y-sequence.
    INT_32 i; 
    INT_32 j;
    INT_32 k;
    INT_32 l;
    INT_32 hashDummy = 0;
    
    struct List temp;
    
    struct hashtable *newVertices_Set;
    
    static void **newVertices_List;
    static INT_32 newVertices_ListSize;
    INT_32 newVertices_List_Index = 0;
    
    static void **newEdges_List;
    static INT_32 newEdges_ListSize;
    INT_32 newEdges_List_Index = 0;
    
    static void **newIllegalEdges_List;
    static INT_32 newIllegalEdges_ListSize;
    INT_32 newIllegalEdges_List_Index = 0;
    
    static void **illegalEdges_List;
    static INT_32 illegalEdges_ListSize;
    INT_32 illegalEdges_List_Index = 0;
    
    void *previous;
    struct Edge *edge;
    
    struct List *sequenceConstraintsCollection;
    INT_32 *sequenceConstraints; 
    
    struct List *startSequenceConstraintsCollection;
    //set up vertex hash 
    newVertices_Set = create_hashtable(vertices->length + SMALL_CHUNK_SIZE, hashtable_intHashKey, hashtable_intEqualKey, NULL, NULL);   
    for (i = 0; i < vertices->length; ++i) { 
        hashtable_insert(newVertices_Set, ((struct GraphMemberHolder *)vertices->list[i])->graphMember, &hashDummy); //value is dummy
    }
    
    //set up illegal edges
    illegalEdges_List = arrayResize(illegalEdges_List, &illegalEdges_ListSize, illegalEdges->length, sizeof(void *));
    for (i = 0; i < illegalEdges->length; ++i) {
        illegalEdges_List[i] = illegalEdges->list[i];
    } 
    illegalEdges_List_Index = illegalEdges->length; 
    previous = NULL; 
    outer:
    while (illegalEdges_List_Index > 0) { 
        edge = getMinEdge(illegalEdges_List, &illegalEdges_List_Index);
        if (edge == previous || hashtable_search(newVertices_Set, &edge->from) == NULL) { //the 
            continue;
        } 
        previous = edge; 
        //for edgeConstraints in aDS.mergedEndConstraints_Edges[edge]:
        sequenceConstraintsCollection = aDS->mergedEndConstraints_Edges[edge->iD];
        for (i=0; i<sequenceConstraintsCollection->length; i++) { 
            sequenceConstraints = sequenceConstraintsCollection->list[i];
            //for sequenceCoordinates in sequenceCoordinatesCollection_Ends:
            for (j=0; j<sequenceCoordinatesCollection_Ends->length; j++) { 
                if (extendEnds_isLegalEnds(aDS, sequenceCoordinatesCollection_Ends->list[j], sequenceConstraints)) { 
                    //is legal 
                    //for startSequenceConstraints in aDS.mergedStartConstraints_Edges[edge]:
                    startSequenceConstraintsCollection = aDS->mergedStartConstraints_Edges[edge->iD]; 
                    newEdges_List = arrayPrepareAppend(newEdges_List, &newEdges_ListSize, newEdges_List_Index + startSequenceConstraintsCollection->length, sizeof(void *)); 
                    l = newEdges_List_Index;
                    for (k=0; k<startSequenceConstraintsCollection->length; k++) { 
                        newEdges_List[newEdges_List_Index++] = startSequenceConstraintsCollection->list[k];
                    }
                    if(newEdges_List_Index > l) { 
                        if (hashtable_search(newVertices_Set, &edge->to) == NULL) { 
                            hashtable_insert(newVertices_Set, &edge->to, &hashDummy); 
                            //for startSequenceConstraints in aDS.mergedStartConstraints_Vertices[edge.to]:
                            startSequenceConstraintsCollection = aDS->mergedStartConstraints_Vertices[edge->to]; 
                            newVertices_List = arrayPrepareAppend(newVertices_List, &newVertices_ListSize, newVertices_List_Index + startSequenceConstraintsCollection->length, sizeof(void *));
                            for (k=0; k<startSequenceConstraintsCollection->length; k++) {
                                newVertices_List[newVertices_List_Index++] = startSequenceConstraintsCollection->list[k];
                            } 
                            //for edge2 in aDS.sequenceGraphY.edgesArrangedByFromVertex[edge.to]:
                            illegalEdges_List = arrayPrepareAppend(illegalEdges_List, &illegalEdges_ListSize, illegalEdges_List_Index + aDS->sequenceGraphY->edgesArrangedByFromVertex[edge->to]->length, sizeof(void *));
                            for (k=0; k<aDS->sequenceGraphY->edgesArrangedByFromVertex[edge->to]->length; k++) {  
                                illegalEdges_List[illegalEdges_List_Index++] = aDS->sequenceGraphY->edgesArrangedByFromVertex[edge->to]->list[k];
                            } 
                        }
                    }
                    goto outer;
                }
            }
        }
        newIllegalEdges_List = arrayPrepareAppend(newIllegalEdges_List, &newIllegalEdges_ListSize, newIllegalEdges_List_Index+1, sizeof(void *));
        newIllegalEdges_List[newIllegalEdges_List_Index++] = edge;
    } 
    temp.list = newVertices_List; temp.length = newVertices_List_Index; temp.destructElement = NULL; 
    (*newVertices) = computeUnion(aDS, &temp, vertices, (INT_32 (*)(struct AlignmentDataStructures *, void *, void *))graphMember_VertexComparator);
    
    temp.list = newEdges_List; temp.length = newEdges_List_Index; 
    (*newEdges)  = computeUnion(aDS, &temp, edges, (INT_32 (*)(struct AlignmentDataStructures *, void *, void *))graphMember_EdgeComparator);
    
    (*newIllegalEdges) = copyConstructList(newIllegalEdges_List, newIllegalEdges_List_Index, NULL); //(void (*)(void *))destructEdge); 
    //final memory clean up
    hashtable_destroy(newVertices_Set, FALSE, FALSE); 
}
 
inline LONG_64 v(struct AlignmentDataStructures *aDS, LONG_64 x, LONG_64 y) {
    return (y + x*aDS->vertexYNo)*stateNo();
}


#define MAX_Z 10000

inline LONG_64 vZ(struct AlignmentDataStructures *aDS, LONG_64 x, LONG_64 y, LONG_64 z) {
    return (y + x*aDS->vertexYNo)*(stateNo()*MAX_Z) + (MAX_Z - z - 1)*stateNo();
}

inline void rVZ(struct AlignmentDataStructures *aDS, LONG_64 c, INT_32 *x, INT_32 *y, INT_32 *z, INT_32 *state) { //reverse vertex number function
    INT_32 vertexNo;

    vertexNo = aDS->vertexYNo;
    *state = c % stateNo();
    *z = MAX_Z - (c % (MAX_Z*stateNo()))/stateNo() - 1;
    c /= MAX_Z*stateNo();
    *y = c % vertexNo;
    *x = c / vertexNo;
}

inline INT_32 rVZ_State(struct AlignmentDataStructures *aDS, LONG_64 c) {
    return c % stateNo();
}

inline LONG_64 rVZ_StatelessVertex(struct AlignmentDataStructures *aDS, LONG_64 c) {
    return c - rVZ_State(aDS, c);
}

inline LONG_64 rVZ_Zless(struct AlignmentDataStructures *aDS, LONG_64 c) {
	INT_32 x, y, z, state;
	rVZ(aDS, c, &x, &y, &z, &state);
	return v(aDS, x, y);
}

inline FLOAT_32 *getTempCell(struct AlignmentDataStructures *aDS, LONG_64 vertex, struct Chunks *matrixChunks, struct hashtable *matrix) {
    LONG_64 *chunk;
    FLOAT_32 *cell;
    INT_32 i;
    
    if ((cell = hashtable_search(matrix, &vertex)) == NULL) {    
        chunk = mallocChunk(matrixChunks); 
        *chunk = vertex; 
        cell = (FLOAT_32 *)(chunk + 1); 
        for(i=0; i<stateNo(); i++) {
            cell[i] = LOG_ZERO; 
        } 
        hashtable_insert(matrix, chunk, cell); 
    }
    return cell;
}

inline FLOAT_32 *getCell(struct AlignmentDataStructures *aDS, LONG_64 vertex) {
    LONG_64 *chunk;
    struct List *chunkList;
    INT_32 i;
       
    //search in the binary list
    chunkList = aDS->matrixChunks->chunkList;
    chunk = chunkList->list[chunkList->length-1];
    if(vertex >= *chunk) {
        chunk = bsearch(&vertex, chunk, aDS->matrixChunks->chunkSize - aDS->matrixChunks->remaining, aDS->matrixChunks->elementSize, (int (*)(const void *, const void *))longComparator_Int);
    }
    else {
        for(i=chunkList->length-2; i>=0; i--) {
            chunk = chunkList->list[i];
            if(vertex >= *chunk) {
                break;
            }
        }
        chunk = bsearch(&vertex, chunk, aDS->matrixChunks->chunkSize, aDS->matrixChunks->elementSize, (int (*)(const void *, const void *))longComparator_Int);
    }
    if(chunk == NULL) {
        return NULL;
    }
    return (FLOAT_32 *)(chunk + 1);
}

inline void addCell(struct AlignmentDataStructures *aDS, LONG_64 vertex, FLOAT_32 *cell, LONG_64 *pVertex) {
    assert(*pVertex < vertex);
    *pVertex = vertex;
    LONG_64 *chunk;
    FLOAT_32 *cell2;
    
    chunk = mallocChunk(aDS->matrixChunks); 
    *chunk = vertex; 
    cell2 = (FLOAT_32 *)(chunk + 1); 
    memcpy(cell2, cell, stateNo()*sizeof(FLOAT_32));
}

inline void assignFromTo_Sum(struct AlignmentDataStructures *aDS, INT_32 fromState, INT_32 toState, FLOAT_32 edgeScore) {
    LOG_PLUS_EQUALS(&aDS->toCell[toState], aDS->fromCell[fromState] + edgeScore);
} 

INT_32 getChanges(INT_32 **vertexSequenceCoordinates, struct Edge *edge, INT_32 *changes, INT_32 leafSeqNo, INT_32 leftMostSeqNo) {
    //function to calculate leaf sequence changes along edge
    INT_32 i; 
    INT_32 *sequenceCoordinatesFrom;
    INT_32 *sequenceCoordinatesTo;
    INT_32 seq;
    static INT_32 *activeLeaves;
    static INT_32 activeLeavesSize;
    
    sequenceCoordinatesFrom = vertexSequenceCoordinates[edge->from];
    sequenceCoordinatesTo = vertexSequenceCoordinates[edge->to];
    
    activeLeaves = arrayResize(activeLeaves, &activeLeavesSize, leafSeqNo, sizeof(INT_32));
    memset(activeLeaves, FALSE, sizeof(INT_32)*leafSeqNo);
    getActiveLeaves(activeLeaves, leftMostSeqNo, edge->treeNode);
    //for seq in xrange(0, aDS.leafSeqNoX) {
    i = 0;
    for (seq=0; seq < leafSeqNo; seq++) {
        if ((sequenceCoordinatesTo[seq] - sequenceCoordinatesFrom[seq] > 0) && activeLeaves[seq]) { 
            changes[i++] = seq;
            changes[i++] = sequenceCoordinatesTo[seq];
        }
    } 
    //assert(i > 0); -- this is not true, as delete-delete node will affect this
    return i;
}

INT_32 filterDuplicatesFunction(struct List *graphMembers, void **list) {
    INT_32 i; 
    INT_32 j;
    void *previous;
    void *current;
    
    previous = NULL;
    j=0;
    //for i in xrange(0, len(graphMembers)) {
    for(i=0; i < graphMembers->length; i++) { 
        current = ((struct GraphMemberHolder *)graphMembers->list[i])->graphMember; 
        if (current != previous) { 
            list[j++] = current;
            previous = current;
        }
    } 
    return j;
}

void printTreeNode(struct TreeNode *treeNode) {
    if(treeNode == NULL) {
        return;
    }
    if(!(treeNode->type & TREE_NODE_EFFECTIVELY_SILENT)) {
        printTreeNode(treeNode->treeNodeX);
        printTreeNode(treeNode->treeNodeY);
    }
    if(treeNode->type == TREE_NODE_LEAF) {
        printf(" active leaf " INT_STRING " \n", treeNode->traversalID->leafNo);
    }
}

void debugSets(struct AlignmentDataStructures *aDS, struct List *newEdgesPointer, struct List *newVerticesPointer, struct List *illegalEdgesPointer, INT_32 toX) {
    if(DEBUG) { //debug code
        INT_32 xx; //debug variables
        INT_32 i;
        //printf(" the vertex is " INT_STRING ", of " INT_STRING " " INT_STRING " " INT_STRING " \n", toX, aDS->sequenceGraphX->vertexNo, aDS->sequenceGraphY->vertexNo, aDS->sequenceGraphXSilentVertices_To[toX]);
        //for(xx=0; xx<aDS->leafSeqNoX; xx++) { printf(" coordinate " INT_STRING " \n", ((INT_32*)aDS->vertexXSequenceCoordinates[toX])[xx]); }
        //for(xx=0; xx<newVerticesPointer->length; xx++) { INT_32 *yy = ((struct GraphMemberHolder *)newVerticesPointer->list[xx])->graphMember; printf(" 1351vert " INT_STRING " " , *yy); } printf("\n"); 
        //for(xx=0; xx<newEdgesPointer->length; xx++) { struct Edge *yy = ((struct GraphMemberHolder *)newEdgesPointer->list[xx])->graphMember; printf(" 1352edge " INT_STRING " " INT_STRING " " INT_STRING " " , yy->from, yy->to, yy->silent); } printf("\n"); 
        //for(xx=0; xx<illegalEdgesPointer->length; xx++) { struct Edge *yy = illegalEdgesPointer->list[xx]; printf(" 1353ill " INT_STRING " " INT_STRING " " INT_STRING " " , yy->from, yy->to, yy->silent); } printf("\n");  printf("\n"); 
        INT_32 xxx = -100000;
        INT_32 yyy = -100000;
        INT_32 *prev = NULL;
        INT_32 *current = NULL;
        struct Edge *pEdge = NULL;
        for(xx=0; xx<newEdgesPointer->length; xx++) {
            struct Edge *yy = (struct Edge *)((struct GraphMemberHolder *)newEdgesPointer->list[xx])->graphMember;
            assert(xxx <= yy->to);
            if (yy->to == xxx) {
                assert(yyy <= yy->from);
                if(edgeComparator(pEdge, yy) == 0) {
                    if(memcmp(prev, ((struct GraphMemberHolder *)newEdgesPointer->list[xx])->sequenceConstraints, aDS->leafSeqNoX*sizeof(INT_32)) == 0) {
                        fprintf(stderr, "New edges the same " INT_STRING " " INT_STRING " " INT_STRING " \n", yy->from, yy->to, yy->silent);
                        assert(FALSE);
                    }
                }
            }
            prev = ((struct GraphMemberHolder *)newEdgesPointer->list[xx])->sequenceConstraints;
            pEdge = yy;
            xxx = yy->to;
            yyy = yy->from;
        }
        xxx = -100000;
        for(xx=0; xx<newVerticesPointer->length; xx++) {
            INT_32 yy = *((INT_32 *)((struct GraphMemberHolder *)newVerticesPointer->list[xx])->graphMember);
            assert(xxx <= yy);
            if(xxx == yy) {
                if(memcmp(prev, ((struct GraphMemberHolder *)newVerticesPointer->list[xx])->sequenceConstraints, aDS->leafSeqNoX*sizeof(INT_32)) == 0) {
                    current = ((struct GraphMemberHolder *)newVerticesPointer->list[xx])->sequenceConstraints;
                    for(i=0; i<aDS->leafSeqNoX; i++) {
                        fprintf(stderr, " Failed " INT_STRING " " INT_STRING " \n ", prev[i], current[i]);
                    }
                    fprintf(stderr, " Failed, vertices the same : " INT_STRING "  \n", xxx);
                    assert(FALSE);
                }
            }
            prev = ((struct GraphMemberHolder *)newVerticesPointer->list[xx])->sequenceConstraints;
            xxx = yy;
        }
        //assert(newVerticesPointer->length > 0); (this will not be true when 
        xxx = -100000;
        yyy = -100000;
        for(xx=0; xx<illegalEdgesPointer->length; xx++) {
            struct Edge *yy = (struct Edge *)illegalEdgesPointer->list[xx];
            assert(xxx <= yy->to);
            if (yy->to == xxx) {
                assert(yyy <= yy->from);
                assert(edgeComparator(pEdge, yy));
            }
            pEdge = yy;
            xxx = yy->to;
            yyy = yy->from;
        }
    }  //#end debug code
}

void computeMatrix(struct AlignmentDataStructures *aDS) {
    //memory that needs cleaning up 
    void **previousVertices;
    void **previousEdges;
    void **previousIllegalEdges;
    INT_32 *rightMostVertices;
    //end
    //core structures computed during each loop
    static struct List newVertices;
    static struct List illegalEdges;
    static struct List newEdges;
    
    static struct List filteredVertices;
    static struct List filteredIllegalEdges;
    static struct List filteredEdges;
    
    static struct List scratchList1;
    static struct List scratchList2;
    
    struct List *newVerticesPointer;
    struct List *newEdgesPointer;
    struct List *illegalEdgesPointer;
    
    static void **scratch;
    static INT_32 scratchSize;
    static void **scratch2;
    static INT_32 scratchSize2;
    
    struct Edge *edgeX; 
    struct Edge *edgeY;
    
    struct hashtable *columnCellsHash;
    struct Chunks *columnCellsChunks;
    FLOAT_32* cell;
    
    INT_32 i;
    INT_32 j;
    INT_32 k;
    INT_32 l;
    INT_32 toX;
    INT_32 fromX;
    INT_32 toY;
    INT_32 vertexY;
    INT_32 state;
    INT_32 rightmostVertex; 
    LONG_64 pVertex; 
    
    scratch = arrayResize(scratch, &scratchSize, aDS->sequenceGraphY->edges->length + 1, sizeof(void *)); //can't have duplicates loaded on, so okay
    scratch2 = arrayResize(scratch2, &scratchSize2, aDS->sequenceGraphY->edges->length + 1, sizeof(void *));
    
    previousVertices = callocLocal(aDS->sequenceGraphX->vertexNo, sizeof(void *));
    previousEdges = callocLocal(aDS->sequenceGraphX->vertexNo, sizeof(void *));
    previousIllegalEdges = callocLocal(aDS->sequenceGraphX->vertexNo, sizeof(void *));
    rightMostVertices = callocLocal(aDS->sequenceGraphX->vertexNo, sizeof(INT_32));
    
    //core structures computed during each loop
    copyList(aDS->mergedStartConstraints_Vertices[0], &newVertices);
    copyList(aDS->sequenceGraphY->edgesArrangedByFromVertex[0], &illegalEdges);
    sortEdges(illegalEdges.list, illegalEdges.length); 
    assert(newEdges.length == 0);
    
    pVertex = -1;
    columnCellsHash = create_hashtable(MEDIUM_CHUNK_SIZE, hashtable_longHashKey, hashtable_longEqualKey, NULL, NULL);
    columnCellsChunks = constructChunks(SMALL_CHUNK_SIZE, sizeof(LONG_64) + sizeof(FLOAT_32)*stateNo());
    cell = getTempCell(aDS, 0, columnCellsChunks, columnCellsHash);
    if(aDS->treeStates == NULL || TRUE) {
        for (i=0; i<stateNo(); i++) {
            cell[i] = aDS->startStates[i];
        } 
    }
    else { 
        for (i=0; i<stateNo(); i++) {
            cell[i] = LOG_ZERO; //defensive (should be unnecessary)
        }
        //cell[0] = LOG_ONE;
        cell[aDS->treeStates[aDS->traversalID->mid]] = LOG_ONE;
    }
    //turnOnDeleteXYLoopCorrection(aDS->model);
    //for toX in xrange(0, aDS.sequenceGraphX.vertexNo):
    for (toX=0; toX < aDS->sequenceGraphX->vertexNo; toX++) {  
        if (aDS->sequenceGraphXSilentVertices_To[toX]) { 
            //horizontal incoming transitions
            //for edgeX in aDS.sequenceGraphX.edgesArrangedByToVertex[toX]: {
            for (i=0; i< aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->length; i++) {
                edgeX = aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->list[i];
                fromX = edgeX->from;    
                deleteEndsScratch(aDS, previousVertices[fromX], previousEdges[fromX], previousIllegalEdges[fromX], 
                                  aDS->mergedEndSequenceCoordinates_Edges[edgeX->iD],
                                  aDS->mergedEndSequenceCoordinates_Vertices[edgeX->from],
                                  &filteredVertices, &filteredEdges, &filteredIllegalEdges);     
                k = filterDuplicatesFunction(&filteredVertices, scratch);
                for (j=0; j<k; j++) {
                    vertexY = *((INT_32 *)scratch[j]); 
                    aDS->fromCell = getCell(aDS, v(aDS, fromX, vertexY));
                    aDS->toCell = getTempCell(aDS, v(aDS, toX, vertexY), columnCellsChunks, columnCellsHash);
                    //for state in xrange(0, aDS.stateNo):
                    for (state=0; state<stateNo(); state++) {
                        assignFromTo_Sum(aDS, state, state, edgeX->edgeScore);
                    }  
                } 
                computeUnionScratch(aDS, &filteredVertices, &newVertices, &scratchList1, (INT_32 (*)(struct AlignmentDataStructures *, void *, void *))graphMember_VertexComparator);
                swapListFields(&newVertices, &scratchList1);
               
                computeUnionScratch(aDS, &filteredEdges, &newEdges, &scratchList1, (INT_32 (*)(struct AlignmentDataStructures *, void *, void *))graphMember_EdgeComparator);
                swapListFields(&newEdges, &scratchList1);
                
                computeUnionScratch(aDS, &filteredIllegalEdges, &illegalEdges, &scratchList1, (INT_32 (*)(struct AlignmentDataStructures *, void *, void *))edgeComparatorStub);
                swapListFields(&illegalEdges, &scratchList1);    
            }   
            diffScratch(aDS, &illegalEdges, &newEdges, &scratchList1, (INT_32 (*)(struct AlignmentDataStructures *, void *, void *))edge_graphMember_EdgeComparator);
            swapListFields(&illegalEdges, &scratchList1);
        
            newVerticesPointer = copyConstructList(newVertices.list, newVertices.length, NULL);
            newEdgesPointer = copyConstructList(newEdges.list, newEdges.length, NULL);
            illegalEdgesPointer = copyConstructList(illegalEdges.list, illegalEdges.length, NULL); 
        }  
        else { 
            //for edgeX in aDS.sequenceGraphX.edgesArrangedByToVertex[toX] {
            for (i=0; i< aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->length; i++) {
                edgeX = aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->list[i];
                fromX = edgeX->from; 
                //horizontal incoming transitions
                aDS->noOfChanges = getChanges(aDS->vertexXSequenceCoordinates, edgeX, aDS->changes, aDS->leafSeqNoX, aDS->leftMostSeqNoX);
                filterScratch(aDS, previousVertices[fromX], &scratchList1, (INT_32 (*)(struct AlignmentDataStructures *, void *))greaterThanConditionFilter); 
                filterScratch(aDS, previousEdges[fromX], &scratchList2, (INT_32 (*)(struct AlignmentDataStructures *, void *))greaterThanConditionFilter);
               
                deleteEndsScratch(aDS, &scratchList1, &scratchList2, previousIllegalEdges[fromX], 
                           aDS->mergedEndSequenceCoordinates_Edges[edgeX->iD],
                           aDS->mergedEndSequenceCoordinates_Vertices[edgeX->from],
                           &filteredVertices, &filteredEdges, &filteredIllegalEdges);
                //for vertexY in filterDuplicatesFunction(filteredVertices) {
                k = filterDuplicatesFunction(&filteredVertices, scratch);  
                for (j=0; j<k; j++) { 
                    vertexY = *((INT_32 *)scratch[j]); 
                    if (!aDS->sequenceGraphYSilentVertices_To[vertexY]) { 
                        //x-gap 
                        aDS->fromCell = getCell(aDS, v(aDS, fromX, vertexY)); 
                        aDS->toCell = getTempCell(aDS, v(aDS, toX, vertexY), columnCellsChunks, columnCellsHash);
                        insertXFn(aDS, aDS->model, edgeX, assignFromTo_Sum); 
                        deleteYFn(aDS, aDS->model, edgeX, assignFromTo_Sum); 
                    }
                }
                //matches                                                     
                //for edgeY in filterDuplicatesFunction(filteredEdges) {
                k = filterDuplicatesFunction(&filteredEdges, scratch);
                for (j=0; j<k; j++) {
                    edgeY = scratch[j];
                    if (!edgeY->silent) { 
                        aDS->fromCell = getCell(aDS, v(aDS, fromX, edgeY->from));
                        aDS->toCell = getTempCell(aDS, v(aDS, toX, edgeY->to), columnCellsChunks, columnCellsHash);
                        matchFn(aDS, aDS->model, edgeX, edgeY, assignFromTo_Sum);
                    }
                }
                computeUnionScratch(aDS, &filteredVertices, &newVertices, &scratchList1, (INT_32 (*)(struct AlignmentDataStructures *, void *, void *))graphMember_VertexComparator);
                swapListFields(&newVertices, &scratchList1);
                
                computeUnionScratch(aDS, &filteredEdges, &newEdges, &scratchList1, (INT_32 (*)(struct AlignmentDataStructures *, void *, void *))graphMember_EdgeComparator);
                swapListFields(&newEdges, &scratchList1);
               
                computeUnionScratch(aDS, &filteredIllegalEdges, &illegalEdges, &scratchList1, (INT_32 (*)(struct AlignmentDataStructures *, void *, void *))edgeComparatorStub);
                swapListFields(&illegalEdges, &scratchList1);    
            } 
            diffScratch(aDS, &illegalEdges, &newEdges, &scratchList1, (INT_32 (*)(struct AlignmentDataStructures *, void *, void *))edge_graphMember_EdgeComparator);
            swapListFields(&illegalEdges, &scratchList1);   
            
            extendEnds(aDS, &newVertices, &newEdges, &illegalEdges, 
                       //aDS.mergedEndSequenceCoordinates_Vertices[toX],
                       aDS->mergedEndSequenceCoordinates_Vertices[toX],
                       &newVerticesPointer, &newEdgesPointer, &illegalEdgesPointer); 
            i = filterDuplicatesFunction(newEdgesPointer, scratch);
            j = 0;
            //for toY in filterDuplicatesFunction(newVertices) {
            k = filterDuplicatesFunction(newVerticesPointer, scratch2); 
            for (l=0; l<k; l++) { 
                toY = *((INT_32 *)scratch2[l]);
                aDS->toCell = getTempCell(aDS, v(aDS, toX, toY), columnCellsChunks, columnCellsHash);
                if (aDS->sequenceGraphYSilentVertices_To[toY]) {
                    while (j < i) {
                        edgeY = scratch[j];
                        assert(edgeY->to >= toY);
                        if (edgeY->to == toY) {
                            aDS->fromCell = getTempCell(aDS, v(aDS, toX, edgeY->from), columnCellsChunks, columnCellsHash);
                            //for state in xrange(0, aDS.stateNo):
                            for (state=0; state < stateNo(); state++) {
                                assignFromTo_Sum(aDS, state, state, edgeY->edgeScore);
                            }
                            j += 1;
                        } else {
                            break;
                        }
                    }
                }
                else {  
                    while (j < i) {
                        edgeY = scratch[j];
                        assert(edgeY->to >= toY); 
                        if (edgeY->to == toY) { 
                            aDS->fromCell = getTempCell(aDS, v(aDS, toX, edgeY->from), columnCellsChunks, columnCellsHash);
                            //for state in xrange(0, aDS.stateNo):
                            insertYFn(aDS, aDS->model, edgeY, assignFromTo_Sum);
                            deleteXFn(aDS, aDS->model, edgeY, assignFromTo_Sum);
                            j += 1;
                        } else {
                            break;
                        }
                    }
                    aDS->fromCell = aDS->toCell;
                    silentFn(aDS, aDS->model, aDS->fromCell, assignFromTo_Sum);
                    //deleteFn(aDS, aDS->model, assignFromTo_Sum);
                }
            }
            assert(j == i); 
        }
        k = filterDuplicatesFunction(newVerticesPointer, scratch2); 
        for (l=0; l<k; l++) { 
            toY = *((INT_32 *)scratch2[l]);
            addCell(aDS, v(aDS, toX, toY), getTempCell(aDS, v(aDS, toX, toY), columnCellsChunks, columnCellsHash), &pVertex);
        }
        //clean up loop memory of previous start vertices if possible
        rightmostVertex = 0;//INT_MIN;
        //for edgeX in aDS.sequenceGraphX.edgesArrangedByFromVertex[toX] {
        for (i=0; i < aDS->sequenceGraphX->edgesArrangedByFromVertex[toX]->length; i++) {
            edgeX = aDS->sequenceGraphX->edgesArrangedByFromVertex[toX]->list[i];
            if (edgeX->to > rightmostVertex) {
                rightmostVertex = edgeX->to;
            }
        } 
        //for fromVertex in Set([edgeX.from for edgeX in aDS.sequenceGraphX.edgesArrangedByToVertex[toX]]) {
        j = 0; 
        sortEdges(aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->list, aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->length);
        for (i=0; i<aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->length; i++) {
            edgeX = aDS->sequenceGraphX->edgesArrangedByToVertex[toX]->list[i];
            assert(edgeX->from >= j);
            if (edgeX->from > j) {
                if (rightMostVertices[edgeX->from] == toX) {
                    destructList(previousVertices[edgeX->from]);
                    previousVertices[edgeX->from] = NULL;
                    destructList(previousEdges[edgeX->from]);
                    previousEdges[edgeX->from] = NULL;
                    destructList(previousIllegalEdges[edgeX->from]);
                    previousIllegalEdges[edgeX->from] = NULL;
                }
                j = edgeX->from;
            }
        } 
        debugSets(aDS, newEdgesPointer, newVerticesPointer, illegalEdgesPointer, toX);
        //end loop memory clean up
        //for the next recursion
        previousVertices[toX] = newVerticesPointer; 
        previousEdges[toX] = newEdgesPointer; 
        previousIllegalEdges[toX] = illegalEdgesPointer; 
        
        newVertices.length = 0;
        newEdges.length = 0;
        illegalEdges.length = 0;
        
        filteredVertices.length = 0;
        filteredEdges.length = 0;
        filteredIllegalEdges.length = 0;
        
        scratchList1.length = 0;
        scratchList2.length = 0;
        
        newVerticesPointer = NULL;
        newEdgesPointer = NULL;
        illegalEdgesPointer = NULL; 
        
        destructChunks(columnCellsChunks);
        hashtable_destroy(columnCellsHash, FALSE, FALSE);
        columnCellsHash = create_hashtable(MEDIUM_CHUNK_SIZE, hashtable_longHashKey, hashtable_longEqualKey, NULL, NULL);
        columnCellsChunks = constructChunks(SMALL_CHUNK_SIZE, sizeof(LONG_64) + sizeof(FLOAT_32)*stateNo());//constructEmptyList(0, free);
    }
    //memory clean up
    for(i=0; i<aDS->sequenceGraphX->vertexNo; i++) {
        if(previousVertices[i] != NULL) {
            assert(((struct List *)previousVertices[i])->destructElement == NULL);
            destructList(previousVertices[i]);
        }
        if(previousEdges[i] != NULL) {
            assert(((struct List *)previousEdges[i])->destructElement == NULL);
            destructList(previousEdges[i]);
        }
        if(previousIllegalEdges[i] != NULL) {
            assert(((struct List *)previousIllegalEdges[i])->destructElement == NULL);
            destructList(previousIllegalEdges[i]);        
        }
    }
    free(previousVertices);
    free(rightMostVertices);
    free(previousEdges);
    free(previousIllegalEdges);
    destructChunks(columnCellsChunks);
    hashtable_destroy(columnCellsHash, FALSE, FALSE);
    //end memory clean up
    //end forward recursion 
    j = 1;
    cell = getCell(aDS, v(aDS, aDS->sequenceGraphX->vertexNo-1, aDS->sequenceGraphY->vertexNo-1));
    if(cell != NULL) {
        logDebug("End state probabilities ");
        for (i = 0; i < stateNo(); ++i) {
            logDebug(" " FLOAT_STRING " ", cell[i]);
            if (cell[i] + aDS->endStates[i] > LOG_ZERO) {
                j = 0;
            }
        }
        logDebug("\n");
    }
    else {
        fprintf(stderr, " No path through alignment possible, so I have no choice but to exit, sorry! \n");
        //logDebug(" No path through alignment possible, so I have no choice but to exit, sorry! \n");
        exit(73);
    }
    if(j) {
    	fprintf(stderr, " Zero prob, so I have no choice but to exit, sorry! \n");
    	//logDebug(" No path through alignment possible, so I have no choice but to exit, sorry! \n");
    	exit(73);
    }
}

inline void assignSampling(struct AlignmentDataStructures *aDS, INT_32 fromState, INT_32 toState, FLOAT_32 edgeScore) {
    FLOAT_32 *fromVertices;
    struct TraceBackEdge *edge;
    
    if(aDS->potentialEdges_Index >= aDS->potentialEdges_Size) {
        //resize potential edges
        aDS->potentialEdges = arrayCopyResize(aDS->potentialEdges, &aDS->potentialEdges_Size, aDS->potentialEdges_Size*2, sizeof(void *));
        aDS->potentialEdges_Size /= 2;
        aDS->potentialEdgeCosts = arrayCopyResize(aDS->potentialEdgeCosts, &aDS->potentialEdges_Size, aDS->potentialEdges_Size*2, sizeof(FLOAT_32));
    } 
    fromVertices = aDS->fromCell;
    aDS->potentialEdgeCosts[aDS->potentialEdges_Index] = fromVertices[fromState] + edgeScore;                
    edge = constructTraceBackEdge(aDS->from + fromState, aDS->to + toState, edgeScore, aDS->edgeX, aDS->edgeY, aDS->silent, aDS->getTreeNode);
    aDS->potentialEdges[aDS->potentialEdges_Index++] = edge;
} 
 
inline void assignSampling_CheckState(struct AlignmentDataStructures *aDS, INT_32 fromState, INT_32 toState, FLOAT_32 edgeScore) {
    //neccesary because model
    //may have multiple states for each type of transition
    if (aDS->state == toState) {
        assignSampling(aDS, fromState, toState, edgeScore);
    }
}
 
struct TreeNode *getTreeNode_insertX(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, INT_32 transitionID) {
    return copyConstructTreeNode(TREE_NODE_INSERT, transitionID, 
                                 aDS->traversalID, traceBackEdge->edgeX->treeNode, NULL, NULL);          
}
  
struct TreeNode *getTreeNode_deleteX(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, INT_32 transitionID, FLOAT_32 *wV) {
	copyWV(traceBackEdge->edgeY->wV, wV);
	//multiplyWV(aDS->subModelX->deletionDistribution, traceBackEdge->edgeY->wV, wV);
    return copyConstructTreeNode(TREE_NODE_INTERNAL, transitionID, 
                            aDS->traversalID, aDS->deleteNodeX, traceBackEdge->edgeY->treeNode, NULL);
} 
                             
struct TreeNode *getTreeNode_silentX(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, INT_32 transitionID) {
    return copyConstructTreeNode(TREE_NODE_PREVIOUSLY_SILENT, transitionID, 
                                 aDS->traversalID, traceBackEdge->edgeX->treeNode, NULL, NULL);
}    

struct TreeNode *getTreeNode_insertY(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, INT_32 transitionID) {
    return copyConstructTreeNode(TREE_NODE_INSERT, transitionID, 
                                 aDS->traversalID, traceBackEdge->edgeY->treeNode, NULL, NULL);
}

struct TreeNode *getTreeNode_deleteY(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, INT_32 transitionID, FLOAT_32 *wV) {
    copyWV(traceBackEdge->edgeX->wV, wV);
	//multiplyWV(traceBackEdge->edgeX->wV, aDS->subModelY->deletionDistribution, wV);
    return copyConstructTreeNode(TREE_NODE_INTERNAL, transitionID, 
                                 aDS->traversalID, traceBackEdge->edgeX->treeNode, aDS->deleteNodeY, NULL);
}                             
                             
struct TreeNode *getTreeNode_silentY(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, INT_32 transitionID) {
    return copyConstructTreeNode(TREE_NODE_PREVIOUSLY_SILENT, transitionID, 
                              aDS->traversalID, traceBackEdge->edgeY->treeNode, NULL, NULL);
}

struct TreeNode *getTreeNode_matchXY(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, INT_32 transitionID, FLOAT_32 *wV) {
    multiplyWV(traceBackEdge->edgeX->wV, traceBackEdge->edgeY->wV, wV);
    return copyConstructTreeNode(TREE_NODE_INTERNAL, transitionID, 
                             aDS->traversalID, traceBackEdge->edgeX->treeNode, traceBackEdge->edgeY->treeNode, NULL);
}
                             
struct TreeNode *getTreeNode_deleteXY(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, INT_32 transitionID, FLOAT_32 *wV) {
	copyWV(aDS->subModelX->deletionDistribution, wV);
	//multiplyWV(aDS->subModelX->deletionDistribution, aDS->subModelY->deletionDistribution, wV);
    return copyConstructTreeNode(TREE_NODE_INTERNAL, transitionID, aDS->traversalID, aDS->deleteNodeX, aDS->deleteNodeY, NULL);
}
                             
struct TreeNode *getTreeNode_silentXY(struct AlignmentDataStructures *aDS, struct TraceBackEdge* traceBackEdge, INT_32 transitionID) {
    return copyConstructTreeNode(TREE_NODE_SILENT, transitionID, aDS->traversalID, NULL, NULL, NULL);
}

void updateIndices(struct TreeNode *treeNode, INT_32 *currentIndices, INT_32 leftMostLeafNo) {
    if (treeNode->left != NULL) { //go left
        updateIndices(treeNode->left, currentIndices, leftMostLeafNo);
    }
    if (treeNode->treeNodeX != NULL) {
        updateIndices(treeNode->treeNodeX, currentIndices, leftMostLeafNo);
    }
    if (treeNode->treeNodeY != NULL) {
        updateIndices(treeNode->treeNodeY, currentIndices, leftMostLeafNo);
    }
    if (treeNode->type == TREE_NODE_LEAF) {
        currentIndices[treeNode->traversalID->leafNo - leftMostLeafNo] += 1;
    }
}

INT_32 **calculateVertexSequenceCoordinates(struct SequenceGraph *sequenceGraph, INT_32 leftMostLeafNo, INT_32 leafSeqNo, struct Chunks *chunks) {
    INT_32 vertex;
    struct Edge *previousEdge;
    INT_32 **sequenceCoordinates;
    INT_32 *currentIndices;
    
    //calculates the coordinates for each vertex of the sequences in the graph
    //sequence coordinates are respresented in arrays offset from left-most leaf no
    sequenceCoordinates = mallocLocal(sizeof(INT_32 *)*sequenceGraph->vertexNo); //{ 0:[-1]*leafSeqNo }
    sequenceCoordinates[0] = mallocChunk(chunks); //callocLocal(leafSeqNo, sizeof(INT_32));
    memset(sequenceCoordinates[0], 0, sizeof(INT_32) * leafSeqNo);
    for(vertex=1; vertex<sequenceGraph->vertexNo; vertex++) { // in xrange(1, sequenceGraph.vertexNo):
        previousEdge = sequenceGraph->edgesArrangedByToVertex[vertex]->list[0];
        currentIndices = memcpy(mallocChunk(chunks), sequenceCoordinates[previousEdge->from], sizeof(INT_32)*leafSeqNo);
        updateIndices(previousEdge->treeNode, currentIndices, leftMostLeafNo);
        sequenceCoordinates[vertex] = currentIndices;
    }
    return sequenceCoordinates;
}

void compactVertex(struct List *toEdges, struct List *fromEdges, struct SequenceGraph *sequenceGraph) {
    struct Edge *edge;
    struct TreeNode *treeNode;
    struct Edge *edgeFrom;
    INT_32 i;
    INT_32 j;
    
    static struct List tempList;
     
    listResize(&tempList, toEdges->length * fromEdges->length);
    tempList.length = 0;
    
    for(i=0; i<toEdges->length; i++) {
        edge = toEdges->list[i];
        assert(edge->silent);
    }
    for(i=0; i<toEdges->length-1; i++) {
        edgeFrom = toEdges->list[i];
        for(j=0; j<fromEdges->length; j++) {
            edge = fromEdges->list[j];
            treeNode = edge->treeNode;
            treeNode = copyConstructTreeNode(treeNode->type, treeNode->transitionID, treeNode->traversalID, treeNode->treeNodeX, treeNode->treeNodeY, edge->treeNode->wV);
            treeNode->left = edgeFrom->treeNode;
            edgeFrom->treeNode->refCount++;
            listAppend(sequenceGraph->edgesArrangedByToVertex[edge->to], copyConstructEdge(edgeFrom->from, edge->to, edgeFrom->edgeScore + edge->edgeScore, edge->insertBranchCost, 
                                                                                           edge->deleteBranchCost, edge->wV, edge->silent, treeNode, INT_32_MAX));
        }
    }
    edgeFrom = toEdges->list[toEdges->length-1];
    for(j=0; j<fromEdges->length; j++) {
        edge = fromEdges->list[j];
        assert(edgeFrom->from < edge->to); 
        edge->from = edgeFrom->from;
        edge->edgeScore += edgeFrom->edgeScore;
        assert(edge->treeNode->left == NULL);
        edge->treeNode->left = edgeFrom->treeNode;
        edgeFrom->treeNode->refCount++;
    }
    for(i=0; i<toEdges->length; i++) {
        destructEdge(toEdges->list[i]); //clean up incoming edges
    }
}

struct SequenceGraph *compactSilentVertices(struct SequenceGraph *sequenceGraph) {
    //simplifies silent-edge/vertice subgraphs by removing greedily vertices whose
    //combined incoming and outgoing edges would be better served by direct connections
    INT_32 vertex;
    struct List *toEdges;
    struct List *fromEdges;
    struct List *newEdges;
    INT_32 vertexShift = 0;
    INT_32 *verticeShifts;
    struct Edge *edge;
    INT_32 i;
    struct SequenceGraph *finalSequenceGraph;
    
    verticeShifts = callocLocal(sequenceGraph->vertexNo, sizeof(INT_32));
    newEdges = sequenceGraph->edges;
    newEdges->length = 0;
    for(vertex=0; vertex <sequenceGraph->vertexNo; vertex++) {
        toEdges = sequenceGraph->edgesArrangedByToVertex[vertex];
        fromEdges = sequenceGraph->edgesArrangedByFromVertex[vertex];
        if(fromEdges->length > 0 && toEdges->length > 0 && ((struct Edge *)toEdges->list[0])->silent &&
        //((struct Edge *)fromEdges->list[0])->silent &&
        fromEdges->length + toEdges->length + 1 >= fromEdges->length * toEdges->length) { //the incoming edges are silent, so we can skip them out if needed
            vertexShift++;
            verticeShifts[vertex] = INT_32_MAX;
            compactVertex(toEdges, fromEdges, sequenceGraph);
        } 
        else {
            listAppendArray(newEdges, toEdges->list, toEdges->length);
            verticeShifts[vertex] = vertexShift;
        }   
    }
    for(i=0; i<newEdges->length; i++) {
        edge = newEdges->list[i];
        edge->iD = i;
        edge->from -= verticeShifts[edge->from];
        edge->to -= verticeShifts[edge->to];
    }
    finalSequenceGraph = constructSequenceGraph(newEdges, sequenceGraph->vertexNo-vertexShift);
    
    //memory clean up
    destructSequenceGraph(sequenceGraph, FALSE);
    free(verticeShifts);
    
    return finalSequenceGraph;
}

void convertTransitionIDToStates(INT_32 stateNo, INT_32 z, INT_32 *fromState, INT_32 *toState) {
    *fromState = z/stateNo;
    *toState = z%stateNo;
}

struct SequenceGraph* convertTraceBackEdgesToSequenceGraph(struct AlignmentDataStructures *aDS, void **newEdges, INT_32 newEdgeNumber) {
    INT_32 i;
    INT_32 transitionID;
    struct Edge *edge;
    struct TraceBackEdge *traceBackEdge;
    struct TreeNode *treeNode;
    FLOAT_32 wV[ALPHABET_SIZE];
    struct SequenceGraph *finalSequenceGraph;
    LONG_64 vertex;
    INT_32 compactVertex;
    struct hashtable *vertexMap;
    struct Chunks *intChunks;
    struct Chunks *longChunks;
    
    vertex=LONG_64_MIN;
    compactVertex=0;
    
    intChunks = constructChunks(newEdgeNumber + 1, sizeof(INT_32));
    longChunks = constructChunks(newEdgeNumber + 1, sizeof(LONG_64));
    vertexMap = create_hashtable(newEdgeNumber*2 + 1, hashtable_longHashKey, hashtable_longEqualKey, NULL, NULL);
    hashtable_insert(vertexMap, constructChunkLong(0, longChunks), constructChunkInt(0, intChunks));
    
    for(i=0; i<newEdgeNumber; i++) {
        traceBackEdge = newEdges[i];
        transitionID = stateNo()*rVZ_State(aDS, traceBackEdge->from) + rVZ_State(aDS, traceBackEdge->to);
        if(traceBackEdge->to != vertex) { 
            assert(!hashtable_search(vertexMap, &traceBackEdge->to)); 
            hashtable_insert(vertexMap, constructChunkLong((vertex = traceBackEdge->to), longChunks), constructChunkInt(++compactVertex, intChunks));
            traceBackEdge->to = compactVertex;
        }
        else {
            traceBackEdge->to = compactVertex;
        }
        traceBackEdge->from = *((INT_32 *)hashtable_search(vertexMap, &traceBackEdge->from));
        if(traceBackEdge->silent) {
            treeNode = ((struct TreeNode *(*)(struct AlignmentDataStructures *, struct TraceBackEdge *, INT_32))traceBackEdge->getTreeNode)(aDS, traceBackEdge, transitionID); 
            edge = copyConstructEdge(traceBackEdge->from, traceBackEdge->to, traceBackEdge->edgeScore, 
            LOG_ZERO, LOG_ZERO, NULL, TRUE, treeNode, i);      
        }
        else {
            treeNode = ((struct TreeNode *(*)(struct AlignmentDataStructures *, struct TraceBackEdge *, INT_32, FLOAT_32 *))traceBackEdge->getTreeNode)(aDS, traceBackEdge, transitionID, wV); 
            normaliseWV(wV, wV); //do adjustment for next ancestor, as is convenient at this point to do so
            edge = copyConstructEdge(traceBackEdge->from, traceBackEdge->to, traceBackEdge->edgeScore, 
            LOG_ZERO, LOG_ZERO, //these values not yet set
            wV, FALSE, treeNode, i); 
        }
        destructTraceBackEdge(traceBackEdge);
        newEdges[i] = edge;
    } 
    i = hashtable_count(vertexMap);
    finalSequenceGraph = constructSequenceGraph(copyConstructList((void **)newEdges, newEdgeNumber, (void (*)(void *))destructEdge), i);
    logDebug("Total number of edges and vertices (without compaction) : " INT_STRING " , " INT_STRING " \n", finalSequenceGraph->edges->length, finalSequenceGraph->vertexNo);
    finalSequenceGraph = compactSilentVertices(finalSequenceGraph);
    logDebug("Total number of edges and vertices (after compaction) : " INT_STRING " , " INT_STRING " \n", finalSequenceGraph->edges->length, finalSequenceGraph->vertexNo);
    //memory clean up
    hashtable_destroy(vertexMap, FALSE, FALSE);
    destructChunks(intChunks); 
    destructChunks(longChunks); 
    return finalSequenceGraph;
}

struct SequenceGraph *traceBackMatrix(struct AlignmentDataStructures *aDS) {
    INT_32 i;
    INT_32 j;
    LONG_64 iL;
    INT_32 state;
    INT_32 pathWeight;
    INT_32 pathWeight2;
    INT_32 toX;
    INT_32 toY;
    INT_32 toZ;
    INT_32 fromX;
    INT_32 fromY;
    INT_32 *temp;
    INT_32 hashDummy;
    
    static void **newEdges; //do not clean up returned array
    static INT_32 newEdgesSize;
    INT_32 newEdgesIndex = 0;
    
    struct hashtable *pathWeightsHash;
    struct heap *vertexHeap;
    LONG_64 finalVertex;
    LONG_64 terminationVertex;
    FLOAT_32 *finalVertices;
    FLOAT_32 *endProbs;
    struct hashtable *startStateNos;
    LONG_64 previous;
    INT_32 *choices;
    
    struct List *edgesX;
    struct List *edgesY;
    struct TraceBackEdge *traceBackEdge;
    struct Chunks *intChunks; 
    struct Chunks *longChunks; 
    
    newEdges = arrayResize(newEdges, &newEdgesSize, (aDS->sequenceGraphY->edges->length + aDS->sequenceGraphX->edges->length + stateNo() + 2)*10, sizeof(void *)); //estimated max
    
    finalVertex = vZ(aDS, aDS->sequenceGraphX->vertexNo-1, aDS->sequenceGraphY->vertexNo-1, 0); //starting vertices
    terminationVertex = finalVertex + 100*stateNo(); //vZ(aDS, aDS->sequenceGraphX->vertexNo-1, aDS->sequenceGraphY->vertexNo-1, 0); //starting vertices
    //special vertex used for creating extra vertices for double-deletes
    finalVertices = getCell(aDS, rVZ_Zless(aDS, finalVertex));
    endProbs = mallocLocal(sizeof(FLOAT_32)*stateNo());
    for (state = 0; state < stateNo(); ++state) {
        endProbs[state] = finalVertices[state] + aDS->endStates[state]; 
    }
    //for state, pathWeight in randomChoices(endProbs, aDS.numberOfSamples) {
    vertexHeap = heap_create(aDS->numberOfSamples+MEDIUM_CHUNK_SIZE);
    pathWeightsHash = create_hashtable(MEDIUM_CHUNK_SIZE, hashtable_longHashKey, hashtable_longEqualKey, NULL, NULL);
    choices = randomChoices(endProbs, stateNo(), aDS->numberOfSamples);
    
    intChunks = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(INT_32));
    longChunks = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(LONG_64));
    
    for (state = 0; state < stateNo(); state++) {
        pathWeight = choices[state];
        if (pathWeight > 0) {
            iL = finalVertex + state; 
            heap_insert(vertexHeap, iL);  //heap.pushOnHeap(iL); 
            hashtable_insert(pathWeightsHash, constructChunkLong(iL, longChunks), constructChunkInt(pathWeight, intChunks)); //pathWeightsHash[i] = pathWeight;
            newEdges[newEdgesIndex++] = constructTraceBackEdge(iL, terminationVertex, LOG_ONE, NULL, NULL, TRUE, getTreeNode_silentXY); //end vertex is 0 state
        }
    } 
   
    startStateNos = create_hashtable(stateNo(), hashtable_intHashKey, hashtable_intEqualKey, NULL, NULL);
    //for state in xrange(0, aDS.stateNo) {
    if(aDS->treeStates == NULL || TRUE) {
        for (state=0; state < stateNo(); state++) { 
            if (aDS->startStates[state] > LOG_ZERO) { 
                hashtable_insert(startStateNos, constructChunkInt(state, intChunks), &hashDummy); 
            } 
        }
    }
    else {
        //hashtable_insert(startStateNos, constructChunkInt(0), &hashDummy); 
        hashtable_insert(startStateNos, constructChunkInt(aDS->treeStates[aDS->traversalID->mid], intChunks), &hashDummy); 
    }
    //main loop 
    previous = LONG_64_MAX; 
    while (!heap_empty(vertexHeap)) {
        aDS->toCombined = heap_extract(vertexHeap); // heap.popTopOfHeap()
        //while not heap.empty() and heap.peekTopOfHeap() == aDS.toCombined {
        while (!heap_empty(vertexHeap) && heap_peek(vertexHeap) == aDS->toCombined) {
            heap_extract(vertexHeap); //heap.popTopOfHeap()
        }  
        //uglyf(" current, " LONG_INT_STRING "  previous " LONG_INT_STRING " \n", aDS->toCombined, previous);
        assert(aDS->toCombined <= previous);
        previous = aDS->toCombined;
        aDS->to = rVZ_StatelessVertex(aDS, aDS->toCombined); 
        rVZ(aDS, aDS->toCombined, &toX, &toY, &toZ, &aDS->state);
        assert(getCell(aDS, rVZ_Zless(aDS, aDS->to))[aDS->state] != LOG_ZERO);
        temp = hashtable_remove(pathWeightsHash, &aDS->toCombined, FALSE);
        pathWeight = *temp;
        assert(pathWeight >= 1); 
        if (rVZ_Zless(aDS, aDS->to) == 0 && hashtable_search(startStateNos, &aDS->state) != NULL) {
            continue;
        }
        //for safety
        aDS->silent = FALSE;
        aDS->getTreeNode = NULL;
        aDS->edgeX = NULL;
        aDS->edgeY = NULL;
        aDS->from = INT_32_MAX; 
        //end
        edgesX = aDS->sequenceGraphX->edgesArrangedByToVertex[toX];
        edgesY = aDS->sequenceGraphY->edgesArrangedByToVertex[toY];
        if (aDS->sequenceGraphXSilentVertices_To[toX]) { 
            aDS->silent = TRUE;
            aDS->getTreeNode = getTreeNode_silentX;
            //for edgeX in edgesX {
            for (i = 0; i < edgesX->length; ++i) {
                aDS->edgeX = edgesX->list[i];
                fromX = aDS->edgeX->from;
                aDS->from = vZ(aDS, fromX, toY, 0); //toZ);
                aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                if(aDS->fromCell != NULL) {
                    assignSampling(aDS, aDS->state, aDS->state, 
                                   aDS->edgeX->edgeScore);
                }
            } 
        }
        else if (aDS->sequenceGraphYSilentVertices_To[toY]) { 
            aDS->silent = TRUE;
            aDS->getTreeNode = getTreeNode_silentY;
            //for edgeY in edgesY {
            for (i = 0; i < edgesY->length; ++i) {
                aDS->edgeY = edgesY->list[i];
                fromY = aDS->edgeY->from;
                aDS->from = vZ(aDS, toX, fromY, 0); //toZ);
                aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                if(aDS->fromCell != NULL) {
                    assignSampling(aDS, aDS->state, aDS->state, 
                                   aDS->edgeY->edgeScore);
                }
            } 
        }
        else if (isSilent(aDS->state)) { 
            aDS->silent = TRUE;
            aDS->getTreeNode = getTreeNode_silentXY;
            aDS->from = aDS->to;
            aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
            if(aDS->fromCell != NULL) { 
                silentFn_TraceBack(aDS, aDS->model, assignSampling_CheckState); 
            }
        }
        else if (isXInsert(aDS->state)) { 
            aDS->silent = TRUE;
            aDS->getTreeNode = getTreeNode_insertX;
            //for edgeX in edgesX {
            for (i = 0; i < edgesX->length; ++i) {
                aDS->edgeX = edgesX->list[i];
                fromX = aDS->edgeX->from;
                aDS->from = vZ(aDS, fromX, toY, 0); //toZ);
                aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                if(aDS->fromCell != NULL) {
                    insertXFn(aDS, aDS->model, aDS->edgeX, assignSampling_CheckState);
                }
            } 
        }
        else if (isYDelete(aDS->state)) { 
            aDS->silent = FALSE;
            aDS->getTreeNode = getTreeNode_deleteY;
            //for edgeX in edgesX {
            for (i = 0; i < edgesX->length; ++i) {
                aDS->edgeX = edgesX->list[i];
                fromX = aDS->edgeX->from;
                aDS->from = vZ(aDS, fromX, toY, 0); //toZ);
                aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                if(aDS->fromCell != NULL) {
                    deleteYFn(aDS, aDS->model, aDS->edgeX, assignSampling_CheckState);
                }
            } 
        }
        else if (isYInsert(aDS->state)) { 
            aDS->silent = TRUE;
            aDS->getTreeNode = getTreeNode_insertY;
            //for edgeY in edgesY {
            for (i = 0; i < edgesY->length; ++i) {
                aDS->edgeY = edgesY->list[i];
                fromY = aDS->edgeY->from;
                aDS->from = vZ(aDS, toX, fromY, 0); //toZ);
                aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                if(aDS->fromCell != NULL) {
                    insertYFn(aDS, aDS->model, aDS->edgeY, assignSampling_CheckState);
                }
            } 
        }
        else if (isXDelete(aDS->state)) { 
            aDS->silent = FALSE;
            aDS->getTreeNode = getTreeNode_deleteX;
            //for edgeY in edgesY {
            for (i = 0; i < edgesY->length; ++i) {
                aDS->edgeY = edgesY->list[i];
                fromY = aDS->edgeY->from;
                aDS->from = vZ(aDS, toX, fromY, 0); //toZ);
                aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                if(aDS->fromCell != NULL) {
                    deleteXFn(aDS, aDS->model, aDS->edgeY, assignSampling_CheckState);
                }
            } 
        }
        else if (isMatch(aDS->state)) {  //is match
            aDS->silent = FALSE;
            aDS->getTreeNode = getTreeNode_matchXY;
            //for edgeX in edgesX {
            for (i = 0; i < edgesX->length; ++i) {
                aDS->edgeX = edgesX->list[i];
                fromX = aDS->edgeX->from;
                //for edgeY in edgesY {
                for (j = 0; j < edgesY->length; ++j) {
                    aDS->edgeY = edgesY->list[j];
                    fromY = aDS->edgeY->from;
                    aDS->from = vZ(aDS, fromX, fromY, 0); //toZ);
                    aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
                    if(aDS->fromCell != NULL) {
                        matchFn(aDS, aDS->model, aDS->edgeX, aDS->edgeY, assignSampling_CheckState);
                    } 
                }
            } //is match
        }
        else {  
            //logDebug("delete-delete state " INT_STRING " \n", aDS->state);
            //exit(0);
            assert(isXYDelete(aDS->state));
            aDS->silent = FALSE;
            aDS->getTreeNode = getTreeNode_deleteXY;
            aDS->from = vZ(aDS, toX, toY, toZ+1); 
            if(toZ+1 >= MAX_Z) {
            	fprintf(stderr, "Reached max z-threshold\n");
            	exit(1);
            } 
            aDS->fromCell = getCell(aDS, rVZ_Zless(aDS, aDS->from));
            deleteFn_TraceBack(aDS, aDS->model, assignSampling_CheckState);
        }
        choices = randomChoices(aDS->potentialEdgeCosts, aDS->potentialEdges_Index, pathWeight);
        for (i = 0; i < aDS->potentialEdges_Index; i++) {
            pathWeight2 = choices[i]; 
            traceBackEdge = aDS->potentialEdges[i];
            if (pathWeight2 > 0) {
                newEdges = arrayPrepareAppend(newEdges, &newEdgesSize, newEdgesIndex+1, sizeof(void *));
                //remove the bias introduced by the correction score
                newEdges[newEdgesIndex++] = traceBackEdge;
                heap_insert(vertexHeap, traceBackEdge->from); //heap.pushOnHeap(j);
                
                if ((temp = hashtable_search(pathWeightsHash, &traceBackEdge->from)) != NULL) { 
                    (*temp) += pathWeight2; 
                    //uglyf(" from to3 " LONG_INT_STRING " " LONG_INT_STRING " \n", aDS->from, aDS->to);
                    //pathWeightsHash[aDS.from] += pathWeight2
                }
                else { 
                	 //uglyf(" from to2 " LONG_INT_STRING " " LONG_INT_STRING " \n", aDS->from, aDS->to);
                    hashtable_insert(pathWeightsHash, constructChunkLong(traceBackEdge->from, longChunks), constructChunkInt(pathWeight2, intChunks));
                    //pathWeightsHash[aDS.from] = pathWeight2
                } 
            }
            else { //clean up edge
                aDS->potentialEdges[i] = NULL;
                destructTraceBackEdge(traceBackEdge);
            }
        }
        aDS->potentialEdges_Index = 0;
    }
    //correct start and end conditions
    //for edge in newEdges {
    for (i=0; i<newEdgesIndex; i++) {
        traceBackEdge = newEdges[i];
        j = rVZ_State(aDS, traceBackEdge->from);
        if (rVZ_Zless(aDS, traceBackEdge->from) == 0 && hashtable_search(startStateNos, &j) != NULL) {
            traceBackEdge->from = 0;
        }
    }
    //convert edges to real edges
    for(i=0; i<newEdgesIndex/2; i++) { //alignment.reverse()
        traceBackEdge = newEdges[i];
        newEdges[i] = newEdges[newEdgesIndex-1-i];
        newEdges[newEdgesIndex-1-i] = traceBackEdge; 
    }
    //final memory clean up
    //hashtables
    //hashtable_destroy(newVerticesHash, FALSE, FALSE); 
    hashtable_destroy(pathWeightsHash, FALSE, FALSE); 
    hashtable_destroy(startStateNos, FALSE, FALSE); 
    //heaps
    heap_destroy(vertexHeap);
    //pieces of array
    free(endProbs);
    //chunks arrays
    destructChunks(intChunks);
    destructChunks(longChunks);
    //now convert newEdges to proper edges -- not implemented yet in C
    return convertTraceBackEdgesToSequenceGraph(aDS, newEdges, newEdgesIndex);
}

INT_32 wrapInGraphHolders_Length;
int wrapInGraphHolders_intComp(INT_32 **iAA, INT_32 **iAA2) {
    INT_32 i;
    INT_32 *iA = *iAA;
    INT_32 *iA2 = *iAA2;
    
    for(i=0; i<wrapInGraphHolders_Length; i++) {
        if(iA[i] < iA2[i]) {
            return -1;
        }
        if(iA[i] > iA2[i]) {
            return 1;
        }
    }
    assert(FALSE);
}

void wrapInGraphHolders(struct AlignmentDataStructures *aDS, struct List *list, void *key, INT_32(*comp)(struct AlignmentDataStructures *, struct GraphMemberHolder *, struct GraphMemberHolder *)) {
    INT_32 i;
    
    wrapInGraphHolders_Length = aDS->leafSeqNoX;
    qsort(list->list, list->length, sizeof(void *), (int (*)(const void *, const void *))wrapInGraphHolders_intComp);
    for(i=0; i<list->length; i++) {
        list->list[i] = constructGraphMember(key, list->list[i], NULL);
    }
    if(DEBUG) {
        for(i=1; i<list->length; i++) {
            assert(comp(aDS, list->list[i-1], list->list[i]) < 0);
        }
    }
}

void prepareGraphForAlignment(struct SequenceGraph *sequenceGraphX, struct SubModel **subModels, struct TraversalID *traversalIDX, struct TraversalID *traversalIDY) {
    INT_32 i;
    struct Edge *edge;
    struct SubModel *subModelX;
    struct SubModel *subModelY;
    
    subModelX = subModels[traversalIDX->mid];
    subModelY = subModels[traversalIDY->mid];
    for(i=0; i<sequenceGraphX->edges->length; i++) {
        edge = sequenceGraphX->edges->list[i];
        if(!edge->silent) {
            transformWVByDistance(edge->wV, subModelX->subMatrixBackward, edge->wV);
            edge->insertBranchCost = combineWV(edge->wV, subModelX->insertionDistribution);
            edge->deleteBranchCost = edge->insertBranchCost;
            //edge->deleteBranchCost = combineWV(edge->wV, subModelY->deletionDistribution);
        }
    }
}

struct SequenceGraph *computeEdgeGraph(struct SequenceGraph *sequenceGraphX, struct SequenceGraph *sequenceGraphY, 
                                      struct CombinedTransitionModel *model, 
                                      struct Constraints ***allConstraints,
                                      struct SubModel **subModels, 
                                      struct TraversalID *traversalID,
                                      struct TraversalID *traversalIDX, 
                                      struct TraversalID *traversalIDY,
                                      INT_32 leftMostSeqNo, INT_32 leafSeqNoX, INT_32 leafSeqNoY,
                                      INT_32 numberOfSamples, INT_32 *treeStates) { 
    INT_32 i;
    struct Chunks *mergedChunksX;  
    struct Chunks *mergedChunksY; 
    struct Chunks *intChunks; 
    INT_32 **vertexYSequenceCoordinates;        
    struct SequenceGraph *newSequenceGraph;
    struct AlignmentDataStructures *aDS;
    //computes a specified number of sample alignments using the forward algorithm
    logDebug("Number of edgesX: " INT_STRING " Number of edgesY: " INT_STRING " Number of verticesX: " INT_STRING " Number of verticesY: " INT_STRING "\n",
           sequenceGraphX->edges->length, sequenceGraphY->edges->length, sequenceGraphX->vertexNo, sequenceGraphY->vertexNo);
    
    aDS = mallocLocal(sizeof(struct AlignmentDataStructures));
    
    aDS->sequenceGraphX = sequenceGraphX;
    aDS->sequenceGraphY = sequenceGraphY;
    prepareGraphForAlignment(sequenceGraphX, subModels, traversalIDX, traversalIDY);
    prepareGraphForAlignment(sequenceGraphY, subModels, traversalIDY, traversalIDX);
    
    //struct Edge *edge;
    //for(i=0; i<aDS->sequenceGraphX->edges->length; i++) {
    //    edge = aDS->sequenceGraphX->edges->list[i];
    //    fprintf(stderr, "edgeX " INT_STRING " " INT_STRING " " INT_STRING " " FLOAT_STRING " " INT_STRING " \n", edge->from, edge->to, edge->silent, edge->edgeScore, edge->treeNode);
    //    printTreeNode(edge->treeNode);
    //}
    //for(i=0; i<aDS->sequenceGraphY->edges->length; i++) {
    //    edge = aDS->sequenceGraphY->edges->list[i];
    //    fprintf(stderr, "edgeY " INT_STRING " " INT_STRING " " INT_STRING " " FLOAT_STRING " " INT_STRING " \n", edge->from, edge->to, edge->silent, edge->edgeScore, edge->treeNode);
    //    printTreeNode(edge->treeNode);
    //}
    //void *allConstraints;
    aDS->model = model;
    aDS->subModel = subModels[traversalID->mid];
    
    aDS->startStates = startStates(model);
    aDS->endStates = endStates(model);
    
    //aDS->traversalIDX = traversalIDX;
    //aDS->traversalIDY = traversalIDY;
    aDS->traversalID = traversalID;
    aDS->leftMostSeqNoX = leftMostSeqNo;
    aDS->leftMostSeqNoY = leftMostSeqNo + leafSeqNoX;
    aDS->leafSeqNoX = leafSeqNoX;
    aDS->leafSeqNoY = leafSeqNoY;
    
    aDS->numberOfSamples = numberOfSamples;
    aDS->matrixChunks = constructChunks(LARGE_CHUNK_SIZE, sizeof(LONG_64) + sizeof(FLOAT_32)*stateNo());//constructEmptyList(0, free);
    mergedChunksX = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(INT_32)*aDS->leafSeqNoX);
    mergedChunksY = constructChunks(MEDIUM_CHUNK_SIZE, sizeof(INT_32)*aDS->leafSeqNoY);
    aDS->vertexXSequenceCoordinates = calculateVertexSequenceCoordinates(sequenceGraphX, aDS->leftMostSeqNoX, aDS->leafSeqNoX, mergedChunksX);
    vertexYSequenceCoordinates = calculateVertexSequenceCoordinates(sequenceGraphY, aDS->leftMostSeqNoY, aDS->leafSeqNoY, mergedChunksY);
   
    aDS->sequenceGraphXSilentVertices_To = calculateSilentVertices(sequenceGraphX);//sequenceGraphXSilentVertices_To;
    aDS->sequenceGraphYSilentVertices_To = calculateSilentVertices(sequenceGraphY);//sequenceGraphYSilentVertices_To; 
    //these datastructures are used to compute the constraINT_32 envelope of the alignments during matrix computation
    
    //used to allocate memory
    calculateMergedConstraints_RightToLeft(vertexYSequenceCoordinates,
                                           aDS->leftMostSeqNoY, aDS->leftMostSeqNoX, aDS->leafSeqNoY, aDS->leafSeqNoX, 
                                           allConstraints, sequenceGraphY, 
                                           &aDS->mergedStartConstraints_Vertices, &aDS->mergedStartConstraints_Edges, mergedChunksX);  
    intChunks = constructChunks(sequenceGraphY->vertexNo, sizeof(INT_32));                                                                                             
    for (i=0; i < sequenceGraphY->vertexNo; i++) {
        ((struct List *)aDS->mergedStartConstraints_Vertices[i])->destructElement = (void (*)(void *))destructGraphMember;
        wrapInGraphHolders(aDS, aDS->mergedStartConstraints_Vertices[i], constructChunkInt(i, intChunks), graphMember_VertexComparator); //leaf numbers must be initialised at this point
    }
    for (i=0; i < sequenceGraphY->edges->length; i++) {
        ((struct List *)aDS->mergedStartConstraints_Edges[i])->destructElement = (void (*)(void *))destructGraphMember;
        wrapInGraphHolders(aDS, aDS->mergedStartConstraints_Edges[i], sequenceGraphY->edges->list[i], graphMember_EdgeComparator);
    }
    calculateMergedConstraints_LeftToRight(vertexYSequenceCoordinates,
                                           aDS->leftMostSeqNoY, aDS->leftMostSeqNoX, aDS->leafSeqNoY, aDS->leafSeqNoX, 
                                           allConstraints, sequenceGraphY, 
                                           &aDS->mergedEndConstraints_Vertices, &aDS->mergedEndConstraints_Edges, mergedChunksX);
    calculateMergedSequenceCoordinates_RightToLeft(aDS->vertexXSequenceCoordinates, sequenceGraphX, aDS->leftMostSeqNoX, aDS->leafSeqNoX,
                                                   &aDS->mergedEndSequenceCoordinates_Vertices, &aDS->mergedEndSequenceCoordinates_Edges, 
                                                   mergedChunksX);                                       
    //added to avoid repeated access
    aDS->vertexXNo = sequenceGraphX->vertexNo;
    aDS->vertexYNo = sequenceGraphY->vertexNo;
     
    //used in scanning 
    aDS->newVertices_Size = MEDIUM_CHUNK_SIZE;
    aDS->newVertices = mallocLocal(sizeof(INT_32) * aDS->newVertices_Size);
    aDS->noOfNewVertices = 0;
    aDS->changes = mallocLocal(sizeof(INT_32) * leafSeqNoX * 2);
    aDS->noOfChanges = 0;
    aDS->sequenceCoordinatesCollection = NULL;
    aDS->mergedEndConstraints_GraphMembers = NULL;
    aDS->getIDFromGraphMember = NULL;
    
    //transition starts and ends, used in compute matrix and traceback
    aDS->to = LONG_64_MAX;
    aDS->from = LONG_64_MAX;
    aDS->toCombined = LONG_64_MAX;
    
    //following used by the traceback methods
    //aDS->potentialEdges;
    //aDS->potentialEdgeCosts;
    aDS->potentialEdges_Size = MEDIUM_CHUNK_SIZE;
    aDS->potentialEdges = mallocLocal(sizeof(void *)*aDS->potentialEdges_Size);
    aDS->potentialEdgeCosts = mallocLocal(sizeof(FLOAT_32)*aDS->potentialEdges_Size);
    aDS->potentialEdges_Index = 0;
    //branches for composing tree, used in traceback
    aDS->deleteNodeX = copyConstructTreeNode(TREE_NODE_DELETE, INT_32_MAX, traversalIDX, NULL, NULL, NULL);
    aDS->subModelX = subModels[traversalIDX->mid]; 
    aDS->deleteNodeY = copyConstructTreeNode(TREE_NODE_DELETE, INT_32_MAX, traversalIDY, NULL, NULL, NULL);
    aDS->subModelY = subModels[traversalIDY->mid]; 
    aDS->state = INT_32_MAX;
    aDS->edgeX = NULL;
    aDS->edgeY = NULL;
    aDS->treeStates = treeStates;
    
    //the important calls
    computeMatrix(aDS);
    newSequenceGraph = traceBackMatrix(aDS);
    //memory clean up
    destructChunks(aDS->matrixChunks);
    //destructList(aDS->matrixChunks);
    
    for (i=0; i < sequenceGraphY->vertexNo; i++) {
        destructList(aDS->mergedEndConstraints_Vertices[i]);
        destructList(aDS->mergedStartConstraints_Vertices[i]);
    }
    free(aDS->mergedStartConstraints_Vertices);
    free(aDS->mergedEndConstraints_Vertices);
    destructChunks(intChunks);
    
    for (i=0; i < sequenceGraphY->edges->length; i++) {
        destructList(aDS->mergedEndConstraints_Edges[i]);
        destructList(aDS->mergedStartConstraints_Edges[i]);
    }
    free(aDS->mergedStartConstraints_Edges); 
    free(aDS->mergedEndConstraints_Edges);
  
    cleanUpSequenceCoordinates(aDS->sequenceGraphX, aDS->mergedEndSequenceCoordinates_Vertices,
                               aDS->mergedEndSequenceCoordinates_Edges);
                               
    free(aDS->sequenceGraphXSilentVertices_To);
    free(aDS->sequenceGraphYSilentVertices_To); 
    
    free(vertexYSequenceCoordinates);
    free(aDS->vertexXSequenceCoordinates);
     
    destructChunks(mergedChunksX); 
    destructChunks(mergedChunksY); 
    
    destructTreeNode(aDS->deleteNodeX);
    destructTreeNode(aDS->deleteNodeY);
                
    free(aDS->startStates);
    free(aDS->endStates);
    free(aDS->newVertices);
    free(aDS->changes);
    free(aDS->potentialEdges);
    free(aDS->potentialEdgeCosts);
    free(aDS);
 
    //end memory clean up
    return newSequenceGraph;
}

struct List *viterbi(struct SequenceGraph *sequenceGraph, FLOAT_32 *finalScore) {
    //get the viterbi alignment from the sequence graph
    //assumes alignment starts from 0 vertex and must end in highest 
    //numbered vertex
    //all in order
    FLOAT_32 *viterbiMatrix;
    void **pointers;
    INT_32 i;
    INT_32 j;
    struct Edge *edge;
    FLOAT_32 score;
    void **alignment;
    INT_32 vertex;
    struct List *list; 
    
    viterbiMatrix = mallocLocal(sizeof(FLOAT_32)*sequenceGraph->vertexNo);
    for(i=0; i<sequenceGraph->vertexNo; i++) {
        viterbiMatrix[i] = LOG_ZERO;
    }
    pointers = callocLocal(sequenceGraph->vertexNo, sizeof(INT_32));
    alignment = mallocLocal(sizeof(void *)*sequenceGraph->vertexNo);
    
    //zero is terminator
    pointers[0] = 0;
    viterbiMatrix[0] = LOG_ONE;
    //for edge in sequenceGraph.edges:
    for(i=0; i<sequenceGraph->edges->length; i++) {
        edge = sequenceGraph->edges->list[i];
        score = viterbiMatrix[edge->from] + edge->edgeScore;
        if(score > viterbiMatrix[edge->to]) {
            viterbiMatrix[edge->to] = score;
            pointers[edge->to] = edge;
        }
    } 
    vertex = sequenceGraph->vertexNo-1;
    i=0;
    while(pointers[vertex] != 0) {
        edge = pointers[vertex];
        alignment[i++] = edge;
        vertex = edge->from;
    } 
    *finalScore = viterbiMatrix[sequenceGraph->vertexNo-1];
    for(j=0; j<i/2; j++) { //alignment.reverse()
        edge = alignment[j];
        alignment[j] = alignment[i-1-j];
        alignment[i-1-j] = edge; 
    } 
    list = copyConstructList(alignment, i, NULL);
    //memory clean up
    free(viterbiMatrix);
    free(pointers);
    free(alignment);
    return list;
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//graph iterative/em thing
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

FLOAT_32 *calcForwardMatrix(struct SequenceGraph *sequenceGraph) {
    /*calculates the forward values for the vertices and edges
    *in the graph
    */
    //all in order
    INT_32 i;
    FLOAT_32 *forwardMatrix;
    struct Edge *edge;
    
    forwardMatrix = mallocLocal(sequenceGraph->vertexNo*sizeof(FLOAT_32));
    for(i=0; i<sequenceGraph->vertexNo; i++) {
        forwardMatrix[i] = LOG_ZERO;
    }
    //zero vertex is always terminator
    forwardMatrix[0] = LOG_ONE;
    
    for(i=0; i<sequenceGraph->edges->length; i++) {
        edge = sequenceGraph->edges->list[i];
        //for edge in sequenceGraph.edges:
        LOG_PLUS_EQUALS(&forwardMatrix[edge->to], forwardMatrix[edge->from] + edge->edgeScore);
    }
    return forwardMatrix;
}

FLOAT_32 *calcBackwardMatrix(struct SequenceGraph *sequenceGraph) {
    /*calculates the forward values for the vertices and edges
    *in the graph
    */
    //all in order
    INT_32 i;
    FLOAT_32 *backwardMatrix;
    struct Edge *edge;
    
    backwardMatrix = mallocLocal(sequenceGraph->vertexNo*sizeof(FLOAT_32));
    for(i=0; i<sequenceGraph->vertexNo; i++) {
        backwardMatrix[i] = LOG_ZERO;
    }
    //zero vertex is always terminator
    backwardMatrix[sequenceGraph->vertexNo-1] = LOG_ONE;
    
    for(i=sequenceGraph->edges->length-1; i>=0; i--) {
        edge = sequenceGraph->edges->list[i];
        //for edge in sequenceGraph.edges:
        LOG_PLUS_EQUALS(&backwardMatrix[edge->from], backwardMatrix[edge->to] + edge->edgeScore);
    }
    return backwardMatrix;
}

FLOAT_32 *posteriorProbabilities(FLOAT_32 *forwardMatrix, FLOAT_32 *backwardMatrix, struct SequenceGraph *sequenceGraph) {
    /*
    *calculates the posterior probabilities of edges in the
    *sequence graph
    */
    INT_32 i;
    FLOAT_32 total;
    FLOAT_32 *edgeProbs;
    struct Edge *edge;
    
    total = forwardMatrix[sequenceGraph->vertexNo-1];
    assert(total <= backwardMatrix[0] + 0.00001);
    assert(total >= backwardMatrix[0] - 0.00001);
    edgeProbs = mallocLocal(sequenceGraph->edges->length*sizeof(FLOAT_32));
    for(i=0; i<sequenceGraph->edges->length; i++) {
    //for(edge in sequenceGraph.edges) {
        edge = sequenceGraph->edges->list[i];
        edgeProbs[edge->iD] = forwardMatrix[edge->from] + backwardMatrix[edge->to] + edge->edgeScore - total;
    }
    return edgeProbs;
}

INT_32 tranCoord(INT_32 from, INT_32 to, INT_32 node, INT_32 stateNo, INT_32 nodeNo) {
    return nodeNo * stateNo * stateNo * node + stateNo * from + to;
}

/*void cummulateTreeNode(struct TreeNode *treeNode, FLOAT_32 *transitionProbs, FLOAT_32 prob) {
    LOG_PLUS_EQUALS(&transitionProbs[treeNode->transitionID], prob);
    if(treeNode->left != NULL) {
        cummulateTreeNode(treeNode->left, transitionProbs, prob);
    }
    if(treeNode->treeNodeX != NULL) {
        cummulateTreeNode(treeNode->left, transitionProbs, prob);
    }
    if(treeNode->treeNodeY != NULL) {
        cummulateTreeNode(treeNode->left, transitionProbs, prob);
    }
}*/

/*FLOAT_32 *cummulateTransitionProbs(FLOAT_32 *edgeProbs, INT_32 nodeNumber, INT_32 stateNo, struct SequenceGraph *sequenceGraph) {
    INT_32 i;

    struct Edge *edge;
    FLOAT_32 *transitionProbs;
    FLOAT_32 total;
    
    j = nodeNumber*stateNo*stateNo;
    transitionProbs = mallocLocal(j*sizeof(FLOAT_32));
    for(i=0; i<j; i++) {
        transitionProbs[i] = LOG_ZERO;
    }
    for(i=0; i<sequenceGraph->edges->length; i++) {
    //for(edge in sequenceGraph.edges) {
        edge = sequenceGraph->edges->list[i];
        cummulateTreeNode(edge->treeNode, transitionProbs, edgeProbs[edge->iD]);
    }
    return transitionProbs;
}*/

/*FLOAT_32 *trainParams(struct SequenceGraph *sequenceGraph, INT_32 nodeNumber, INT_32 stateNo) {
    FLOAT_32 *forwardMatrix;
    FLOAT_32 *backwardMatrix;
    FLOAT_32 *posteriorProbs;
    FLOAT_32 *transitionProbs;
    
    INT_32 i;
    INT_32 j;
    INT_32 k;
    
    forwardMatrix = calcForwardMatrix(sequenceGraph);
    backwardMatrix = calcBackwardMatrix(sequenceGraph);
    posteriorProbabilities(forwardMatrix, backwardMatrix, sequenceGraph);
    transitionProbs = cummulateTransitionProbs(posteriorProbs, nodeNumber, stateNo, sequenceGraph);
    free(forwardMatrix);
    free(backwardMatrix);
    free(posteriorProbabilities);
    
    //normalise outgoing probs
    for(i=0; i<nodeNumber; i++) {
        for(j=0; j<stateNo; j++) {
            //total = transitionProbs[tranCoord(i, j, 0, nodeNumber, stateNo)];
            for(k=1; k<stateNo; k++) {
                LOG_PLUS_EQUALS(&total, transitionProbs[tranCoord(i, j, k, nodeNumber, stateNo)]);
            }
            for(k=0; k<stateNo; k++) {
                transitionProbs[tranCoord(i, j, k, nodeNumber, stateNo)] -= total;
            }   
        }
    }
    
    //do regression
    
    
    
    return transitionProbs;
}*/

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//felsensteins
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

FLOAT_32 *upPass(struct TreeNode *treeNode, FLOAT_32 *working, struct SubModel **subModels, FLOAT_32 *totalP);

void branchUp(struct TreeNode *treeNode, FLOAT_32 *working, FLOAT_32 *currentResult, struct SubModel **subModels, FLOAT_32 *totalP) {
    FLOAT_32 *i;
    
    if (treeNode->type != TREE_NODE_DELETE) {
        //return subModels[treeNode->traversalID->mid].transformForwardInTime(upPass(treeNode)); --python error!
        i = upPass(treeNode, working, subModels, totalP);
        transformWVByDistance(i, subModels[treeNode->traversalID->mid]->subMatrixBackward, currentResult);
    }
    else {
        i = subModels[treeNode->traversalID->mid]->deletionDistribution;
        memcpy(currentResult, i, sizeof(FLOAT_32)*ALPHABET_SIZE);
        totalP[treeNode->traversalID->mid] = 0.0f;
    }
    //return subModels[treeNode->traversalID->mid].deletionDistribution 
}

//calculates the un-normalised probabilties of each non-gap residue position
FLOAT_32 *upPass(struct TreeNode *treeNode, FLOAT_32 *working, struct SubModel **subModels, FLOAT_32 *totalP) {
    FLOAT_32 *i;
    FLOAT_32 *j;
    FLOAT_32 *k;
    FLOAT_32 l;
    
    if(treeNode->type == TREE_NODE_INTERNAL) { //is internal treeNode
        i = working + treeNode->traversalID->mid*3*ALPHABET_SIZE;
        j = i + ALPHABET_SIZE;
        k = j + ALPHABET_SIZE;
        branchUp(treeNode->treeNodeX, working, i, subModels, totalP);
        branchUp(treeNode->treeNodeY, working, j, subModels, totalP);
        multiplyWV(i, j, k);
        //uglyf("%i %f %f %f \n", treeNode->traversalID->mid, sumWV(k), totalP[treeNode->treeNodeX->traversalID->mid], totalP[treeNode->treeNodeY->traversalID->mid]);
        totalP[treeNode->traversalID->mid] = LOG(sumWV(k)) + totalP[treeNode->treeNodeX->traversalID->mid] + totalP[treeNode->treeNodeY->traversalID->mid];
        normaliseWV(k, k);
        return k;
    }
    assert(treeNode->type != TREE_NODE_DELETE);
    assert((treeNode->type & TREE_NODE_EFFECTIVELY_SILENT) == 0); //strong than above
    totalP[treeNode->traversalID->mid] = 0.0f;
    l = sumWV(treeNode->wV);
    assert(l < 1.01f);
    assert(l > 0.99f);
    return treeNode->wV;
}

void downPass(struct TreeNode *treeNode, FLOAT_32 *ancestorProb, FLOAT_32 *working, FLOAT_32 *results, struct SubModel **subModels, FLOAT_32 *totalP, FLOAT_32 totalAncP, FLOAT_32 globalTotalP);

void branchDown(struct TreeNode *treeNode, FLOAT_32 *ancestorProbs, FLOAT_32 *working, FLOAT_32 *results, struct SubModel **subModels, FLOAT_32 *totalP, FLOAT_32 totalAncP, FLOAT_32 globalTotalP) {
    transformWVByDistance(ancestorProbs, subModels[treeNode->traversalID->mid]->subMatrixForward, ancestorProbs);
    downPass(treeNode, ancestorProbs, working, results, subModels, totalP, totalAncP, globalTotalP);
}

void downPass(struct TreeNode *treeNode, FLOAT_32 *ancestorProbs, FLOAT_32 *working, FLOAT_32 *results, 
			  struct SubModel **subModels, FLOAT_32 *totalP, FLOAT_32 totalAncP, FLOAT_32 globalTotalP) {
    FLOAT_32 *i;
    FLOAT_32 *j;
    FLOAT_32 *k;
    FLOAT_32 l;
    FLOAT_32 *m;
    
    m = results + treeNode->traversalID->mid*ALPHABET_SIZE;
    if (treeNode->type == TREE_NODE_INTERNAL) { //is internal treeNode
        i = working + treeNode->traversalID->mid*3*ALPHABET_SIZE;
        j = i + ALPHABET_SIZE;
        k = j + ALPHABET_SIZE;
        multiplyWV(ancestorProbs, k, m);
        l = totalP[treeNode->traversalID->mid] + totalAncP + LOG(sumWV(m));
        normaliseWV(m, m);
        //uglyf("%i %f %f \n", treeNode->traversalID->mid, l, globalTotalP);
        assert(l < globalTotalP + 0.01f);
        assert(l > globalTotalP - 0.01f);
        
        if(treeNode->treeNodeX->type != TREE_NODE_DELETE) {
            multiplyWV(ancestorProbs, j, j);
            l = LOG(sumWV(j)) + totalAncP + totalP[treeNode->treeNodeY->traversalID->mid];
            normaliseWV(j, j);
            branchDown(treeNode->treeNodeX, j, working, results, subModels, totalP, l, globalTotalP);
        }
        if(treeNode->treeNodeY->type != TREE_NODE_DELETE) {
            multiplyWV(ancestorProbs, i, i);
            l = LOG(sumWV(i)) + totalAncP + totalP[treeNode->treeNodeX->traversalID->mid];
            normaliseWV(i, i);
            branchDown(treeNode->treeNodeY, i, working, results, subModels, totalP, l, globalTotalP);
        }
    }
    else {
        assert(treeNode->type != TREE_NODE_DELETE);
        assert((treeNode->type & TREE_NODE_EFFECTIVELY_SILENT) == 0); //stronger than above
        memcpy(m, treeNode->wV, sizeof(FLOAT_32)*ALPHABET_SIZE);
        l = sumWV(m);
        //uglyf("%i %f %f \n", treeNode->traversalID->mid, l, globalTotalP);
        assert(l < 1.01f);
        assert(l > 0.99f);
    }
} 

FLOAT_32 *felsensteins(struct TreeNode *treeNode, struct SubModel **subModels, INT_32 nodeNumber) {
    FLOAT_32 *ancestorProbs;
    static FLOAT_32 *working;
    static INT_32 workSize;
    FLOAT_32 *results;
    INT_32 i;
    INT_32 j;
    FLOAT_32 *totalP;
    FLOAT_32 totalAncP;
    FLOAT_32 globalTotalP;
    FLOAT_32 fA[ALPHABET_SIZE];
    FLOAT_32 *fA2;
    
    i = nodeNumber*ALPHABET_SIZE;
    
    working = arrayResize(working, &workSize, i*3, sizeof(FLOAT_32));
    results = mallocLocal(sizeof(FLOAT_32)*i);
    totalP = mallocLocal(sizeof(FLOAT_32)*nodeNumber);
    for(j=0; j<i; j++) { //this is 'gap value'
        results[j] = -1.0f;
    }
    if(treeNode->type == TREE_NODE_INSERT) {
        treeNode = treeNode->treeNodeX;
    }
    totalAncP = 0.0f;
    assert((treeNode->type & TREE_NODE_EFFECTIVELY_SILENT) == 0);
    ancestorProbs = ((struct SubModel *)subModels[treeNode->traversalID->mid])->insertionDistribution;
    fA2 = upPass(treeNode, working, subModels, totalP);
    multiplyWV(ancestorProbs, fA2, fA);
    globalTotalP = totalP[treeNode->traversalID->mid] + LOG(sumWV(fA));
    assert(globalTotalP > -1000000);
    assert(globalTotalP < 1000000);
    downPass(treeNode, ancestorProbs, working, results, subModels, totalP, totalAncP, globalTotalP);
    free(totalP);
    /*for(i=0; i<nodeNumber; i++) {
    	fA2 = results + i*ALPHABET_SIZE;
    	for(j=0; j<ALPHABET_SIZE; j++) {
    		if(fA2[j] >= -0.5) {
    			goto end;
    		}
    	}
    }
    assert(FALSE);
    end:*/
    return results;
}

FLOAT_32 *calculateEdgeInsertionSurvivalBias(FLOAT_32 *ancestorSurvivalMatrix, INT_32 *sequenceLengths, 
                                             INT_32 internalNode, INT_32 alignmentLength, INT_32 seqNo, INT_32 leftMostSeqNo,
                                             INT_32 **alignmentCoordinates, INT_32 **vertexSequenceCoordinates, struct SequenceGraph *sequenceGraph) {
    INT_32 i;
    INT_32 j;
    INT_32 k;
    INT_32 *l;
    INT_32 seq;
    struct Edge *edge;
    INT_32 *changes;
    INT_32 offset;
    FLOAT_32 total;
    FLOAT_32 *ancestorSurvival;
    
    ancestorSurvival = mallocLocal(sizeof(FLOAT_32)*sequenceGraph->edges->length);
    changes = mallocLocal(sizeof(INT_32)*2*seqNo);
    //for each edge, look up edge coordinates
    offset = (internalNode/2)*alignmentLength;
    for(i=0; i<sequenceGraph->edges->length; i++) {
        //do a binary search for each column in the alignment
        //lookup the survival/insertion values, and cummulate them
        edge = sequenceGraph->edges->list[i];
        j = getChanges(vertexSequenceCoordinates, edge, changes, seqNo, leftMostSeqNo);
        total = 0.0;
        for(k=0; k<j; k+=2) {
            seq = changes[i] + leftMostSeqNo;
            l = bsearch(&(changes[i+1]), alignmentCoordinates[seq], sequenceLengths[seq], sizeof(INT_32), (int (*)(const void *, const void *))intComparator_Int);
            total += ancestorSurvivalMatrix[offset + *l];
        }
        total /= k; //(k/2) * 2, because of null model;
        ancestorSurvival[edge->iD] = total;
    }
    //memory clean up
    free(changes);
    //return cummulated tables
    return ancestorSurvival;
}


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//control scripts
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

struct SequenceGraph *align_BottomUpScript(struct BinaryTree *binaryTree, struct SequenceGraph **sequenceGraphs, struct Constraints ***constraints,
                                           char **nodeNames, struct SubModel **subModels, struct CombinedTransitionModel **combinedTransitionModels,
                                           INT_32 numberOfSamples, INT_32 leftMostSeqNo, INT_32 *leafNo, INT_32 *treeStates) {
    struct TraversalID *traversalID;
    char *nodeName;
    struct CombinedTransitionModel *combinedTransitionModel;

    struct BinaryTree *binaryTreeX;
    struct TraversalID *traversalIDX;
    struct SubModel *subModelX;
    char *nodeNameX;
    struct SequenceGraph *sequenceGraphX;
    
    struct BinaryTree *binaryTreeY;
    struct TraversalID *traversalIDY;
    struct SubModel *subModelY;
    char *nodeNameY;
    struct SequenceGraph *sequenceGraphY;
    
    INT_32 leafSeqNoX;
    INT_32 leafSeqNoY;
    
    struct SequenceGraph *sequenceGraph;
    
    traversalID = binaryTree->traversalID;
    nodeName = nodeNames[traversalID->mid];
    logInfo("Starting recursion to create node : %s \n", nodeName);
    if (binaryTree->internal) {
        binaryTreeX = binaryTree->left;
        traversalIDX = binaryTreeX->traversalID;
        subModelX = subModels[traversalIDX->mid];
        nodeNameX = nodeNames[traversalIDX->mid];
        
        binaryTreeY = binaryTree->right;
        traversalIDY = binaryTreeY->traversalID;
        subModelY = subModels[traversalIDY->mid];
        nodeNameY = nodeNames[traversalIDY->mid];
        
        leafSeqNoX = *leafNo;
        sequenceGraphX = align_BottomUpScript(binaryTreeX, sequenceGraphs, constraints, nodeNames, subModels, 
        combinedTransitionModels, numberOfSamples, leftMostSeqNo, leafNo, treeStates);
        leafSeqNoX = *leafNo - leafSeqNoX;
        
        leafSeqNoY = *leafNo;
        sequenceGraphY = align_BottomUpScript(binaryTreeY, sequenceGraphs, constraints, nodeNames, subModels, 
        combinedTransitionModels, numberOfSamples, leftMostSeqNo + leafSeqNoX, leafNo, treeStates);
        leafSeqNoY = *leafNo - leafSeqNoY;
      
        logInfo("Node %s is internal \n", nodeName);
        logInfo("Child x sequence : %s \n", nodeNames[traversalIDX->mid]);
        logInfo("X branch has length " FLOAT_STRING " \n", binaryTreeX->distance);
        logInfo("Child y sequence : %s \n", nodeNames[traversalIDY->mid]);
        logInfo("Y branch has length " FLOAT_STRING " \n", binaryTreeY->distance);
        combinedTransitionModel = combinedTransitionModels[traversalID->mid];
        sequenceGraph = computeEdgeGraph(sequenceGraphX, sequenceGraphY,
                                        combinedTransitionModel, constraints, subModels, 
                                        traversalID, traversalIDX, traversalIDY,
                                        leftMostSeqNo, leafSeqNoX, leafSeqNoY,
                                        numberOfSamples, treeStates);
        logInfo("sequence graph successfull\n");
        //memory clean up of input graphs
        destructSequenceGraph(sequenceGraphX, TRUE);
        destructSequenceGraph(sequenceGraphY, TRUE);
        //end clean up
        return sequenceGraph;                               
    }
    logInfo("Node %s is a leaf sequence \n", nodeName);
    //sequenceGraph = convertSeqToSeqGraph(seqsIt.next(), traversalIDs[binaryTree])
    return sequenceGraphs[(*leafNo)++];
}
