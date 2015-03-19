#ifndef SEQUENCEGRAPHC_H_
#define SEQUENCEGRAPHC_H_

#include <stdlib.h>
#include <limits.h>

#include "fastCMaths.h"
#include "commonC.h"
#include "substitutionC.h"

#include "constraintsC.h"
#include "xyzModelC.h" 

//tree nodes
struct TreeNode {
    struct TreeNode *left;
    INT_32 refCount;
    
    INT_32 type;
    INT_32 transitionID;
    struct TraversalID *traversalID;
    struct TreeNode *treeNodeX; 
    struct TreeNode *treeNodeY;
    FLOAT_32 *wV; 
};  

struct TreeNode *copyConstructTreeNode(INT_32 type, INT_32 transitionID, struct TraversalID *traversalID,
                                   struct TreeNode *treeNodeX, struct TreeNode *treeNodeY, FLOAT_32 *wV);                                 

void destructTreeNode(struct TreeNode *treeNode);

#define TREE_NODE_EFFECTIVELY_SILENT 64

#define TREE_NODE_INTERNAL 1
#define TREE_NODE_INSERT 66 //2 | 64;
#define TREE_NODE_DELETE 68 ////4 | 64;
#define TREE_NODE_LEAF 8
#define TREE_NODE_SILENT 80 //16 | 64;
#define TREE_NODE_PREVIOUSLY_SILENT 96 //32 | 64;

//edge

struct Edge {
    //edge in sequence graph  
    
    INT_32 from;
    INT_32 to;
    FLOAT_32 edgeScore;
    FLOAT_32 insertBranchCost;
    FLOAT_32 deleteBranchCost;
    FLOAT_32 *wV;
    INT_32 silent;
    struct TreeNode *treeNode;
    INT_32 iD;
};

struct Edge *copyConstructEdge(INT_32 from, INT_32 to, FLOAT_32 edgeScore, FLOAT_32 insertBranchCost, 
                               FLOAT_32 deleteBranchCost, FLOAT_32 *wV, INT_32 silent, void *treeNode, INT_32 iD);

void destructEdge(struct Edge *edge);

struct TraceBackEdge {
    LONG_64 from;
    LONG_64 to;
    FLOAT_32 edgeScore;
    struct Edge *edgeX;
    struct Edge *edgeY;
    char silent;
    void *getTreeNode;
}; 

struct TraceBackEdge *constructTraceBackEdge(LONG_64 from, LONG_64 to, FLOAT_32 edgeScore, struct Edge *edgeX, struct Edge *edgeY, char silent, void *getTreeNode);

void destructTraceBackEdge(struct TraceBackEdge *edge);

//sequence graphs

struct SequenceGraph {
    INT_32 vertexNo;
    struct List *edges;
    struct List **edgesArrangedByToVertex;
    struct List **edgesArrangedByFromVertex;
};

struct SequenceGraph *constructSequenceGraph(struct List *edges, INT_32 vertexNo);

void destructSequenceGraph(struct SequenceGraph *sequenceGraph, INT_32 freeEdgeList);

//graph member holders

struct GraphMemberHolder {
    void *graphMember;
    INT_32 *sequenceConstraints;
    void (*destructGraphMember)(void *);
};

struct GraphMemberHolder *constructGraphMember(void *graphMember, INT_32 *sequenceConstraints, void (*destructGraphMember)(void *));

void destructGraphMember(struct GraphMemberHolder *graphMemberHolder);

INT_32  isSilent(INT_32 state);
INT_32  isXInsert(INT_32 state);
INT_32  isXDelete(INT_32 state);
INT_32  isYInsert(INT_32 state);
INT_32  isYDelete(INT_32 state);
INT_32  isMatch(INT_32 state);
INT_32  isXYDelete(INT_32 state);

struct CombinedTransitionModel *constructCombinedTransitionModel(FLOAT_32 DELETE_OPEN_PROB,
                                                                 FLOAT_32 DELETE_CONTINIUE_PROB,
                                                                 FLOAT_32 INSERT_OPEN_PROB,
                                                                 FLOAT_32 INSERT_CONTINUE_PROB,
                                                                 FLOAT_32 distanceX, FLOAT_32 distanceY, 
                                                                 //FLOAT_32 lengthX, FLOAT_32 lengthY, 
                                                                 struct SubModel *subModelX, struct SubModel *subModelY);

void destructCombinedTransitionModel(struct CombinedTransitionModel *model);

//end model
//over all data structure

struct AlignmentDataStructures {
    struct SequenceGraph *sequenceGraphX;
    struct SequenceGraph *sequenceGraphY;
    //void *allConstraints;
    struct CombinedTransitionModel *model;
    FLOAT_32 *startStates;
    FLOAT_32 *endStates;
    //void *traversalIDX;
    //void *traversalIDY;
    struct TraversalID *traversalID; 
    INT_32 leftMostSeqNoX;
    INT_32 leftMostSeqNoY;
    INT_32 leafSeqNoX;
    INT_32 leafSeqNoY;
    struct SubModel *subModel;
    INT_32 numberOfSamples;
    
    struct Chunks *matrixChunks;
    FLOAT_32 *fromCell;
    FLOAT_32 *toCell; 
    
    INT_32 **vertexXSequenceCoordinates;
    //INT_32 **vertexYSequenceCoordinates;
    INT_32 *sequenceGraphXSilentVertices_To;
    INT_32 *sequenceGraphYSilentVertices_To;
    
    //these datastructures are used to compute the constraint envelope of the alignments during matrix computation
    void **mergedStartConstraints_Vertices;
    void **mergedStartConstraints_Edges;
    
    void **mergedEndConstraints_Vertices;
    void **mergedEndConstraints_Edges;
    
    void **mergedEndSequenceCoordinates_Vertices;
    void **mergedEndSequenceCoordinates_Edges;
      
    //added to avoid repeated access
    INT_32 vertexXNo;
    INT_32 vertexYNo;
    
    //used in scanning
    INT_32 *newVertices;
    INT_32 newVertices_Size;
    INT_32 noOfNewVertices;
    INT_32 *changes;
    INT_32 noOfChanges;
    struct List *sequenceCoordinatesCollection;
    void **mergedEndConstraints_GraphMembers;
    INT_32 (*getIDFromGraphMember)(void *graphMember);
    
    //transition starts and ends, used in compute matrix and tracebacl
    LONG_64 to;
    LONG_64 from;
    
    //following used by the traceback methods
    void **potentialEdges;
    FLOAT_32 *potentialEdgeCosts;
    INT_32 potentialEdges_Size;
    INT_32 potentialEdges_Index;
    //branches for composing tree, used in traceback
    struct TreeNode *deleteNodeX;
    struct TreeNode *deleteNodeY; 
    struct SubModel *subModelX;
    struct SubModel *subModelY;
    LONG_64 toCombined;
    INT_32 state;
    void *getTreeNode;
    char silent;
    struct Edge *edgeX;
    struct Edge *edgeY;
    INT_32 *treeStates;
};

inline INT_32 stateNo();
 
FLOAT_32 *startStates(struct CombinedTransitionModel *model);

FLOAT_32 *endStates(struct CombinedTransitionModel *model);

//void turnOnDeleteXYLoopCorrection(struct CombinedTransitionModel *model);
    
//void turnOffDeleteXYLoopCorrection(struct CombinedTransitionModel *model);

inline void silentFn(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, FLOAT_32 *cell,
					  void (*assignFn)(struct AlignmentDataStructures *aDS, INT_32, INT_32, FLOAT_32));
  
//inline void silentFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));
    
//inline void deleteFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));
    
//inline void deleteDeleteFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));

inline void silentFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *aDS, INT_32, INT_32, FLOAT_32));

inline void deleteFn_TraceBack(struct AlignmentDataStructures *aDS, struct CombinedTransitionModel *model, void (*assignFn)(struct AlignmentDataStructures *aDS, INT_32, INT_32, FLOAT_32));

inline void insertXFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge, 
                       void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));
    
inline void insertYFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge, 
                      void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));
    
inline void deleteXFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge, 
                      void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));
    
inline void deleteYFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model, struct Edge *edge, 
                      void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));
    
inline void matchFn(struct AlignmentDataStructures *, struct CombinedTransitionModel *model,
            struct Edge *edgeX, struct Edge *edgeY, void (*assignFn)(struct AlignmentDataStructures *, INT_32, INT_32, FLOAT_32));
            
         
INT_32  *randomChoices(FLOAT_32 *probs, INT_32 sizeA, INT_32 pathWeight);  

struct SequenceGraph *computeEdgeGraph(struct SequenceGraph *sequenceGraphX, struct SequenceGraph *sequenceGraphY, 
                        struct CombinedTransitionModel *model,  
                        struct Constraints ***allConstraints,
                          struct SubModel **subModels, 
                          struct TraversalID *traversalID,
                          struct TraversalID *traversalIDX, 
                          struct TraversalID *traversalIDY,
                          INT_32 leftMostSeqNo, INT_32 leafSeqNoX, INT_32 leafSeqNoY,
                          INT_32 numberOfSamples, INT_32 *treeStates); 

struct SequenceGraph *align_BottomUpScript(struct BinaryTree *binaryTree, struct SequenceGraph **sequenceGraphs, struct Constraints ***constraints,
                                           char **nodeNames, struct SubModel **subModels, 
                                           struct CombinedTransitionModel **combinedTransitionModels,
                                           INT_32 numberOfSamples, INT_32 leftMostSeqNo, INT_32 *leafNo, INT_32 *treeStates);
                                           
FLOAT_32 *calculateAncestorSurvivalMatrix(char **alignment, INT_32 alignmentLength, struct BinaryTree *binaryTree, INT_32 seqNo);

FLOAT_32 *calculateEdgeInsertionSurvivalBias(FLOAT_32 *ancestorSurvivalMatrix, INT_32 *sequenceLengths, 
                                             INT_32 internalNode, INT_32 alignmentLength, INT_32 seqNo, INT_32 leftMostSeqNo,
                                             INT_32 **alignmentCoordinates, INT_32 **vertexSequenceCoordinates, struct SequenceGraph *sequenceGraph);

INT_32 **calculateAlignmentCoordinates(char **alignment, INT_32 alignmentLength, INT_32 *seqLengths, INT_32 seqNo);                                                                                        
                                            
struct List *viterbi(struct SequenceGraph *sequenceGraph, FLOAT_32 *finalScore);  

void convertTransitionIDToStates(INT_32 stateNo, INT_32 z, INT_32 *fromState, INT_32 *toState);                                    
 
//felsensteins classic
FLOAT_32 *felsensteins(struct TreeNode *treeNode, struct SubModel **subModels, INT_32 nodeNumber);

#endif /*SEQUENCEGRAPHC_H_*/
