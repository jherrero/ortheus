#!/usr/local/bin/python2.3

import sys
import os
import re

from benLib.bioio import multiFastaRead
from benLib.bioio import  newickTreeParser
from benLib.bioio import  printBinaryTree

from benLib.tree import binaryTree_depthFirstNumbers
from benLib.tree import binaryTree_leafNo

from benLib.misc import linOriginRegression

MAX_LENGTH = 50
MIN_LENGTH = 0

def counts(alignmentFile, treeString):
    tree = newickTreeParser(treeString, False)
    binaryTree_depthFirstNumbers(tree)
    print printBinaryTree(tree, False)
    
    seqNo = binaryTree_leafNo(tree)
    nodeNo = seqNo * 2 - 1
    
    inserts = [ [] for i in xrange(0, nodeNo) ]
    deletes = [ [] for i in xrange(0, nodeNo) ]
    cInserts = [True]*nodeNo
    cDeletes = [True]*nodeNo
    subs = [1.0]*nodeNo
    cons = [1.0]*nodeNo
    totalLengths = [1.0]*nodeNo
    
    def fn(column, tree, j):
        if tree != None:
            k = tree.traversalID.mid
            i = column[k]
            if tree.internal:
                fn(column, tree.left, i)
                fn(column, tree.right, i)
            if i != '-':
                totalLengths[k] += 1
                cDeletes[k] = True
                if j != '-':
                    cInserts[k] = True
                    if i != j:
                        subs[k] += 1
                    else:
                        cons[k] += 1
                else:
                    if cInserts[k]:
                        inserts[k].append(1)
                        cInserts[k] = False
                    else:
                        inserts[k][-1] += 1
            else:
                if j != '-':
                    cInserts[k] = True
                    if cDeletes[k]:
                        deletes[k].append(1)
                        cDeletes[k] = False 
                    else:
                        deletes[k][-1] += 1
    columnNo = 0      
    for column in multiFastaRead(alignmentFile):
        columnNo += 1
        fn(column, tree, '-')
        
    
     
    insertCounts = []
    deleteCounts = []
    def fn2(tree, ancLength):
        if tree != None:
            i = tree.traversalID.mid
            subCount = subs[i]/(cons[i] + subs[i])
            totalSubs = subs[i]
            insertCount = len([ j for j in inserts[i] if j <= MAX_LENGTH and j >= MIN_LENGTH ])
            deleteCount = len([ j for j in deletes[i] if j <= MAX_LENGTH and j >= MIN_LENGTH ])
            
            insertCounts.append((tree.distance + 0.000001, (insertCount + 0.000001) / ancLength))
            deleteCounts.append((tree.distance + 0.000001, (deleteCount + 0.000001) / ancLength))
            
            outputString = "For branch %s, distance %s, expectedInserts %s, expectedDeletes %s, expectedSubs %s totalSubs %s ancLength %s " % (i, tree.distance, insertCount, deleteCount, subCount, totalSubs, ancLength)
            if tree.internal:
                fn2(tree.left, totalLengths[tree.traversalID.mid])
                print outputString
                fn2(tree.right, totalLengths[tree.traversalID.mid])
            else:
                print outputString
    if tree.left.internal:
        fn2(tree.left.left, totalLengths[tree.left.traversalID.mid])
        fn2(tree.left.right, totalLengths[tree.left.traversalID.mid])
    if tree.right.internal:
        fn2(tree.right.left, totalLengths[tree.right.traversalID.mid])
        fn2(tree.right.right, totalLengths[tree.right.traversalID.mid])
    print "insertCounts", insertCounts
    print "deleteCounts", deleteCounts
    print "regression", linOriginRegression(insertCounts)
    print linOriginRegression(deleteCounts)
    print ((linOriginRegression(insertCounts)[0] + linOriginRegression(deleteCounts)[0])/2)

def main():
    if len(sys.argv) != 3:
        print "alignmentFile, treeString"
        sys.exit(0)
    counts(sys.argv[1], sys.argv[2])

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()

