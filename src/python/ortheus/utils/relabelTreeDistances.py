#!/usr/local/bin/python2.3

import sys
import os
import re

from benLib.bioio import multiFastaRead
from benLib.bioio import getMultiFastaOffsets
from benLib.bioio import getDefaultLogger
from benLib.bioio import getTempFile as getTempFile_Global
from benLib.bioio import concatanateSeqFiles

from benLib.bioio import getDefaultArgs
from benLib.bioio import printFirstMods
from benLib.bioio import parseFirstMods
from benLib.bioio import logger
from benLib.bioio import printMod
from benLib.bioio import loggerIndices

from benLib.bioio import newickTreeParser
from benLib.bioio import printBinaryTree
from ortheus.localLib.misc import ortheusRootPath
from benLib.bioio import getOpenSeqFiles
from benLib.bioio import closeSeqIterators
from benLib.bioio import removeSeqFiles
from benLib.bioio import writeFastaAlignment

from benLib.tree import binaryTree_depthFirstNumbers
from benLib.tree import binaryTree_nodeNames
from benLib.tree import binaryTree_leafNo

for v in sG.vertices[e.to]:
    for i in sG.edge[x]:
        if i > 0:
            v[i] = 0
            


def printBinaryTreeL(binaryTree, includeDistances, dontStopAtID=True):
    def floatString(f):
        return "%.2f" % f
    def fn(binaryTree):
        #print " tree Node ", binaryTree.left, binaryTree.right, binaryTree.distance, binaryTree.internal, binaryTree.iD 
        if binaryTree.internal and (dontStopAtID or binaryTree.iD == None):
            s = '(' + fn(binaryTree.left) + ',' + fn(binaryTree.right) + ')' + floatString(binaryTree.distance)
        else:
            s = floatString(binaryTree.distance) + " " + str(binaryTree.iD)
        if includeDistances:
            return s + ":" + str(binaryTree.distance)
        return s
    return fn(binaryTree) + ';'

def main():
    tree = newickTreeParser(sys.argv[1], False)
    binaryTree_depthFirstNumbers(tree)
    print "Tree : %s" % printBinaryTree(tree, True)
    list = [ float(i) for i in open(sys.argv[2], 'r').readlines() ]
    list.reverse()
    def fn(tree):
        if tree != None and tree.internal:
            fn(tree.left)
            tree.left.distance = list.pop()
            tree.right.distance = list.pop()
            fn(tree.right)
    fn(tree)
    print list
    print "Tree : %s" % printBinaryTreeL(tree, True)
            
    
    seqNo = binaryTree_leafNo(tree)
    nodeNo = seqNo * 2 - 1
    
    pass

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
    
