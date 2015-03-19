#!/usr/local/bin/python2.3

import sys
import os
import re

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


def main():
    tree = newickTreeParser(sys.argv[1], False)
    binaryTree_depthFirstNumbers(tree)
    seqNo = binaryTree_leafNo(tree)
    nodeNo = seqNo * 2 - 1
    m = []
    def fn(t):
        if t != None:
            if t.internal:
                l = fn(t.left)
                r = fn(t.right)
                for i in l:
                    for j in r:
                        m.append((i[0], j[0], (i[1] + j[1]))) 
                        m.append((j[0], i[0], (i[1] + j[1]))) 
                l = l + r
                for i in l:
                    m.append((i[0], t.traversalID.mid, i[1])) 
                    m.append((t.traversalID.mid, i[0], i[1]))
                l = l + ([t.traversalID.mid, 0.0],)
                for k in l:
                    k[1] += t.distance
                return l
            return ([t.iD, t.distance],)
    fn(tree)
    m.sort()
    for i in m:
        print i


def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
