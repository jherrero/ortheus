#!/usr/local/bin/python2.3

import sys
import os
import re
import time

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

from bisect import bisect_right

def getColumnAlignments(alignment, ranges, seqNo, outputFile):
    if ranges == None:
        return [ alignment ]
    l = [ ]
    
    linRanges = []
    for i, j in ranges:
        linRanges.append(i)
        linRanges.append(j)
    def fn(index):
        i = bisect_right(linRanges, index)
        i /= 2
        i *= 2
        if i >= 0 and i < len(linRanges) and index >= linRanges[i] and index < linRanges[i+1]:
            return i/2
        return None
    index = 0
    for column in multiFastaRead(alignment):
        #print column, fn(index), ranges
        if fn(index) != None:
             l.append(column[:])
        if column[0] != '-':  
            index += 1
    writeFastaAlignment(l, [ str(j) for j in xrange(0, seqNo) ], seqNo, outputFile)

#aligner
def addDefaultChopUpArgs(args):
    ortheusPath = ortheusRootPath()
    args.JAVA_PREFIX = "java -server "
    args.ALIGNER_PREFIX =  " -cp %s bp.pecan.Pecan " % os.path.join(ortheusPath, "src/python/ortheus/localLib/pecan.jar")
    #args.ALIGNER_PREFIX =  " bp.pecan.Pecan"
    args.ALIGNMENT_ARGS = " -X -d -q -r 1.0 "
    args.RANGES = None
    
def parseMods(mods, args, indices, skipped):
    mods.reverse()
    args.RANGES = None
    while mods:
        mod = mods.pop()
        if mod == '-' + indices[0]:
            args.SEQ_FILE = mods.pop()
            continue
        if mod == '-' + indices[1]:
            file = mods.pop()
            args.RANGES = [ int(i) for i in open(file, 'r').readlines()[0].split() ]
            args.RANGES = [ (args.RANGES[i*2], args.RANGES[i*2 + 1]) for i in xrange(0, len(args.RANGES)/2) ] 
            if len(args.RANGES) == 1 and (args.RANGES[0][1] - args.RANGES[0][0] > 1000000):
                args.RANGES == None 
            continue
        if mod == '-' + indices[2]:
            args.SEQ_NO = int(mods.pop())
            continue
        if mod == '-' + indices[3]:
            args.OUTPUT_FILE = mods.pop()
            continue
        skipped.append(mod)
    return indices[4:]

def printMods(args, indices):
    printMod(indices[0], '[STRING] Seq file')
    printMod(indices[1], '[STRING] Chop up ranges')
    printMod(indices[2], '[INT] Seq no')
    printMod(indices[3], '[STRING] Output file')
    return indices[4:]

def main():
    startTime = time.time()
    args = getDefaultArgs()
    addDefaultChopUpArgs(args)
    if len(sys.argv) == 1:
        i = loggerIndices
        print "CutUpAlignment.py [MODIFIER_ARGUMENTS]"
        print "Arguments:"
        i = printFirstMods(args, i)
        i = printMods(args, i)
        sys.exit(0)
    mods = sys.argv[1:]
    i = loggerIndices
    l = []
    i = parseFirstMods(mods, args, i, l)
    i = parseMods(l, args, i, mods)
    if len(mods) != 0:
        logger.info(" Ooops, remaining arguments %s ", " ".join(mods))
        assert False  
    logger.info("Arguments received : %s " % " ".join(sys.argv))   
    
    columnAlignment = getColumnAlignments(args.SEQ_FILE, args.RANGES, args.SEQ_NO, args.OUTPUT_FILE)

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
    
    