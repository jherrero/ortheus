#!/usr/local/bin/python2.4

import sys
import os
import os.path
import re
import logging
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

def getTempFile():
    return getTempFile_Global(".comp")
    
def getIndelRates(alignment):
    matches = 0.0
    subs = 0.0
    pGap1 = True
    pGap2 = True
    count1 = 0.0
    count2 = 0.0
    totalCount1 = 0.0
    totalCount2 = 0.0
    length1 = 0.0
    length2 = 0.0
    currentLength1=0
    currentLength2=0
    lengths1=[0]*1001
    lengths2=[0]*1001
    AT = { 'A':0, 'T':0, 'a':0, 't':0 }
    atCount = 0.0
    for column in alignment:
        if column[0] != '-' and column[1] != '-':
            if column[0] == column[1]:
                matches += 1
            else:
                subs += 1
               
        if column[0] == '-':
            currentLength1 += 1
            totalCount1 += 1
        else:
            if column[0] in AT:
                atCount += 1
            length1 += 1
            
        if column[1] == '-':
            currentLength2 += 1
            totalCount2 += 1
        else:
            length2 += 1
        
        if column[0] == '-' and pGap1:
            pGap1 = False
            count1 += 1
        if column[0] != '-' and currentLength1 != 0:
            if currentLength1 < len(lengths1):
                lengths1[currentLength1] += 1
            currentLength1 = 0
        if column[0] != '-':
            pGap1 = True
            
        if column[1] == '-' and pGap2:
            pGap2 = False
            count2 += 1
        if column[1] != '-':
            pGap2 = True
        if column[1] != '-' and currentLength2 != 0:
            if currentLength2 < len(lengths2):
                lengths2[currentLength2] += 1
            currentLength2 = 0
            
    #print "Prop_messed_up: ", (subs*2 + totalCount1 + totalCount2) / (length1 + length2 + 0.00001)
    #print "Gaps: ", (totalCount1 + totalCount2) / (length1 + length2 + 0.00001)
    #print "Subs: ", (subs*2) / (length1 + length2 + 0.0001)
    #print "Correct: ", (matches*2) / (length1 + length2 + 0.0001)
    #return count1, count2, subs, matches  
    return count1, count2, totalCount1, totalCount2, subs, matches, length1, length2, atCount, lengths1, lengths2

#print "stats", subs, matches, length1, length2
#print "Gap_stability:", count1, count2, totalCount1, totalCount2


def compareAncestors(args, seqOne, seqTwo):
    #call pecan
    ortheusPath = ortheusRootPath()
    tempBigAlignment = getTempFile()
    command = "%s %s -F %s %s -E '(a, b);' -G %s %s " % (args.JAVA_PREFIX, args.ALIGNER_PREFIX, seqOne, seqTwo, tempBigAlignment, args.ALIGNMENT_ARGS)
    logger.info("Calling Pecan with : %s", command)
    pecanTime = time.time()
    if os.system(command):
        print "Something went wrong calling aligner, so I've got to go"
        sys.exit(1)
    #print "Completed_alignment %s " % (time.time()-pecanTime)
    i = getIndelRates(multiFastaRead(tempBigAlignment))  
    os.remove(tempBigAlignment)
    return i

def getSubAlignment(alignment, iD):
    seqFile = getTempFile()
    out = open(seqFile, 'w')
    out.write('>\n')
    for column in multiFastaRead(alignment):
        if column[iD] != '-':
            out.write(column[iD])
    out.close()
    return seqFile

def getColumnAlignment(alignment, iDs):
    l = []
    for column in alignment:
        #print "idS", iDs, column
        for iD in iDs:
            if column[iD] != '-':
                l.append([ column[i] for i in iDs ])
                break
    return l

def getColumnAlignments(alignment, ranges, seqNo):
    if ranges == None:
        return [ alignment ]
    l = [ [] for i, j in ranges ]
    
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
             l[fn(index)].append(column[:])
        if column[0] != '-':  
            index += 1
    l2 = []
    for i in l:
        #print "i", i, seqNo
        tempFile = getTempFile()
        writeFastaAlignment(i, [ str(j) for j in xrange(0, seqNo) ], seqNo, tempFile)
        l2.append(tempFile)
    return l2

def compareAncestorsScript(args, alignmentFile1, alignmentFile2, treeString, alignmentScore1, alignmentScore2, selectedSeqs):
    tree = newickTreeParser(treeString, False)
    binaryTree_depthFirstNumbers(tree)
    logger.info("Tree : %s" % printBinaryTree(tree, False))
    
    seqNo = binaryTree_leafNo(tree)
    nodeNo = seqNo * 2 - 1
    
    def fn(tree):
        if tree != None and tree.internal:
            fn(tree.left)
            seqOne = getSubAlignment(alignmentFile1, tree.traversalID.mid)
            seqTwo = getSubAlignment(alignmentFile2, tree.traversalID.mid)
            
            print "Comparing " + printBinaryTree(tree, False)
            count1, count2, totalCount1, totalCount2, subs, matches, length1, length2, atCount, lengths1, lengths2 = compareAncestors(args, seqOne, seqTwo)
            print "Gap_stability:", (count1 + count2), (totalCount1 + totalCount2), "sub_stability", subs, matches, (length1 + length2)
            messedUp = (subs*2 + totalCount1 + totalCount2) / (length1 + length2 + 0.00001)
            gaps = (totalCount1 + totalCount2) / (length1 + length2 + 0.00001)
            subs = (subs*2) / (length1 + length2 + 0.0001)
            correct = (matches*2) / (length1 + length2 + 0.0001)
            print "stats", messedUp, gaps, subs, correct
            
            os.remove(seqOne)
            os.remove(seqTwo)
            fn(tree.right)
    fn(tree)
    
    print "movingOn"
    
    alignments1 = getColumnAlignments(alignmentFile1, selectedSeqs, nodeNo)
    
    def fn2(tree):
        if tree != None and tree.internal:
            fn2(tree.left)
            def fn3(childID, childTree):
                print "Comparing " + printBinaryTree(childTree, False)
                count1 = 0
                count2 = 0
                totalCount1 = 0
                totalCount2 = 0
                subs = 0
                matches = 0
                length1 = 0
                length2 = 0
                atCount = 0
                len1 = [0]*1001
                len2 = [0]*1001
                for l1 in alignments1:
                    #print l1
                    l2 = getColumnAlignment(multiFastaRead(l1), [ childID, tree.traversalID.mid ])
                    #print childID, tree.traversalID.mid, l2
                    count12, count22, totalCount12, totalCount22, subs2, matches2, length12, length22, atCount2, lengths1, lengths2 = getIndelRates(l2)
                    count1 += count12
                    count2 += count22
                    totalCount1 += totalCount12
                    totalCount2 += totalCount22
                    subs += subs2
                    matches += matches2
                    length1 += length12
                    length2 += length22
                    atCount += atCount2 
                    len1 = [ len1[i] + lengths1[i] for i in xrange(0, len(lengths1)) ] 
                    len2 = [ len2[i] + lengths2[i] for i in xrange(0, len(lengths2)) ] 
                print "Gap_stability:", sum(len1[0:11]), sum(len2[0:11]), totalCount1, totalCount2
                print "stats", subs, matches, length1, length2
                #print "lengths1", " ".join([ str(i) for i in lengths1])
                #print "lengths2", " ".join([ str(i) for i in lengths2])
            fn3(tree.left.traversalID.mid, tree.left)
            fn3(tree.right.traversalID.mid, tree.right)
            fn2(tree.right)
    fn2(tree)
    for l1 in alignments1:
        if selectedSeqs != None:
            os.remove(l1)
    
    print "Scores %f %f" % (alignmentScore1, alignmentScore2)
    
def readScore(file):
    i = open(file, 'r')
    j = float(i.readline())
    i.close()
    return j
    
#aligner
def addDefaultCompareAncestorArgs(args):
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
            args.ALIGNMENT_ONE = mods.pop()
            args.ALIGNMENT_SCORE_ONE = readScore(mods.pop())
            continue
        if mod == '-' + indices[1]:
            args.ALIGNMENT_TWO = mods.pop()
            args.ALIGNMENT_SCORE_TWO = readScore(mods.pop())
            continue
        if mod == '-' + indices[2]:
            args.TREE = mods.pop()
            continue
        if mod == '-' + indices[3]:
            file = mods.pop()
            args.RANGES = [ int(i) for i in open(file, 'r').readlines()[0].split() ]
            args.RANGES = [ (args.RANGES[i*2], args.RANGES[i*2 + 1]) for i in xrange(0, len(args.RANGES)/2) ] 
            if len(args.RANGES) == 1 and (args.RANGES[0][1] - args.RANGES[0][0] > 1000000):
                args.RANGES == None 
            continue
        skipped.append(mod)
    return indices[4:]

#read pair of sequences left out of the comparison
#
#read sequence predicted
#run compare ancestors , look at correct, gaps subs etc.

def printMods(args, indices):
    printMod(indices[0], '[STRING] Reconstruction one')
    printMod(indices[1], '[STRING] Reconstruction two')
    printMod(indices[2], '[BRACKETED STRING] Newick tree')
    printMod(indices[3], '[BRACKETED STRING] Sequence ranges')
    return indices[4:]

def main():
    startTime = time.time()
    args = getDefaultArgs()
    addDefaultCompareAncestorArgs(args)
    if len(sys.argv) == 1:
        i = loggerIndices
        print "CompareAncestors.py [MODIFIER_ARGUMENTS]"
        print "A script for comparing reconstructions"
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
    compareAncestorsScript(args, args.ALIGNMENT_ONE, args.ALIGNMENT_TWO, args.TREE, args.ALIGNMENT_SCORE_ONE, args.ALIGNMENT_SCORE_TWO, args.RANGES)
    

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()