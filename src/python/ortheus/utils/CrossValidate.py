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


def align(args, seqOne, seqTwo):
    tempBigAlignment = getTempFile()
    command = "%s %s -F %s %s -E '(a, b);' -G %s %s " % (args.JAVA_PREFIX, args.ALIGNER_PREFIX, seqOne, seqTwo, tempBigAlignment, args.ALIGNMENT_ARGS)
    logger.info("Calling Pecan with : %s", command)
    pecanTime = time.time()
    if os.system(command):
        print "Something went wrong calling aligner, so I've got to go"
        sys.exit(1)
    return tempBigAlignment

def compareAncestors(args, seqOne, seqTwo):
    #call pecan
    ortheusPath = ortheusRootPath()
    tempBigAlignment = align(args, seqOne, seqTwo)
    #print "Completed_alignment %s " % (time.time()-pecanTime)
    i = getIndelRates(multiFastaRead(tempBigAlignment))  
    os.remove(tempBigAlignment)
    return i

def compareTwoAncestors(args, seqOne, seqTwo, seqThree):
    #call pecan
    ortheusPath = ortheusRootPath()
    tempBigAlignment1 = align(args, seqOne, seqTwo)
    tempBigAlignment2 = align(args, seqOne, seqThree)
    def intersectionAlign(one, two):
        l = []
        for i in one:
            if i[0] != '-':
                l.append(i[:])
        l2 = []
        for i in two:
            if i[0] != '-':
                l2.append(i[:])
        assert len(l) == len(l2)
        l3 = []
        l4 = []
        for i in xrange(0, len(l)):
            if ('-' not in l[i]) and ('-' not in l2[i]):
                l3.append(l[i])
                l4.append(l2[i])
        return l3, l4
    align1, align2 = intersectionAlign(multiFastaRead(tempBigAlignment1), multiFastaRead(tempBigAlignment2))
    #print "Completed_alignment %s " % (time.time()-pecanTime)
    i = getIndelRates(align1)  
    j = getIndelRates(align2)  
    os.remove(tempBigAlignment1)
    os.remove(tempBigAlignment2)
    return i, j

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

def compareAncestorsScript(args, seqFile, ancFile, cloFile, idString):
    def fn(seqOne, seqTwo, string):
        print "Comparing " + string
        count1, count2, totalCount1, totalCount2, subs, matches, length1, length2, atCount, lengths1, lengths2 = compareAncestors(args, seqOne, seqTwo)
        print "Gap_stability:", count1, count2, "sub_stability", subs, length1, length2
        #messedUp = (subs*2 + totalCount1 + totalCount2) / (length1 + length2 + 0.00001)
        #gaps = (totalCount1 + totalCount2) / (length1 + length2 + 0.00001)
        insertGaps = totalCount2 / (length1 + length2 + 0.00001)
        deleteGaps = totalCount1 / (length1 + length2 + 0.00001)
        subs = (subs*2) / (length1 + length2 + 0.0001)
        correct = (matches*2) / (length1 + length2 + 0.0001)
        print "stats", insertGaps, deleteGaps, subs, correct
    
    fn(seqFile, ancFile, idString + "ancestor")
    fn(seqFile, cloFile, idString + "closest")
    
    def fn2(comp, string):
        print "Comparing " + string + "_intersection"
        count1, count2, totalCount1, totalCount2, subs, matches, length1, length2, atCount, lengths1, lengths2 = comp
        print "Gap_stability:", (count1 + count2), (totalCount1 + totalCount2), "sub_stability", subs, matches, (length1 + length2)
        #messedUp = (subs*2 + totalCount1 + totalCount2) / (length1 + length2 + 0.00001)
        #gaps = (totalCount1 + totalCount2) / (length1 + length2 + 0.00001)
        insertGaps = totalCount2 / (length1 + length2 + 0.00001)
        deleteGaps = totalCount1 / (length1 + length2 + 0.00001)
        subs = (subs*2) / (length1 + length2 + 0.0001)
        correct = (matches*2) / (length1 + length2 + 0.0001)
        print "stats", insertGaps, deleteGaps, subs, correct
    comp = compareTwoAncestors(args, seqFile, ancFile, cloFile)
    fn2(comp[0], idString + "ancestor_int")
    fn2(comp[1], idString + "closest_int")
    
    print "Scores %f %f" % (10, 10)

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
            args.SEQ_FILE = mods.pop()
            continue
        if mod == '-' + indices[1]:
            args.ANC_FILE = mods.pop()
            continue
        if mod == '-' + indices[2]:
            args.CLO_FILE = mods.pop()
            continue
        if mod == '-' + indices[3]:
            args.ID_STRING = mods.pop()
            continue
        skipped.append(mod)
    return indices[3:]

def printMods(args, indices):
    printMod(indices[0], '[STRING] Seq file')
    printMod(indices[1], '[STRING] Anc file')
    printMod(indices[2], '[STRING] Closest extant file')
    printMod(indices[3], '[STRING] ID string')
    return indices[4:]

def main():
    startTime = time.time()
    args = getDefaultArgs()
    addDefaultCompareAncestorArgs(args)
    if len(sys.argv) == 1:
        i = loggerIndices
        print "CrossValidate.py [MODIFIER_ARGUMENTS]"
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
    
    compareAncestorsScript(args, args.SEQ_FILE, args.ANC_FILE, args.CLO_FILE, args.ID_STRING)
    

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()