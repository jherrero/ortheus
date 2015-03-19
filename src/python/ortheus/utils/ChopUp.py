#!/usr/local/bin/python2.3

import sys
import os
import re

import sys
import os
import os.path
import re
import logging
import time
import tempfile

from benLib.bioio import multiFastaRead
from benLib.bioio import getDefaultLogger
from benLib.bioio import getTempFile as getTempFile_Global
from benLib.bioio import concatanateSeqFiles
from benLib.bioio import logger
from benLib.bioio import printMod
from ortheus.localLib.misc import ortheusRootPath
from benLib.bioio import getOpenSeqFiles
from benLib.bioio import closeSeqIterators
from benLib.bioio import removeSeqFiles

global tempFile_seqNo
tempFile_seqNo = -1

def getTempFile_Seq():
    global tempFile_seqNo
    tempFile_seqNo += 1
    handle, file = tempfile.mkstemp(".fa", prefix="%03i." % tempFile_seqNo)
    os.close(handle)
    return file

def getTempFile_Align():
    return getTempFile_Global(".mfa")

def getTempFile_Ranges():
    return getTempFile_Global(".ranges")

def getNextAlignmentChunk(alignment, size, seqNo, labels, previousCount, ranges):
    """
    outputs fragment of multiple alignment, and individual sequence files
    """
    global tempFile_seqNo
    tempFile_seqNo = -1
    seqFiles, seqIterators = getOpenSeqFiles(seqNo/2 + 1, getTempFile_Seq)
    for seq in xrange(0, seqNo, 2):
        seqIterators[seq/2].write(">\n")
    alignmentFiles, alignmentIterators = getOpenSeqFiles(seqNo, getTempFile_Seq)
    columnCount = 0
    subRanges = []
    for column in alignment:
        assert column != None
        assert len(column) == seqNo
        #assert len(column) != ['-']*seqNo
        for seq in xrange(0, seqNo):
            residue = column[seq]
            alignmentIterators[seq].write(residue)
            if column[seq] != '-' and (seq % 2) == 0:
                seqIterators[seq/2].write(residue)
        if column[0] != '-':
            columnCount += 1
        if columnCount >= size:
            break
    tempRangesFile = getTempFile_Ranges()
    tempRangesOut = open(tempRangesFile, 'w')
    for i, j in ranges:
        if i < columnCount + previousCount and j >= previousCount:
            tempRangesOut.write(" %i " % (i - previousCount))
            tempRangesOut.write(" %i " % (j - previousCount))
    tempRangesOut.close()
    closeSeqIterators(alignmentIterators, seqNo)
    closeSeqIterators(seqIterators, seqNo/2 + 1)
    tempAlignmentFile = getTempFile_Align()
    concatanateSeqFiles(alignmentFiles, tempAlignmentFile, seqNo, labels)
    removeSeqFiles(alignmentFiles, seqNo)
    if columnCount > 0:
        return seqFiles, tempAlignmentFile, tempRangesFile, columnCount + previousCount
    removeSeqFiles(seqFiles, seqNo/2 + 1)
    os.remove(tempAlignmentFile)
    os.remove(tempRangesFile)
    return None, None, None, previousCount

def parseBed(file, chrom, start, end):
    j = open(file, 'r')
    ranges = []
    for i in j.readlines():
        c, s, e, id = i.split()
        if c == chrom and int(s) >= start and int(e) < end:
            ranges.append((int(s) - start, int(e) - start))
    ranges.sort()
    j.close()
    return ranges

def main():
    if len(sys.argv) == 1:
        print "Chop up alignment into chunks [alignment file] [chunk size] [seqNo in alignment]"
        sys.exit(0)
    alignmentFile = sys.argv[1]
    chunkSize = int(sys.argv[2])
    seqNo = int(sys.argv[3])
    alignIter = multiFastaRead(alignmentFile)
    labels = [ str(i) for i in xrange(0, seqNo) ]
    totalChunkNo = int(sys.argv[5])
    
    if sys.argv[4] != 'None':
        ranges = parseBed(sys.argv[4], sys.argv[5], int(sys.argv[6]), int(sys.argv[7]))
    else:
        ranges = [ (-100000000, 100000000) ]
    
    print ranges
    seqFiles, alignmentFile, rangesFile, columnCount = getNextAlignmentChunk(alignIter, chunkSize, seqNo, labels, 0, ranges)
    chunkNo = 0
    while seqFiles != None and chunkNo < totalChunkNo:
        os.system("mkdir chunks")
        os.system("mv %s chunks" % alignmentFile)
        os.system("mv %s chunks" % rangesFile)
        for i in seqFiles:
            os.system("mv %s chunks" % i)
        if os.system("tar -cvf chunk%s.tar chunks" % chunkNo):
            logger.info("Failed to tar files")
            sys.exit(1)
        os.system("rm -r chunks")
        seqFiles, alignmentFile, rangesFile, columnCount = getNextAlignmentChunk(alignIter, chunkSize, seqNo, labels, columnCount, ranges)
        chunkNo += 1
        tempFile_seqNo = -1
        
    
def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
