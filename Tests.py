#!/usr/local/bin/python2.3

import sys
import os
import re
import random

import unittest

from benLib.bioio import newickTreeParser

from benLib.bioio import logger

from benLib.bioio import printBinaryTree
from benLib.bioio import multiFastaRead
from benLib.bioio import writeFastaAlignment
from benLib.bioio import getDefaultLogger
from benLib.bioio import getTempFile as getTempFile_Global
from ortheus.localLib.misc import ortheusRootPath

from benLib.tree import binaryTree_depthFirstNumbers
from benLib.tree import BinaryTree

import ortheus.scripts.Nester
PYTHON_STRING = "python"

#JAVA_STRING="~/bin/jre1.6.0/bin/java"
#PECAN_ARGS="-J /lustre/work1/ensembl/bpaten/exonerate-1.4.0-i686/bin/exonerate"

JAVA_STRING="java -server -Xmx500m "
PECAN_ARGS=""

ORTHEUS_ARGS=""
ORTHEUS_PATH = ortheusRootPath()

class TestCase(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
    
    """
    def testENm001(self):
        return True
        #treeString = '(((((((((((((human:0.006969,chimp:0.009727):0.025291,((baboon:0.008968):0.011019):0.024581):0.023649):0.066673):0.018405,((rat:0.081244,mouse:0.072818):0.238435):0.021892):0.02326,(((cow:0.164728,(cat:0.109852,dog:0.107805):0.049576):0.004663):0.010883):0.033242):0.028346):0.016015):0.226853):0.063898):0.126639):0.119814):0.16696);'
        treeString = '((((human:0.006969,chimp:0.009727):0.025291,baboon:0.044568):0.108727,(rat:0.081244,mouse:0.072818):0.260327):0.02326,(cow:0.164728,(cat:0.109852,dog:0.107805):0.049576):0.048788):0.749525;'
        seqFiles = [ "human.ENm001.fa", "chimp.ENm001.fa", "baboon.ENm001.fa", "rat.ENm001.fa", "mouse.ENm001.fa", "cow.ENm001.fa", "cat.ENm001.fa", "dog.ENm001.fa" ]
        seqFiles = [ ORTHEUS_PATH + "/test/ENm001/" + i for i in seqFiles ]
        outputFile = ORTHEUS_PATH + "/test/ENm001/output.mfa"
        command = "%s %s/Ortheus.py -l %s -k '%s' -e '%s' -c '%s' -m %s -i -d '-e'" % \
        (PYTHON_STRING, ORTHEUS_PATH, " ".join(seqFiles), treeString, JAVA_STRING, PECAN_ARGS, outputFile)
        print "running command", command
        if(os.system(command)):
            assert False 
    """   
        
    def testSimulation(self):
        return True
        treeString = '(((((((((((((human:0.006969,chimp:0.009727):0.025291,((baboon:0.008968):0.011019):0.024581):0.023649):0.066673):0.018405,((rat:0.081244,mouse:0.072818):0.238435):0.021892):0.02326,(((cow:0.164728,(cat:0.109852,dog:0.107805):0.049576):0.004663):0.010883):0.033242):0.028346):0.016015):0.226853):0.063898):0.126639):0.119814):0.16696);'
        seqFiles = [ "HUMAN", "CHIMP", "BABOON", "RAT", "MOUSE", "COW", "CAT", "DOG" ]
        seqFiles = [ ORTHEUS_PATH + "/test/" + i for i in seqFiles ]
        outputFile = ORTHEUS_PATH + "/test/output.mfa"
        command = "%s %s/Ortheus.py -m %s -l '%s' -e '%s' -c '%s' -n %s -j -d '-e'" % \
        (PYTHON_STRING, ORTHEUS_PATH, " ".join(seqFiles), treeString, JAVA_STRING, PECAN_ARGS, outputFile)
        print "running command", command
        if(os.system(command)):
            assert False
        
    def testRandom(self):
        MAX_SEQS = 20
        for test in xrange(0, 10):
            print "test no : %i " % test
            #seqNo
            binaryTree = randomTree()
            middleSeq = randomDNASeq(int(random.random()*500))
            seqs = []
            getTreeSeqs(binaryTree, middleSeq, seqs)
            seqNo = len(seqs)
            seqLengths = [ len(i) for i in seqs ]
            #if len(seqs) <= MAX_SEQS and len(seqs) > 2:
            if len(seqs) <= MAX_SEQS and len(seqs) > 2:
                seqFiles = []
                for i in xrange(0, seqNo):
                    seqFiles.append(getTempFile())
                    writeFastaAlignment(seqs[i], str(i), 1, seqFiles[i])
                    print "next sequence ", seqFiles[i], seqLengths[i], seqs[i]
                print "seqs ", seqFiles
                print "seqLengths ", seqLengths
                treeString = printBinaryTree(binaryTree, True)
                print "For tree ", treeString
                #align seqs and check no failure
                outputFile = getTempFile()
                treeFile = getTempFile()
                command = "%s %s/Ortheus.py -e %s -f %s -g %s -m '%s' -l '#%s' -k '#%s' -b" % (PYTHON_STRING, ORTHEUS_PATH, " ".join(seqFiles), outputFile, treeFile, JAVA_STRING, ORTHEUS_ARGS, PECAN_ARGS)
                print "command to call", command
                if(os.system(command)):
                    assert False
                #check alignment is complete
                treeString, seqOrderString = open(treeFile, 'r').readlines()
                orderedSeqFiles = seqOrderString.split()
                print seqFiles, " boo ", orderedSeqFiles
                orderedSeqs = [ seqs[seqFiles.index(orderedSeqFiles[i])] for i in xrange(0, seqNo) ]   
                alignment = [ i[:] for i in multiFastaRead(outputFile) ]
                #print "alignment", alignment
                checkAlignment(alignment, orderedSeqs)
                #check alignment is constrained
                pass
                #clean up
                os.remove(outputFile)
                os.remove(treeFile)
                for i in xrange(0, seqNo):
                    os.remove(seqFiles[i])
                print "test no is finished : %i " % test

def randomTree():
    leafNo = [-1]
    def fn():
        if random.random() > 0.6:
            return BinaryTree(random.random()*0.8, True, fn(), fn(), None)
        else:
            leafNo[0] += 1
            return BinaryTree(random.random()*0.8, False, None, None, str(leafNo[0]))
    return BinaryTree(random.random(), True, fn(), fn(), None)

def getTreeSeqs(binaryTree, seq, l):
    seq = mutateDNASeq(seq, binaryTree.distance)
    if binaryTree.internal:
        getTreeSeqs(binaryTree.left, seq, l)
        getTreeSeqs(binaryTree.right, seq, l)
    else:
        l.append(seq)

def getTempFile():
    return getTempFile_Global(".test")

def randomDNASeq(length):
    return [ random.choice([ 'A', 'C', 'T', 'G' ]) for i in xrange(0, length) ]

def expLength(i=0, prob=0.95):
    if random.random() >= prob:
        return expLength(i+1)
    return i

def mutateDNASeq(seq, distance):
    subProb=distance
    inProb=0.1*distance
    deProb=0.1*distance
    contProb=0.9
    l = []
    bases = [ 'A', 'C', 'T', 'G' ]
    i=0
    while i < len(seq):
        if random.random() < subProb:
            l.append(random.choice(bases))
        else:
            l.append(seq[i])
        if random.random() < inProb:
            l += randomDNASeq(expLength(0, contProb))
        if random.random() < deProb:
            i += int(expLength(0, contProb))
        i += 1
    return l

def checkAlignment(align, seqs):
    #for i in align:
    #    assert i != [GAP]*len(seqs)
    i = [0]*len(seqs)
    for j in align:
        for k in xrange(0, len(seqs)):
            if j[k*2] != '-':
                #print "grr", j[k*2], seqs[k][i[k]]
                assert j[k*2] == seqs[k][i[k]]
                i[k] += 1
    for j in xrange(0, len(seqs)):
        assert i[j] == len(seqs[j])
        
if __name__ == '__main__':
    unittest.main()

def main():
    pass

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()