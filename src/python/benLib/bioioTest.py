import unittest

import os
import random

from bioio import pWMParser
from bioio import writePWMFile
from bioio import getTempFile
from bioio import newickTreeParser
from bioio import printBinaryTree

class TestCase(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testPWMParser(self):
        for test in xrange(0, 1000):
            pWM = getRandomPWM()
            file = getTempFile()
            writePWMFile(pWM, file)
            pWM2 = pWMParser(file)
            for i in xrange(0, len(pWM)):
                pWM[i] == pWM2[i]
            os.remove(file)
            
    def testNewickTreeParser(self):
        d = '((human,baboon),chimp);'
        e = newickTreeParser(d)
        f = printBinaryTree(e, False)
        print d, f
        assert d == f
    
    def testNewickTreeParser_UnaryNodes(self):
        #tests with unary nodes
        for test in xrange(0, 1000):
            tree = getRandomTreeString()
            print "tree to try\t", tree
            tree2 = newickTreeParser(tree, reportUnaryNodes=True)
            tree3 = printBinaryTree(tree2, True)
            print "tree found\t", tree3
            assert tree == tree3
    
def getRandomTreeString():
    def iDFn():
        return random.choice([ "one", "1", "", "he44"])
    def dFn():
        #if random.random() > 0.5:
        return ":%.6f" % random.random()
        #return '' 
    def fn3():
        if random.random() > 0.5:
            if random.random() > 0.5:
                #is unary
                return '(' + fn3() + ')' + iDFn() + dFn()
            else:
                return '(' + fn3() +  ","  + fn3() + ')' + iDFn() + dFn()
        else:
            return iDFn() + dFn()
    return fn3() + ';'
    
def getRandomPWM(length=-1):
    if length == -1:
        length = 1 + int(random.random()*10)
    def fn():
        l = [ random.random()*100 for i in xrange(0, 4) ]
        i = sum(l)
        return [ j/i for j in l ]
    return [ fn() for i in xrange(0, length) ]
        
if __name__ == '__main__':
    unittest.main()