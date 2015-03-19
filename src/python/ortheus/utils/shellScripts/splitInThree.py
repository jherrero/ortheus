#!/usr/local/bin/python2.3

import sys
import os
import re


def main():
    file = sys.argv[1]
    lines = open(file, 'r').readlines()
    outFile1 = open("%s.1" % file, 'w')
    outFile2 = open("%s.2" % file, 'w')
    outFile3 = open("%s.3" % file, 'w')
    
    for line in lines:
        if "human" in line:
            frags = line.split()
            #print "len", len(frags)
            #print frags
            assert len(frags) == 563
            outFile1.write("%s\n" % " ".join(frags[0:203]))
            l = []
            l2 = []
            for i in xrange(0, 20):
                l.append("1")
                l2.append("1")
                l += frags[203 + i*9: 203 + (i+1)*9]
                l2 += frags[383 + i*9: 383 + (i+1)*9]
            outFile2.write("%s\n" % " ".join(l))
            outFile3.write("%s\n" % " ".join(l2))
            
            #print frags
            #initials = frags[3:203]
            #l = []
            #for i in xrange(0, 20):
            #    l2 = initials[i*10:(i+1)*10]
            #    l += l2[0:8] + l2[9:]
            #outFile1.write("%s\n" % " ".join(frags[0:3] + l))
            #outFile2.write("%s\n" % " ".join(frags[203:383]))
            #outFile3.write("%s\n" % " ".join(frags[383:563]))
            #pass
        else:
            outFile1.write(line)
            outFile2.write(line)
            outFile3.write(line)
    
    outFile1.close()
    outFile2.close()
    outFile3.close()

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
