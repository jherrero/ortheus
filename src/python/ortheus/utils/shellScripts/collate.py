#!/usr/local/bin/python2.3

import sys
import os
import re

def main():
    fragmentNo = int(sys.argv[1])
    prefix = sys.argv[2]
    suffix= sys.argv[3]
    lineNo = int(sys.argv[4])
    files = []
    score1 = 0.0
    score2 = 0.0
    subtract = 0
    time1 = 0.0
    for fragment in xrange(0, fragmentNo):
        try:
            #file =  open("%s%s.%s" % (prefix, fragment, suffix), 'r').readlines()
            #file = []
            #if len(file) != 0:
            os.system("grep stats %s%s.%s > one" % (prefix, fragment, suffix))
            os.system("grep Gap_st %s%s.%s > two" % (prefix, fragment, suffix))
            os.system("paste one two > three")
            os.system("sed s/stats// three > two")
            #print "boo1"
            os.system("sed s/Gap_stability:// two > three")
            #print "boo2"
            os.system("sed s/sub_stability// three > two")
            #print "boo3"
            os.system("sed s/\\t// two > one")
            #print "boo4"
            os.system("sed s/\\n// one > three")
            #print "boo5"
            lines = [ i.split() for i in open("three", 'r').readlines() ]
            j = lines[lineNo-1]
            files.append(lines)
            
            os.system("grep Scores %s%s.%s > one" % (prefix, fragment, suffix))
            os.system("sed s/Scores// one > two ")
            line = open("two", 'r').readlines()[0]
            #print "line", line
            j = [ float(i) for i in line.split() ]
            score1 += j[0]
            score2 += j[1]
            
            os.system("grep total_time %s%s.%s > one" % (prefix, fragment, suffix))
            os.system("sed s/total_time// one > two ")
            line = open("two", 'r').readlines()[0]
            #print "line", line
            j = [ float(i) for i in line.split() ]
            time1 += j[0]
        except IOError:
            subtract += 1
        except IndexError:
            subtract += 1
    fragmentNo -= subtract
    ancestors = []
    for fragment in xrange(0, fragmentNo):
        try:
            file =  open("%s%s.%s" % (prefix, fragment, suffix), 'r').readlines()
        except IOError:
            file = []
        if len(file) != 0:
            os.system("grep 'Comparing' %s%s.%s > one" % (prefix, fragment, suffix))
            for line in open("one", 'r').readlines():
                s = line.split()
                #print s
                ancestors.append(s[1])
            break
    s =  " %f %f %f " % (score1, score2, time1)
    for line in xrange(0, lineNo):
        j = [ float(i) for i in files[0][line] ]
        for file in files[1:]:
            k = [ float(i) for i in file[line] ]
            assert len(j) == len(k)
            for l in xrange(0, len(j)):
                j[l] += k[l]
        j = [ str(i / fragmentNo) for i in j ]
        s += " ".join(j) + " " + ancestors[line] + " "
    os.system("rm one two three")
    print s

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
