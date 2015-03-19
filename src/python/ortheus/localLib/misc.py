#!/usr/local/bin/python2.3

import sys
import os
import re


#########################################################
#########################################################
#########################################################
#extenral file functions
#########################################################
#########################################################
#########################################################

def ortheusRootPath():
    """
    function for finding external location
    """
    import dummy
    i = os.path.abspath(dummy.__file__)
    return os.path.split(os.path.split(os.path.split(os.path.split(os.path.split(i)[0])[0])[0])[0])[0]

def main():
    pass

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
