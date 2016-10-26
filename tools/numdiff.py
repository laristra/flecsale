#!/usr/bin/python
#####################################################################
# File: numdiff.py
#
# Description:
#
# Author: Marc R.J. Charest
#
# Date: Friday, April  3 2015
#####################################################################

#!/usr/bin/env python

import sys
import re
from optparse import OptionParser
from math import fabs

# verbosity levels
VERBOSE = 1
VERY_VERBOSE = 2

splitPattern = re.compile(r',|\s+|;')

################################################################################
# The failure exception
################################################################################
class FailObject(object):
    def __init__(self, options):
        self.options = options
        self.failure = False

    def fail(self, brief):
        print( ">>>> " + brief )
        self.failure = True


    def exit(self):
        if (self.failure):
            print( "FAILURE" )
            sys.exit(1)
        else:
            print( "SUCCESS" )
            sys.exit(0)

################################################################################
# split lines into numbers
################################################################################
def numSplit(line):
    list = splitPattern.split(line)
    if list[-1] == "":
        del list[-1]

    numList = [float(a) for a in list]
    return numList

################################################################################
# check for numerical equivalence
################################################################################
def softEquiv(ref, target, relative_tolerance, absolute_tolerance):


    isEquiv = True
    err = fabs(target - ref)
    ref_abs = fabs(ref)

    # if this is greater than machine zero
    try:
        rel_err = err / ref_abs
    except ZeroDivisionError:
        rel_err = float('Inf')
    
    # check relative tolerance
    # .. ignore if negative tolerance
    if ( err > fabs(ref) * relative_tolerance and
            relative_tolerance > 0 ) :
        isEquiv = False

    # check absolute tolerance
    # .. ignore if negative tolerance
    elif ( err > absolute_tolerance and 
            absolute_tolerance > 0 ) :
        isEquiv = False

    # return the results
    return [isEquiv, err, rel_err]

################################################################################
# Compare two strings, they may be text or numbers at this point
################################################################################
def compareStrings(f, options, expLine, actLine, lineNum):


    ### check that they're a bunch of numbers
    try:
        exp = numSplit(expLine)
        act = numSplit(actLine)
    except ValueError as e:

        if options.checkText:

            if (expLine != actLine):
                if options.verbosity:
                    print( "-" * 16 ) 
                    print( "##%-8d<==%s" % (lineNum, expLine) )
                    print( "##%-8d==>%s" % (lineNum, actLine) )
                    print( "@ Text lines do not match!" )
                    print( "-" * 16 )
                f.fail( "Text did not match in line %d" % lineNum )
            else:
                if options.verbosity > VERBOSE:
                    print( "-" * 16 ) 
                    print( "##%-8d<==%s" % (lineNum, expLine) )
                    print( "##%-8d==>%s" % (lineNum, actLine) )
                    print( "@ Text lines match" )
                    print( "-" * 16 )
        
        else:
            if options.verbosity > VERBOSE:
                print( "-" * 16 ) 
                print( "##%-8d<==%s" % (lineNum, expLine) )
                print( "##%-8d==>%s" % (lineNum, actLine) )
                print( "@ Ignoring text lines" )
                print( "-" * 16 )
            
        return

    ### check the ranges
    if len(exp) != len(act):
        if options.verbosity:
            print( "-" * 16 )
            print( "##%-8d<==%s" % (lineNum, expLine) )
            print( "##%-8d==>%s" % (lineNum, actLine) )
            print( "@ Number of columns doesn't match!" )
            print( "-" * 16 )
        f.fail( "Wrong number of columns in line %d" % lineNum )
        return

    ### soft equiv on each value
    for col in range(0, len(exp)):
        expVal = exp[col]
        actVal = act[col]
        [isEquiv, abs_err, rel_err] = softEquiv(expVal, actVal, options.rel_tol, options.abs_tol)

        # for very verbose, always print errors
        if options.verbosity > VERBOSE:
            print( "-" * 16 )
            print( "##%-8d#:%-4d<==%15.8e" % (lineNum, col+1, expVal) )
            print( "##%-8d#:%-4d==>%15.8e" % (lineNum, col+1, actVal) )
            print( "@ Absolute error = %15.8e, Relative error = %15.8e" % (abs_err, rel_err) )
            print( "-" * 16 ) 

        # ERROR
        if not isEquiv:
            if options.verbosity == VERBOSE:
                print( "-" * 16 ) 
                print( "##%-8d#:%-4d<==%15.8e" % (lineNum, col+1, expVal) )
                print( "##%-8d#:%-4d==>%15.8e" % (lineNum, col+1, actVal) )
                print( "@ Absolute error = %15.8e, Relative error = %15.8e" % (abs_err, rel_err) )
                print( "-" * 16 ) 
            f.fail( "Non-equivalence in line %d, column %d" % (lineNum, col) )
    return

################################################################################
# Main Driver
################################################################################
def run(expectedFileName, actualFileName, options):
    # message reporter
    f = FailObject(options)

    expected  = open(expectedFileName)
    actual    = open(actualFileName)
    lineNum   = 0

    while True:

        lineNum += 1
        expLine = expected.readline().rstrip()
        actLine = actual.readline().rstrip()

        if options.verbosity > VERY_VERBOSE:
            print 
            print ( "=" * 16 )
            print ( "Checking line %8d" % (lineNum) )
            print ( "##%-8d<==%s" % (lineNum, expLine) ) 
            print ( "##%-8d==>%s" % (lineNum, actLine) )
            print ( "=" * 16 )
            print 

        ## check that the files haven't ended,
        #  or that they ended at the same time
        if expLine == "":
            if actLine != "":
                f.fail("Tested file ended too late.")
            break
        if actLine == "":
            f.fail("Tested file ended too early.")
            break


        compareStrings(f, options, expLine, actLine, lineNum)

    f.exit()

################################################################################
# Main Function
################################################################################
if __name__ == '__main__':

    parser = OptionParser(usage = "usage: %prog [options] ExpectedFile NewFile")

    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="Don't print status messages to stdout (turns off verbosity).")

    parser.add_option("-v", "--verbose",
                      action="count", dest="verbosity", default=0,
                      help="Set verboseness (-v,-vv,-vvv).")

    parser.add_option("-c", "--check-text",
                      action="store_true", dest="checkText", default=False,
                      help="Verify that lines of text match exactly.")

    parser.add_option("-t", "--tolerance",
                      action="store", type="float", dest="tol", default=1.e-15,
                      help="Relative/absolute error when comparing doubles.")

    parser.add_option("-r", "--relative",
                      action="store", type="float", dest="rel_tol", default=-1.,
                      help="Relative error when comparing doubles.")

    parser.add_option("-a", "--absolute",
                      action="store", type="float", dest="abs_tol", default=-1.,
                      help="Absolute error when comparing doubles.")

    (options, args) = parser.parse_args()

    # print usage
    if len(args) != 2:
        parser.print_help()
        sys.exit(1)


    # check if relative or absolute tolerance was specified
    if options.rel_tol < 0 and options.abs_tol < 0:
        options.abs_tol = options.tol
        options.rel_tol = options.tol

    # quiet option overrides verbosity
    if not options.verbose:
        options.verbosity = 0

    # run diff
    run(args[0], args[1], options)
