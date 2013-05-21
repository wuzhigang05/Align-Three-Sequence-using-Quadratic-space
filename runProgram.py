#!/usr/bin/env python2.7
import sys
import pdb
import argparse
import random 
import os

if __name__ == '__main__': 
  o = sys.stdout
  e = sys.stderr
  parser= argparse.ArgumentParser(description="Run program in Screen")
  parser.add_argument("program", help="the name of the program you wanna run under screen" )
  parser.add_argument("file", help="the files you wanna do sequende alignment. " +
      "The number of files has to be multiple of three because each three files will " +
      "be compared with each other.",
      nargs = '*')
  args = parser.parse_args()

  if not len(args.file) % 3 == 0:
    print "Error! number of input files have to be multiples of three"
    sys.exit(1)

  for i in range(0, len(args.file), 3):
    A = args.file[i]
    B = args.file[i + 1]
    C = args.file[i + 2]  
    os.system("screen -S T%r -d -m %r/%r %r %r %r" % (i, os.getcwd(), args.program, A, B, C))
#    os.system("screen -S T%r -d -m %r/alignmentThreeSeq.py %r %r %r" % (i, os.getcwd(), A, B, C))
