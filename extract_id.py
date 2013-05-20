#!/usr/bin/env python2.7
import sys
import pdb
import argparse
import random 
import re

def line_yielder(file):
  if not hasattr(file, "read"):
    handle = open(file)
  else:
    handle = file

  while True:
    line = handle.readline()
    if not line.endswith('\n'):
      break
    else:
      yield line
  handle.close()

def get_accesion(file, lineWrapping):
  result = []
  accession_pattern = re.compile('([NX]M_\d+\.*\d*)\.')
  for line in line_yielder(file):
#    pdb.set_trace()
    m = accession_pattern.search(line)
    if m:
      result.append(m.group(1))
  if lineWrapping: 
    for i in range(0, len(result), 3):
      print " ".join(result[i:i+3])
  else:
    print " ".join(result)
if __name__ == '__main__': 
  o = sys.stdout
  e = sys.stderr
  parser= argparse.ArgumentParser(description="")
  parser.add_argument("file", help="")
  parser.add_argument('--lineWrapping', action='store_true', 
  help="print three accessions in a single line. default=[False]", default=False)
  args = parser.parse_args()
  
  get_accesion(args.file, args.lineWrapping)
