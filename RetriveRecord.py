#!/usr/bin/env python2.7
import sys
import pdb
import argparse

from Bio import Entrez

def retrieve(id, format):
  Entrez.email = "Your.Name.Here@example.org"
  handle = Entrez.efetch(db="nucleotide", id=id, rettype=format, retmode="text")
  suffix = "." + format
  with open(id+suffix, 'w') as OUT:
    for line in handle.readlines():
      OUT.write(line)

if __name__ == '__main__':
  o = sys.stdout
  e = sys.stderr
  parser= argparse.ArgumentParser(description="")
  parser.add_argument("id", nargs = '+', help="the geneid list")
  parser.add_argument("--format", help="format default=[Fasta]", default="fasta")
  args = parser.parse_args() 

  for id in args.id:
    retrieve(id, args.format)
