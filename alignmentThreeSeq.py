#!/usr/bin/env python2.7
import sys
import pdb
import argparse
import random 
import numpy as np
import copy
from Bio import SeqIO
from itertools import permutations, combinations, repeat
from StringIO import StringIO
from guppy import hpy
import time
import multiprocessing

start = time.time()
h = hpy()

mismatch = -4
indel    = -8
match    = 5


def Prep_Score_Hash(mismatch=-4, indel=-8, match=5):
  ScoreTwoChar = {}
  #store all possible mismatch and indel scores into ScoreTwoChar Hash table
  alphabet = ['A', 'C', 'G', 'T', '-']
  for a in alphabet:
    for b in alphabet:
      if a == b :
        if not a == '-':
          ScoreTwoChar[(a, b)] = match
        else:
          ScoreTwoChar[(a, b)] = 0
      else:
        if a == '-' or b == '-':
          ScoreTwoChar[(a, b)] = indel
        else:
          ScoreTwoChar[(a, b)] = mismatch
  ScoreThreeChar = {}
  for a in alphabet:
    for b in alphabet:
      for c in alphabet:
        if a == b and b == c and c == '-':
          continue
        else:
          ScoreThreeChar[(a, b, c)] = sum([ ScoreTwoChar[i] for i in combinations([a, b, c], 2)])
  return ScoreThreeChar

score = Prep_Score_Hash()

def printScoreMatrix():
  Score = Prep_Score_Hash()
  for k, v in sorted(Score.items(), key=lambda x: x[0]):
    print "%r => %r " % (k, v)

def scoreTwoSeq(x, y):
  """This function returns a 2D score matrix.
  even though the name of function says this is a two sequence scoring function
  actually in this function  the score is caluclated as if calculating three sequences
  """
  global score

  new_indel = score[('A', '-', '-')]
  new_match = score[('A', 'A', '-')]
  new_mismatch = score[('A', 'C', '-')]

  M = np.zeros((len(x) + 1, len(y) + 1))
  Path = np.empty((len(x) + 1, len(y) + 1), dtype=object)
  for i in range(1, len(y) + 1):
    M[0][i] = i * new_indel
    Path[0][i] = (0, i-1)
  for j in range(1, len(x) + 1):
    M[j][0] = j * new_indel 
    Path[j][0] = (j-1, 0)
  
  for i in range(1, len(x) + 1):
    for j in range(1, len(y) + 1):
      if x[i-1] == y[j-1]:
        M[i][j] = max(M[i-1][j-1] + new_match, M[i-1][j] + new_indel, M[i][j-1] + new_indel)
        if M[i][j] == M[i-1][j-1] + new_match:
          Path[i][j] =  (i-1, j-1)
        elif M[i][j] == M[i-1][j] + new_indel:
          Path[i][j] = (i-1, j)
        else:
          Path[i][j] = (i, j-1)
      else:
        M[i][j] = max(M[i-1][j-1] + new_mismatch, M[i-1][j] + new_indel, M[i][j-1] + new_indel)
        if M[i][j] == M[i-1][j-1] + new_mismatch:
          Path[i][j] =  (i-1, j-1)
        elif M[i][j] == M[i-1][j] + new_indel:
          Path[i][j] = (i-1, j)
        else:
          Path[i][j] = (i, j-1)
# backtracing
#  row = []
#  column= []
#  middle = []
#  i = len(x)
#  j = len(y)
#  while Path[i][j]:
#    new_i, new_j = Path[i][j]
#    if i - new_i > 0 and j - new_j > 0 :
#      column.insert(0, x[i-1])
#      row.insert(0, y[j-1])
#      if x[i-1] == y[j-1]:
#        middle.insert(0, '|')
#      else:
#        middle.insert(0, ':') 
#    elif j - new_j > 0:
#      row.insert(0, y[j-1])
#      column.insert(0, '-')
#      middle.insert(0, 'x') 
#    elif i - new_i > 0:
#      column.insert(0, x[i-1])
#      row.insert(0, '-')
#      middle.insert(0, 'x') 
#    i = new_i
#    j = new_j

#  return row, column, middle
  return M, Path

def scoreThreeSeq(A, B, C):
  """
  This function returns the last surface of the alignment 
  between three sequences.
  """
  global score
  new_indel = score[('A', '-', '-')]
  prev, Path = scoreTwoSeq(B, C)
  
#  pdb.set_trace()
  new = np.zeros((len(B) + 1, len(C) + 1)) 
  
  for a in range(1, len(A) + 1):
    new[0, 0] += new_indel

  # fill in the dumpy row when B == 0
    for c in range(1, len(C) + 1):
      new[0, c] = max(new[0, c-1] + new_indel, prev[0, c] + new_indel, 
          prev[0, c-1] + score[(A[a-1], C[c-1], '-')])
  # fill in the dumpy row when C == 0
    for b in range(1, len(B) + 1):
      new[b, 0] = max(new[b-1, 0] + new_indel, prev[b, 0] + new_indel, 
          prev[b-1, 0] + score[(A[a-1], B[b-1], '-')])
     
    for b in range(1, len(B) + 1):
      for c in range(1, len(C) + 1):
        #1-3 edges, 4-6 diagnoal edges, 7, cube diagnoal edge
        seven = prev[b-1, c-1] + score[(A[a-1], B[b-1], C[c-1])]
        six   = prev[b-1, c] + score[(A[a-1], B[b-1], '-')]
        five  = prev[b, c-1] + score[(A[a-1], '-', C[c-1])]
        three = prev[b, c] + score[(A[a-1], '-', '-')]
        four  = new[b-1, c-1] + score[('-', B[b-1], C[c-1])]
        two   = new[b-1, c] + score[('-', B[b-1], '-')]
        one   = new[b, c-1] + score[('-', '-', C[c-1])]
        new[b, c] =  max([one, two, three, four, five, six, seven])
    prev = copy.deepcopy(new)
  return new


def alignThreeSeq(A, B, C):
  global score

  M = np.zeros((len(A) + 1, len(B) + 1, len(C) + 1))
#  np.ndarray.fill(M[0, :, :] , -100)

  P = np.empty((len(A) + 1, len(B) + 1, len(C) + 1), dtype=object) 
# A number is assigned to P[a, b, c] to traceback.
# please take a look at the fig here 3DAlign.pdf if you don't know which corner 
# of the cube does a specific number refers to 
  ################################################
  # initialization

# initialize three dummy plane
  M[0, :, : ], P_tmp = scoreTwoSeq(B, C)
  # below magic line converts an two dimentional tuple into a three dimensional tuple
  # without breaking the shape of original array
  P[0, :, :] = np.array([(0, x[0], x[1]) if x else x for x in P_tmp.flatten()]).reshape(P_tmp.shape)
  M[:, 0, : ], P_tmp = scoreTwoSeq(A, C)
  P[:, 0, :]= np.array([(x[0], 0, x[1]) if x else x for x in P_tmp.flatten()]).reshape(P_tmp.shape)
  M[:, :, 0 ], P_tmp = scoreTwoSeq(A, B)
  P[:, :, 0] = np.array([(x[0], x[1], 0) if x else x for x in P_tmp.flatten()]).reshape(P_tmp.shape)
#  num_exact_equal = 0 # stores the number of postions, in which nucleotides are identical
  for a in range(1, len(A) + 1):
    for b in range(1, len(B) + 1):
      for c in range(1, len(C) + 1):
#        if a == 1 and b == 1 and c == 1: continue
        seven = M[a-1, b-1, c-1] + score[(A[a-1], B[b-1], C[c-1])]
        six   = M[a-1, b-1, c] + score[(A[a-1], B[b-1], '-')]
        five  = M[a-1, b, c-1] + score[(A[a-1], '-', C[c-1])]
        three = M[a-1, b, c] + score[(A[a-1], '-', '-')]
        four  = M[a, b-1, c-1] + score[('-', B[b-1], C[c-1])]
        two   = M[a, b-1, c] + score[('-', B[b-1], '-')]
        one   = M[a, b, c-1] + score[('-', '-', C[c-1])]
        M[a, b, c] =  max([one, two, three, four, five, six, seven])

        if M[a, b, c] == seven:
          P[a, b, c] = (a-1, b-1, c-1)
        elif M[a, b, c] == six:
          P[a, b, c] = (a-1, b-1, c)
        elif M[a, b, c] == four:
          P[a, b, c] = (a, b-1, c-1)
        elif M[a, b, c] == five:
          P[a, b, c] = (a-1, b, c-1)
        elif M[a, b, c] == one:
          P[a, b, c] = (a, b, c-1)
        elif M[a, b, c] == two:
          P[a, b, c] = (a, b-1, c)
        elif M[a, b, c] == three:
          P[a, b, c] = (a-1, b, c)
        else:
          raise ValueError
#  pdb.set_trace()
  #Get the alignment str
  A_aln, B_aln, C_aln = [], [], []
#  if trace_alignment:
  a, b, c = len(A), len(B), len(C)
  while P[a, b, c]:
    new_a, new_b, new_c = P[a, b, c]
    if a - new_a > 0:
      A_aln.insert(0, A[a-1])
    else:
      A_aln.insert(0, '-')

    if b - new_b > 0:
      B_aln.insert(0, B[b-1])
    else:
      B_aln.insert(0, '-')

    if c - new_c > 0:
      C_aln.insert(0, C[c-1])
    else:
      C_aln.insert(0, '-')
    a, b, c = new_a, new_b, new_c
  max_score = M[len(A), len(B), len(C)]
#  return A_aln, B_aln, C_aln, max_score, num_exact_equal
  return A_aln, B_aln, C_aln, max_score

def paritionBC(M_upper, M_down):
  # see post: http://stackoverflow.com/questions/42519/how-do-you-rotate-a-two-dimensional-array
  reversed_M_down = np.fliplr(np.transpose(np.fliplr(np.transpose(M_down))))
  sum_array = M_upper + reversed_M_down
  index = np.unravel_index(sum_array.argmax(), sum_array.shape)
  return index

def recursive_call(A, B, C):
  if len(A) <= 1 or len(B) <= 1 or len(C) <= 1:
#    upper, middle, down, score, num_exact_equal = alignThreeSeq(A, B, C, trace_alignment)
    upper, middle, down, score = alignThreeSeq(A, B, C)
  else:
    xmid = len(A)/2
    matrix_upper = scoreThreeSeq(A[:xmid], B, C)
    matrix_down = scoreThreeSeq(A[xmid:][::-1], B[::-1], C[::-1])
    b, c = paritionBC(matrix_upper, matrix_down)

#    upper_l, middle_l, down_l, score_l, num_l = recursive_call(A[:xmid], B[:b], C[:c], trace_alignment)
#    upper_r, middle_r, down_r, score_r, num_r = recursive_call(A[xmid:], B[b:], C[c:], trace_alignment)
    upper_l, middle_l, down_l, score_l = recursive_call(A[:xmid], B[:b], C[:c])
    upper_r, middle_r, down_r, score_r = recursive_call(A[xmid:], B[b:], C[c:])

    upper = upper_l + upper_r 
    middle = middle_l + middle_r
    down = down_l + down_r
    score = score_l + score_r
#    num_exact_equal = num_l + num_r

#  return upper, middle, down, score, num_exact_equal
  return upper, middle, down, score

def testScoreTwoSeq():
  Xs = ["AGTACGCA", "hello", "T", "T", "T"]
  Ys = ["TATGC", "hllo", "C", "T", ""]

  for x, y in zip(Xs, Ys):
    u, d, m = scoreTwoSeq(x, y) 
    print u
    print m
    print d

def testAlignThreeSeq():
  for x, y, z in zip( ( 'GC', 'GC', 'C'), ('C', 'GC', 'GC' ), ('C', 'GC', 'GC' )):
    u, m, d, s = alignThreeSeq(x, y, z)
    print "score: %r" % s
    print u
    print m
    print d

def get_seq(file):
  result = ""
  for seq in SeqIO.parse(file, "fasta"):
    result += str(seq.seq)
  return result

def test_recursive_call():
  As = ['AGC' , 'AG' , 'AG' , 'AGT']
  Bs = ['AGC' , 'CG' , 'CG' , 'CG']
  Cs = ['AC'  , 'TG' , 'TC' , 'TC']

  for i, (A, B, C) in enumerate(zip(As, Bs, Cs)):
    u, m, d, s = recursive_call(A, B, C)
    print "#" * 8, "alignment %r" % i, '#' * 8
    print "score: %r" % s
    print u
    print m
    print d

def write_align(stream, u, m, d, match_str, width=80):
  """
  u, m, d are equal length string
  """
  for i in range(0, len(u), width):
    for label, j in zip(['A', 'B', 'C', 'M'], [u, m, d, match_str]):
      stream.write("%r:  %s\n" % (label, j[i:i+width]))
    stream.write("\n") 

def getJobDone(fileA, fileB, fileC):
  out_name = "_".join(map(lambda x: x.split('.')[0], [fileA, fileB, fileC]))
  print "In process: %r" % out_name
  A = get_seq(fileA)
  B = get_seq(fileB)
  C = get_seq(fileC)  
  with open(out_name, 'w') as OUT:
#      u, m, d, s, num_exact_equal = recursive_call(A, B, C, args.outputAlign)
    u, m, d, s = recursive_call(A, B, C)
#      if args.outputAlign:
#    pdb.set_trace()
    num_exact_equal= [1 if a == b and b == c else 0 for a, b, c in zip(u, m, d)]
    tmp_list = "".join(map(str, num_exact_equal))
    num_exact_equal = sum(num_exact_equal)
    tmp_list = tmp_list.replace('1','*')
    tmp_list = tmp_list.replace('0', ' ')
    u = "".join(u)
    m = "".join(m)
    d = "".join(d)
#    OUT.write("#" * 8 + " alignment %r " % i + '#' * 8 + "\n")
    OUT.write("score: %r\tlength:%r\tperfect match nts: %r\n" % (s, len(u), num_exact_equal))
    write_align(OUT, u, m, d, tmp_list)
    OUT.write("Memory_usage: %r (bytes)\t Running_time: %.2f (s)\n" % (h.heap().size, time.time() - start))

if __name__ == '__main__': 
  o = sys.stdout
  e = sys.stderr
  parser= argparse.ArgumentParser(description="")
  parser.add_argument("file", help="", nargs = '+')
  parser.add_argument("--outputAlign", 
      help="logical if given alignment be write to file [Default: False]", action= 'store_true' )
  args = parser.parse_args()
  printScoreMatrix()
#  testAlignThreeSeq()
  
##  testScoreTwoSeq()
#  if not len(args.file) % 3 == 0:
#    print "Error! number of input files have to be multiples of three"
##    pdb.set_trace()
#    sys.exit(1)
#    
#  getJobDone(args.file[0], args.file[1], args.file[2])
##  for i in range(0, len(args.file), 3):
###    A = get_seq(args.file[i])
###    B = get_seq(args.file[i + 1])
###    C = get_seq(args.file[i + 2])  
##    p = multiprocessing.Process(target=getJobDone, args=(args.file[i], args.file[i+1], args.file[i+2]))
##    p.start()
#
