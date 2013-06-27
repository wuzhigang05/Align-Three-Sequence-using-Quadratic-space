Align three sequence using Quadratic space
===========================================

Space, not running time, is often the factor limiting the application of dynamic programming. 
Often cubic space is required to align three sequences (3D dynamic programming). Given this 
space requirement, alignment of three genes with a length of more than 1000 characters is 
almost infeasible in most of our modern computers. Here, by using a divide conque technique, 
I implemented alignment of three sequences in quadratic space. The idea of this implementation 
is similar to [Hirschberg\'s algorithm](http://en.wikipedia.org/wiki/Hirschberg's_algorithm), which 
I have an implmentation in Python, see https://github.com/wuzhigang05/Dynamic-Programming-Linear-Space.git

Requirements
=============
  1. Python version >2.7
  
  2. package: [numpy] (http://www.numpy.org/)

  3. package: [Biopython] (http://biopython.org/wiki/Main_Page)

Example Usage
=============
Below command will write the alignment of three sequences stored in file test1, test2 and test3, respectively
   to a file named with test1_test2_test3 to the current directory. In addition to the alignment, the file also 
   contains the memory and time usage information. 

    python alignmentThreeSeq.py test1 test2 test3

To display help

    python AlignTwoStringInLinearSpace.py -h

Send Bugs/Commnents to
======================
Zhigang Wu (zhigang.wu@email.ucr.edu)



LICENSE
=========
Copyright (c) <2013>, <Zhigang Wu>
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, 
       this list of conditions and the following disclaimer.
    
    2. Redistributions in binary form must reproduce the above copyright 
       notice, this list of conditions and the following disclaimer in the
       documentation and/or other materials provided with the distribution.

    3. Neither the name of Django nor the names of its contributors may be used
       to endorse or promote products derived from this software without
       specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

