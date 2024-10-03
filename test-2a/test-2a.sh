#!/usr/bin/env bash

# To run the test, make sure that your script (`2a.py`) 
# and the files needed (`test-2a-1.fasta`, `test-2a-2.fasta`, `test-2a-3.fasta`) are in the working directory.
# Or modify the commands below to accomodate for their path

# In the terminal, execute `bash test-2a/test-2a.sh` to run all the test cases.
# The folder `out_2a` contains the correct output for all test cases
# Correct Viterbi probabilities of all these test cases are written in order in the file `test-2b-out.txt`.


# TEST 1
python 2a.py -f test-2a/test-2a-1.fasta -mu 0.01 -out test-2a/my_out_2a/viterbi-intervals-1.txt
echo -e "----------------------------------------------\n"
# The correct viterbi interval of this test case is in `viterbi-intervals-1.txt`
# The correct log probability of each state at each position of this test case is in `m-1.json` (THE CONTENT OF THIS FILE IS NOT GRADED, JUST PROVIDED AS A DEBUGGING GUIDE)

# TEST 2
python 2a.py -f test-2a/test-2a-1.fasta -mu 0.5 -out test-2a/my_out_2a/viterbi-intervals-2.txt
echo -e "----------------------------------------------\n"
# The correct viterbi interval of this test case is in `viterbi-intervals-2.txt`
# The correct log probability of each state at each position of this test case is in `m-2.json` (THE CONTENT OF THIS FILE IS NOT GRADED, JUST PROVIDED AS A DEBUGGING GUIDE)

# TEST 3
python 2a.py -f test-2a/test-2a-2.fasta -mu 0.6 -out test-2a/my_out_2a/viterbi-intervals-3.txt
echo -e "----------------------------------------------\n"
# The correct viterbi interval of this test case is in `viterbi-intervals-3.txt`

# TEST 4
python 2a.py -f test-2a/test-2a-2.fasta -mu 0.08 -out test-2a/my_out_2a/viterbi-intervals-4.txt
echo -e "----------------------------------------------\n"
# The correct viterbi interval of this test case is in `viterbi-intervals-4.txt`

# TEST 5
python 2a.py -f test-2a/test-2a-3.fasta -mu 0.5 -out test-2a/my_out_2a/viterbi-intervals-5.txt
# The correct viterbi interval of this test case is in `viterbi-intervals-5.txt`




