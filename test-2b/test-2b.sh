#!/usr/bin/env bash

# Parameters to execute 2b.py for the test output:
# test-2b-1
echo "**** test-2b-1 ******************************"
python 2b.py -f test-2b/test-2b-1.fasta -mu 0.01 

# test-2b-2
echo -e "\n"
echo "**** test-2b-2 ******************************"
python 2b.py -f test-2b/test-2b-2.fasta -mu 0.01 