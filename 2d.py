#!/usr/bin/env python3

'''Script for computing marginalized log-likelihoods.
Arguments:
    -f: file containing the sequence (fasta file)
    -out: file to output log-likelihoods

Outputs:
    - CSV with log-likelihoods
    - log-likelihoods plot (2c.png)
    
Example Usage:
    python 2d.py -f hmm-sequence.fasta -out ll.csv
'''

import argparse
import numpy as np
import math
import matplotlib.pyplot as plt


'''Computes the log(exp(a) + exp(b)) with better numerical stability'''
def sumLogProbs(a,b):
  if a > b:
    return a + np.log1p(math.exp(b-a))
  else:
    return b + np.log1p(math.exp(a-b))
    

'''Reads the fasta file and outputs the sequence to analyze.
Arguments:
	filename: name of the fasta file
Returns:
	s: string with relevant sequence
'''
def read_fasta(filename):
    with open(filename, "r") as f:
        s = ""
        for l in f.readlines()[1:]:
            s += l.strip()
    return s


''' Outputs the forward and backward probabilities of a given observation.
'''
def forward_backward(obs, trans_probs, emiss_probs, init_probs):
    ''' The same as your implementation from 2b.
        Copy other necessary helper functions you created for 2b to 2d.py.
    '''
    pass
    
    
''' Computes the log-likelihood of the FASTA observation given mu
Arguments:
    mu: the parameter given for computing the log-likelihood of the observation.
Returns:
    the log-likelihood
'''
def log_likelihood(fasta_file, mu):
    ''' Complete this function. '''
    pass


def main():
    parser = argparse.ArgumentParser(
        description='Compute marginalized log-likelihoods for a given sequence.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-out', action="store", dest="out", type=str, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    ll_file = args.out

    a = np.arange(0.001, 0.1, 0.001)
    b = np.arange(0.1, 0.25, 0.005)
    x = np.concatenate((a, b))
    y = np.array([log_likelihood(fasta_file, mu) for mu in x])
    np.savetxt(ll_file, np.vstack((x, y)), delimiter=",")

    ll = np.loadtxt(ll_file, delimiter = ",")
    ll_x = ll[0, :]
    ll_y = ll[1, :]
    plt.plot(ll_x, ll_y, 'r-')
    plt.title('mu vs. log(P(x | mu))')
    plt.savefig('2c.png')
    
    mu_mle = ll_x[np.argmax(ll_y)]
    print("Approximate MLE of mu: {:.3f}".format(mu_mle))

if __name__ == "__main__":
    main()
