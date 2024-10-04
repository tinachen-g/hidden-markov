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
    N = len(obs)                                       # Number of emissions (N)
    K = len(init_probs)                                # Number of hidden states (K)
    F = {state: np.full(N, None) for state in init_probs}   # Forward log-probabilities
    B = {state: np.full(N, None) for state in init_probs}   # Backward log-probabilities
    R = {state: np.full(N, None) for state in init_probs}   # Posterior probabilities (non-log scale)

    # forward initialization
    for state in init_probs:
        F[state][0] = init_probs[state] + emiss_probs[state][obs[0]]
    
    # forward iteration
    for t in range(1, N):
        for state in init_probs:
            x_t = obs[t]
            # summation = 0
            summation = -np.inf
            for i in init_probs:
                summation = sumLogProbs(summation, F[i][t-1] + trans_probs[i][state])
                # summation = summation + F[i][t-1] + trans_probs[i][state]
                # summation = sumLogProbs(summation, F[i][t-1] + trans_probs[i][state])
            F[state][t] = summation + emiss_probs[state][x_t]

    # like_f = 0
    like_f = -np.inf
    for state in init_probs:
        # like_f = like_f + F[state][N-1]
        like_f = sumLogProbs(like_f, F[state][N-1])

    # bckward initializae
    for state in init_probs:
        # log(1) = 0
        B[state][N-1] = 0

    # backward iteration
    for t in range(N-2, -1, -1):
        for state in init_probs:
            summation = -np.inf 
            for i in init_probs:
                summation = sumLogProbs(summation, B[i][t+1] + trans_probs[state][i] + emiss_probs[i][obs[t+1]])
            B[state][t] = summation
    
    # like_b = 0
    like_b = -np.inf
    for state in init_probs:
        # like_b = like_b + init_probs[state] + emiss_probs[state][obs[0]] + B[state][0]
        like_b = sumLogProbs(like_b, init_probs[state] + emiss_probs[state][obs[0]] + B[state][0])


    # posterior
    assert np.isclose(like_f, like_b,)
    assert np.isclose(math.e**like_f, math.e**like_b,)

    for t in range(N):
        local_estimate = -np.inf
        for k in init_probs:
            local_estimate = sumLogProbs(local_estimate, F[k][t] + B[k][t])
        
        local_estimate_exp= np.exp(local_estimate)
        assert np.isclose(math.exp(like_b), local_estimate_exp)
        assert np.isclose(math.exp(like_f), local_estimate_exp)

        for k in init_probs:
            R[k][t] = np.exp(F[k][t] + B[k][t] - local_estimate) 

    return F, like_f, B, like_b, R
    
    
''' Computes the log-likelihood of the FASTA observation given mu
Arguments:
    mu: the parameter given for computing the log-likelihood of the observation.
Returns:
    the log-likelihood
'''
def log_likelihood(fasta_file, mu):
    ''' Complete this function. '''
    obs = read_fasta(fasta_file)

    init_probs = {'h': np.log(0.5), 'l': np.log(0.5)}
    
    trans_probs = {
        'h': {'h': np.log(1 - mu), 'l': np.log(mu)},
        'l': {'h': np.log(mu), 'l': np.log(1 - mu)}
    }

    emiss_probs = {
        'h': {'A': np.log(0.13), 'C': np.log(0.37), 'G': np.log(0.37), 'T': np.log(0.13)},
        'l': {'A': np.log(0.32), 'C': np.log(0.18), 'G': np.log(0.18), 'T': np.log(0.32)}
    }

    _, like_f, _, _, _ = forward_backward(obs, trans_probs, emiss_probs, init_probs)
    
    return like_f


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
