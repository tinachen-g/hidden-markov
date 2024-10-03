#!/usr/bin/env python3

'''Script for computing posterior probabilities of hidden states at each
   position of a given sequence.
Arguments:
    -f: file containing the sequence (fasta file)
    -mu: the probability of switching states

Outputs:
    posteriors.csv - a KxL matrix outputted as a CSV with the posterior
                     probability of each state at each position

Example Usage:
    python 2b.py -f hmm-sequence.fasta -mu 0.01
'''

import argparse
import numpy as np
import math


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
Arguments:
    obs (string): A sequence of observed emissions.
    trans_probs (dict of dicts): Transition log-probabilities between hidden states. 
        The keys of the outer dictionary represent the current hidden state, and the inner keys represent the next hidden state.
        For example, `trans_probs[i][j]` gives the log-probability of transitioning from state `i` to state `j`.
    emiss_probs (dict of dicts): Emission log-probabilities for each hidden state. 
        The outer keys are hidden states, and the inner keys are observed states (emissions). 
        For example, `emiss_probs[j][k]` gives the log-probability of emitting state `k` from hidden state `j`.
    init_probs (dict): Initial log-probabilities for each hidden state.

Returns:
    F (dict): A dictionary of forward log-probabilities. Each key is a hidden state, and each value is a NumPy array of shape `(N)`, where `N` is the length of the observed sequence.
    like_f (float): The log-likelihood of the observation sequence as computed by the forward algorithm, i.e., `log P(obs)`.
    B (dict): A dictionary of backward log-probabilities. Each key is a hidden state, and each value is a NumPy array of shape `(N)`.
    like_b (float): The log-likelihood of the observation sequence as computed by the backward algorithm, i.e., `log P(obs)`.
    R (dict): A dictionary of posterior probabilities (in normal scale, not log-scale). Each key is a hidden state, and each value is a NumPy array of shape `(N)`.

Requirements:
    - F (forward log-probabilities), B (backward log-probabilities), and R (posterior probabilities) are dictionaries with keys as hidden states and values as numpy arrays of shape (N), where N represent the number of emissions (i.e., the length of the observed sequence).
'''
def forward_backward(obs, trans_probs, emiss_probs, init_probs):
    ''' Complete this function. 
    '''
    # Please follow the data format design below
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


def main():
    parser = argparse.ArgumentParser(
        description='Compute posterior probabilities at each position of a given sequence.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-mu', action="store", dest="mu", type=float, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    mu = args.mu

    obs_sequence = read_fasta(fasta_file)
    transition_probabilities = {
        'h': {'h': np.log(1 - mu), 'l': np.log(mu)},
        'l': {'h': np.log(mu), 'l': np.log(1 - mu)}
    }
    emission_probabilities = {
        'h': {'A': np.log(0.13), 'C': np.log(0.37), 'G': np.log(0.37), 'T': np.log(0.13)},
        'l': {'A': np.log(0.32), 'C': np.log(0.18), 'G': np.log(0.18), 'T': np.log(0.32)}
    }
    initial_probabilities = {'h': np.log(0.5), 'l': np.log(0.5)}
    F, like_f, B, like_b, R = forward_backward(obs_sequence,
                                              transition_probabilities,
                                              emission_probabilities,
                                              initial_probabilities)
    
    R_arr = np.array([R[state] for state in ['h', 'l']])
    np.savetxt("posteriors.csv", R_arr, delimiter=",", fmt='%.4e')
    print("Backward log-likelihood: {:.2f}".format(like_b))
    print("Forward log-likelihood: {:.2f}".format(like_f))


if __name__ == "__main__":
    main()
