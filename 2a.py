#!/usr/bin/env python3

'''Script for computing GC-rich and GC-poor intervals in a given sequence.
Arguments:
    -f: file containing the sequence (fasta file)
    -mu: the probability of switching states
    -out: file to output intervals to (1 interval per line)

Outputs:
    File with list of intervals (a_i, b_i) such that bases a_i to b_i are
    classified as GC-rich.
    
Example Usage:
    python 2a.py -f hmm-sequence.fasta -mu 0.01 -out test-2a/my_out_2a/viterbi-intervals.txt
'''

import argparse
import numpy as np

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


''' Outputs the Viterbi decoding of a given observation.
Arguments:
	obs: observed sequence of emitted states (list of emissions)
	trans_probs: transition log-probabilities (dictionary of dictionaries)
	emiss_probs: emission log-probabilities (dictionary of dictionaries)
	init_probs: initial log-probabilities for each hidden state (dictionary)
Returns:
	l: list of most likely hidden states at each position
        (list of hidden states)
	p: log-probability of the returned hidden state sequence
'''
def viterbi(obs, trans_probs, emiss_probs, init_probs):
    ''' Complete this function. '''

    # m is a dictionary that has two keys, 'h' and 'l', the values of each key being a list with the length of the observations.
    # The value list for key 'h' should be storing the log-probabilities of the optimal Viterbi path to state 'h' up to observation x_t for each t.
    # The value list for key 'l' should be storing the log-probabilities of the optimal Viterbi path to state 'l' up to observation x_t for each t.
    # We provided the correct m in this format (in JSON format) for test case 1 & 2 to help you debug, but m will not be part of the grade.
    # which means that you don't have to use m as specified here if you find it more natural using other implementations.
    m = {
            'h' : [None for i in range(len(obs))],
            'l' : [None for i in range(len(obs))],
        }
    
    t = len(obs)
    # initialization
    for states in m:
        m[states][0] = init_probs[states] + emiss_probs[states][obs[0]]
        
    tb = np.full((2, t+1,), -1)
    H, L = 0, 1
    for i in range(1, t):
        x_t = obs[i]
        v_h_from_h = m['h'][i-1] + emiss_probs['h'][x_t] + trans_probs['h']['h']
        v_l_from_h = m['l'][i-1] + emiss_probs['h'][x_t] + trans_probs['l']['h']
        v_max_from_h = max(v_h_from_h, v_l_from_h)
        m['h'][i] = v_max_from_h

        tb[H][i] = H if np.isclose(v_h_from_h, v_max_from_h) else L

        v_h_from_l = m['h'][i-1] + emiss_probs['l'][x_t] + trans_probs['h']['l']
        v_l_from_l = m['l'][i-1] + emiss_probs['l'][x_t] + trans_probs['l']['l']
        v_max_from_l = max(v_h_from_l, v_l_from_l)
        m['l'][i] = v_max_from_l

        tb[L][i] = H if (v_h_from_l == v_max_from_l) else L

    path_prob = max(m['h'][-1], m['l'][-1])

    curr = H if (path_prob == m['h'][-1]) else L
    path = ['h' if curr == H else 'l']
    time = t-1
    while (time > 0): 
        curr = tb[curr][time]
        path.insert(0, 'h' if curr == H else 'l')
        time -= 1

    return path, path_prob


''' Returns a list of non-overlapping intervals describing the GC rich regions.
Arguments:
	sequence: list of hidden states
Returns:
	intervals: list of tuples (i, j), 1 <= i <= j <= len(sequence), that
                describe GC rich regions in the input list of hidden states.
'''
def find_intervals(sequence):
    '''Finds contiguous intervals of 'h' in the sequence.'''
    intervals = []
    i = 0
    while i < len(sequence):
        if sequence[i] == 'h':
            j = i
            while j < len(sequence) and sequence[j] == 'h':
                j += 1
            intervals.append((i + 1, j))
            i = j
        else:
            i += 1

    return intervals



def main():
    parser = argparse.ArgumentParser(
        description='Parse a sequence into GC-rich and GC-poor regions using Viterbi.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-mu', action="store", dest="mu", type=float, required=True)
    parser.add_argument('-out', action="store", dest="out", type=str, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    mu = args.mu
    intervals_file = args.out

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
    sequence, p = viterbi(obs_sequence, transition_probabilities,
                        emission_probabilities, initial_probabilities)
    intervals = find_intervals(sequence)
    with open(intervals_file, "w") as f:
        f.write("\n".join([("%d,%d" % (start, end)) for (start, end) in intervals]))
        f.write("\n")
    print("Viterbi probability in log scale: {:.2f}".format(p))


if __name__ == "__main__":
    main()
