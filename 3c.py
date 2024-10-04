def read_intervals(file):
    intervals = []
    with open(file, 'r') as f:
        for line in f:
            start, end = map(int, line.strip().split(','))
            intervals.append((start, end))
    return intervals

def read_fasta(file):
    sequence = ""
    with open(file, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence

def is_within_interval(index, intervals):
    for start, end in intervals:
        if start <= index <= end:
            return True
    return False

def count_nucleotides(sequence, intervals):
    g_h, g_l, a_h, a_l = 0, 0, 0, 0

    for i, nucleotide in enumerate(sequence):
        if is_within_interval(i + 1, intervals):
            if nucleotide in 'GC':
                g_h += 1
            elif nucleotide in 'AT':
                a_h += 1
        else:
            if nucleotide in 'GC':
                g_l += 1
            elif nucleotide in 'AT':
                a_l += 1

    return g_h, g_l, a_h, a_l

interval_file = 'test-2a/my_out_2a/viterbi-intervals.txt' 
fasta_file = 'hmm-sequence.fasta'

intervals = read_intervals(interval_file)
sequence = read_fasta(fasta_file)

g_h, g_l, a_h, a_l = count_nucleotides(sequence, intervals)

print(f"G or C within intervals (g_h): {g_h}")
print(f"G or C outside intervals (g_l): {g_l}")
print(f"A or T within intervals (a_h): {a_h}")
print(f"A or T outside intervals (a_l): {a_l}")


# checking my previous implementation with a simpler way/hard coding it
def check(seq, start, end):
    g, a = 0, 0
    for j, nucleotide in enumerate(sequence):
        i = j + 1
        if ((i>=start) and (i<=end)):
            if nucleotide in 'GC':
                g += 1
            if nucleotide in 'AT':
                a += 1
    # print(f'G/C: {g}, A/T: {a}')
    return g, a
                
# high
# (66,415)
# (527,720)
# (951,1000)

# low 
# (1, 65)
# (416, 526)
# (721, 950)
  
# print('high!')
g1, a1 = check(sequence, 66, 415)
g2, a2 = check(sequence, 527, 720)
g3, a3 = check(sequence, 951, 1000)

high_g = g1 + g2 + g3
high_a = a1 + a2 + a3

# print('\nlow!')
c1, t1 = check(sequence, 1, 65)
c2, t2 = check(sequence, 416, 526)
c3, t3 =check(sequence, 721, 950)

low_g = c1 + c2 + c3
low_a = t1 + t2 + t3

print(f'high_g: {high_g}')
print(f'low_g: {low_g}')
print(f'high_a: {high_a}')
print(f'low_a: {low_a}')

assert high_g == g_h, 'not the same G/C numbers for high intervals'
assert high_a == a_h, 'not the same A/T numbers for high intervals'
assert low_g == g_l, 'not the same G/C numbers for low intervals'
assert low_a == a_l, 'not the same G/C numbers for low intervals'
    
    
    
