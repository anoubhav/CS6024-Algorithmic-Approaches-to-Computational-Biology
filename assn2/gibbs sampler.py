import sys
import os
import random
from bisect import bisect
from time import time
def lin_cong_generator(X = 1,a = 1140671485, c = 128201163, m = 2**24, num=1):
    """ linear congruential generator i.e., a psuedo random number generator.
    Input: X0, a, c, m
    Output: list of Xn, last Xn """
    rnd_nums = []
    for _ in range(num):
        Xn = (a*X+c)%m
        rnd_nums.append(Xn)
        X = Xn
    return rnd_nums, rnd_nums[-1]

def read_fasta(FNAME):
    sequences = []
    with open(FNAME) as file:
        for line in file.readlines():
            if '>' in line: continue
            sequences.append(line.strip().upper())
    return sequences
    
def get_count_matrix(motifs):
    """ Returns the count matrix given a collection of motifs """
    # Initialize the count matrix. Numpy matrix is preferred due to speed of execution.
    count_matrix = [[0 for x in range(len(motifs[0]))] for x in range(4)]
    for col_ind, col in enumerate(zip(*motifs)):
        count_matrix[0][col_ind] = col.count('A')
        count_matrix[1][col_ind] = col.count('T')
        count_matrix[2][col_ind] = col.count('G')
        count_matrix[3][col_ind] = col.count('C')
    return count_matrix

def get_score(motifs):
    """ Obtains the score of the motifs which is the sum of the (hamming) distance between the motifs and the consensus motif """

    score = 0
    count_matrix = get_count_matrix(motifs)
    k = len(count_matrix[0])
    t = len(motifs)
    for i in range(k):
        num_deviations = t - max(list(zip(*count_matrix))[i])  # number of deviations from consensus
        score += num_deviations

    return score
    
def get_profile(motifs):
    """ Returns the profile matrix (after pseudocount) of a collection of motifs """
    count_matrix = get_count_matrix(motifs)

    PSEUDOCOUNT = 1
    for row in range(len(count_matrix)):
        for col in range(len(count_matrix[0])):
            count_matrix[row][col] += PSEUDOCOUNT
            count_matrix[row][col] = count_matrix[row][col]/len(motifs) # count matrix converted to profile

    return count_matrix

def get_probability_vector(pattern, profile):
    """ Returns the probability vector of the kmers in the deleted sequence given profile matrix """
    p_vector = list()
    k = len(profile[0])
    letter_to_rowind_map = {'A': 0, 'T':1, 'G':2, 'C':3}

    for start in range(len(pattern) - k):
        kmer = pattern[start: start + k] 
        proba = 1
        # calculate probability of occurence of the give k-mer
        for ind in range(k):
            proba = proba*profile[letter_to_rowind_map[kmer[ind]]][ind]
        p_vector.append(proba)
    return p_vector

def roll_die(p_vector, seed):
    """ Picks a kmer index by **sampling** from a weighted probability vector """
    # random.seed(seed)    # Makes the results reproducible as the seed is tied to the LNG
    # values, weights = list(range(len(p_vector))), p_vector
    # total = 0
    # cum_weights = []
    # for w in weights:
    #     total += w
    #     cum_weights.append(total)
    # x = random.random() * total
    # i = bisect(cum_weights, x)
    # return i

    return random.choices(range(len(p_vector)), p_vector)[0]

def gibbs_sampler(sequences, k, t, n, random_start_seed):
    motif_start_indices, seed = lin_cong_generator(X = random_start_seed, num = t)
    # To ensure motif start indice is in correct range i.e. (0, len(sequences-k) take the modulus. Below code only works if are t sequences are of same length.

    temp = [len(seq)-k for seq in sequences]
    motif_start_indices = [motif_start_indices[ind]%temp[ind] for ind in range(t)]

    # Get all motifs
    motifs = list()
    for ind in range(t):
        motifs.append(sequences[ind][motif_start_indices[ind]: motif_start_indices[ind] + k])

    best_motifs_indices = motif_start_indices
    best_motifs_score = get_score(motifs)

    for j in range(n):
        rnd_num, _ = lin_cong_generator(X = seed, num = 1)
        # To ensure i lies in range (0, t) take modulus
        i = rnd_num[0] % t

        # Generate profile of all motifs except motif i.
        profile = get_profile(motifs[:i]+motifs[i+1:])
        
        # Probabilistically obtain a motif in the ignored sequence i (by sampling)
        pattern = sequences[i]
        p_vector = get_probability_vector(pattern, profile)

        # Normalize the p_vector so that probabilities add to 1
        p_vector = [p/sum(p_vector) for p in p_vector]
        new_motif_index = roll_die(p_vector, rnd_num[0])

        # Update motif i in motifs list
        motifs[i] = sequences[i][new_motif_index: new_motif_index + k]
        
        # Compare the score between the current motif and best motif
        motif_score = get_score(motifs)
        if motif_score < best_motifs_score:
            new_motif_indices = []
            for seq, motif in zip(sequences, motifs):
                new_motif_indices.append(seq.index(motif))
            best_motifs_indices = new_motif_indices
            best_motifs_score = motif_score

        seed = rnd_num[0]
    
    return best_motifs_indices, best_motifs_score

def get_consensus(sequences, motif_start_indices, k):
    motifs = []
    t = len(motif_start_indices)

    for ind in range(t):
        # print(sequences[ind][motif_start_indices[ind]: motif_start_indices[ind] + k])
        motifs.append(sequences[ind][motif_start_indices[ind]: motif_start_indices[ind] + k])
    count_matrix = get_count_matrix(motifs)

    l = list(zip(*count_matrix))

    consensus = ''
    num_to_letter_map = {0:'A', 1:'T', 2:'G', 3:'C'}

    for i in range(k):
        
        f = lambda x: l[i][x]
        consensus += num_to_letter_map[max(range(4), key=f)] 
    return consensus    

def main():
    FNAME = ''
    ktn= input()
    k, t, n = ktn.split(' ')
    k, t, n = int(k), int(t), int(n)

    # Enter file location in FNAME
    # FNAME = 'P53.fa'

    # If sequences are given directly to terminal
    if len(FNAME) == 0:
        sequences = list()
        for sequence in sys.stdin:
            if sequence == '\n': break
            sequence = sequence.strip()
            sequences.append(sequence)

    # If fasta file is given for the sequences
    else:

        sequences = read_fasta(FNAME)

    # print('Running')
    NUM_RANDOM_STARTS = 2

    # These two variables store the best score and motifs over multiple number of random runs
    overall_best_score, overall_best_motif_indices = k*t, []

    for random_start in range(NUM_RANDOM_STARTS):
        print('Random start ', random_start)
        # random_start_seed or its function will determine Xo in our LCG.
        best_motifs_indices, best_motifs_score = gibbs_sampler(sequences, k, t, n, random_start+82)

        if best_motifs_score < overall_best_score:
            overall_best_motif_indices = best_motifs_indices
            overall_best_score = best_motifs_score


            # os.system('cls')
            # print('-'*50)
            print('run: ',random_start,'; Score: ', overall_best_score)
            # print('consensus string:', get_consensus(sequences, overall_best_motif_indices, k))
            # print('-'*50)

    print('-'*50)
    print('(k, t, n, NUM_STARTS, Score) : ({}, {}, {}, {}, {})'.format(k, t, n, NUM_RANDOM_STARTS, overall_best_score))
    for ind in range(t):
        print(sequences[ind][overall_best_motif_indices[ind]: overall_best_motif_indices[ind] + k])

if __name__ == '__main__':
    t1 = time()
    main()
    print('Time taken (in mins): ', round((time()-t1)/60, 2))

    ## Get motif score for a collection of input sequences
    # sequences = list()
    # for sequence in sys.stdin:
    #     if sequence == '\n': break
    #     sequence = sequence.strip()
    #     sequences.append(sequence)
    # # print(get_score(sequences))
    # print(get_consensus(sequences, [0]*721, 11))
    # print(get_profile(sequences)[3])
