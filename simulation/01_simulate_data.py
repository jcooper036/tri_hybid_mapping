#! /usr/bin/env python3

import random
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks_cwt


LENGTH = 50000000
N_TRIALS = 350
selection_strength = 0.746

replicates = 10


out_pre = '/Volumes/Jacob_2TB_storage/sim_sec_recombination_mapping/simulation/'


class Sample_pool():
    """Holds samples, which need to be Sequence objects"""

    def __init__(self):
        self.pool = {}
        self.keys = []
    
    def add_sample(self, sequence):
        if self.keys:
            newkey = max(self.keys) + 1
        else:
            newkey = 0
        self.pool[newkey] = sequence
        self.keys = list(self.pool.keys())

    def pick_sample(self):
        """Selects a random sample"""
        sample = random.choice(self.keys)
        return self.pool[sample]

    def sample_freq(self, position):
        """Gathers a sample allele frequency for a given position"""
        reps = random.randint(36,40)
        mel = 0
        sim = 0
        sec = 0
        for i in range(reps):
            sample = self.pick_sample()
            if sample[position] == 1:
                mel += 1
            elif sample[position] == 2:
                sim += 1
            elif sample[position] == 3:
                sec += 1

        melf = mel / (mel + sim + sec)
        simf = sim / (mel + sim + sec)
        secf = sec / (mel + sim + sec)

        return melf, simf, secf

def bool_chance(num):
    x = random.random()
    if x < num:
        return True
    else:
        return False

def generate_parents():
    """
    Generate a set of 3 parents where all positions are 0, if mel specific 1,
    if sim specific 2, and if sec specific 3.
    Return Dicitonary of parents
    """

    everyN = 6000
    
    # get a list of the positions that would be differentiating spots
    id_spots = [random.randint(0,LENGTH) for i in range(int(LENGTH/everyN))]
    id_spots.sort()

    # make a list of the same length for each species
    mel_seq = [1] * len(id_spots)
    sim_seq = [2] * len(id_spots)
    sec_seq = [3] * len(id_spots)

    parents = {
        'mel' : mel_seq,
        'sim' : sim_seq,
        'sec' : sec_seq,

    }

    return parents, id_spots

def pick_selected_index(id_spots):
    """
    Takes the list of different sites between the three sequences
    Picks a random spot in the middle half of the sequences to be the selected site
    Returns the index of the differentiated spots list that is closest to the 
        selected point, and the selected point
    """

    q1 = int(LENGTH / 8)
    q3 = int(7 * LENGTH / 8)
    sel_spot = random.randint(q1, q3)

    hold = 0
    for idx, x in enumerate(id_spots):
        if x < sel_spot:
            hold = x
        else:
            if (sel_spot - hold) <= (x - sel_spot):
                return (idx-1), sel_spot
            else:
                return idx, sel_spot

def reorder(one, two):
    if one >= two:
        return two, one
    else:
        return one, two

def gamma_model():
    model = {}
    for z in range(1,98):
        z /= 100
        summation = 0
        for n in range(1,10000):
            term1 = 1 / (5*(n) - 1)
            term2 = 1 ** (5*n)
            term3 = z ** (5*(n) - 1)
            term4 = np.e ** (-(z))
            summation += term1 * term2 * term3 * term4        
        
        model[z] = 5 * summation
    model[98] = 1
    return model

def co_math(seq1, seq2, cos):
    if len(cos) > 1:
        req_seq = seq1[:cos[0]] + seq2[cos[0]:cos[1]] + seq1[cos[1]:]
        return co_math(req_seq, seq2, cos[2:])
    else:
        if len(cos) == 1:
            return seq1[:cos[0]] + seq2[cos[0]:]
        else:
            return seq1
        

def recombine(seq1, seq2, model):
    """Make a recombinant sequence from sequence 1 and 2"""

    idxs = [idx for idx, x in enumerate(seq1)]

    cents = int(len(idxs) / 100)

    cos = {}
    cos[0] = random.choice(idxs)

    # check to the left
    last_co = cos[0]
    search = True
    while search:
        for i in range(1,120):
            if (last_co - i*cents) < 0:
                search = False
                break
            if bool_chance(model[i/100]):
                new_entry = max(list(cos.keys()))+1
                center = last_co - (i*cents)
                cos[new_entry] = random.randint(max([int(center - cents/2),0]), min([int(center + cents/2),len(idxs)]))
                last_co = cos[new_entry]
                break

    # check to the right
    last_co = cos[0]
    search = True
    while search:
        for i in range(1,120):
            if (last_co + i*cents) > len(idxs):
                search = False
                break
            if bool_chance(model[i/100]):
                new_entry = max(list(cos.keys()))+1
                center = last_co + (i*cents)
                cos[new_entry] = random.randint(max([int(center - cents/2),0]), min([int(center + cents/2),len(idxs)]))
                last_co = cos[new_entry]
                break

    crossovers = [cos[x] for x in cos]
    crossovers.sort()

    if bool_chance(0.5):
        rec_seq = co_math(seq1, seq2, crossovers)
    else:
        rec_seq = co_math(seq2, seq1, crossovers)

    return rec_seq

def filter(seq, sel_idx, sel_strength):
    """Checks if a sequence contains the selected site"""
    if seq[sel_idx] == 2:
        return True
    else:
        if bool_chance(1-sel_strength):
            return True
        else:
            return False

def make_pool(parents, sel_idx, sex, sel_strength, model):
    """
    Input: mel, sim, sec parent sequence  , site to be selected on
    Returns: sample pool of recombinant sequences and mel sequences
    """

    pool = Sample_pool()
    for i in range(N_TRIALS):
        survive = False
        while not survive:
            recseq = recombine(parents['sim'], parents['sec'], model)
            survive = filter(recseq, sel_idx, sel_strength)
            if sex == 'female':
                survive = True
        pool.add_sample(recseq)
        # because diploids, we add one mel sequence for each recombinant sequence
        pool.add_sample(parents['mel'])

    return pool

def simulate_sequencing(pool, id_spots):
    table = []
    for idx, x in enumerate(id_spots):
        mel, sim, sec = pool.sample_freq(idx)
        table.append([x,mel,sim,sec])
    return table


def write_to_tsv(samples, replciate, outfile_prefix, sex, selected_spot):
    for subrep in samples:
        table = samples[subrep]
        filename = outfile_prefix + '/data/simulated_data/' + str(replciate) + '_' + str(subrep) + '_' + sex + '.tsv'
        with open(filename, 'w') as f:
            f.write('@@@@header:' + str(sel_spot) + '\n')
            f.write('CHROM\tPOS\tmel\tsim\tsec\n')
            for line in table:
                line = [str(x) for x in line]
                f.write('simulated\t')
                f.write('\t'.join(line))
                f.write('\n')


parents, id_spots = generate_parents()
model = gamma_model()
differences = []
for i2 in range(replicates):
    sel_idx, sel_spot = pick_selected_index(id_spots)
    male_reps = {}
    female_reps = {}
    # replicate 3 times
    for i in range(3):
        male_pool = make_pool(parents, sel_idx, 'male', selection_strength, model)
        female_pool = make_pool(parents, sel_idx, 'female', selection_strength, model)
        male_reps[i] = simulate_sequencing(male_pool, id_spots)
        female_reps[i] = simulate_sequencing(female_pool, id_spots)

    write_to_tsv(male_reps, i2, out_pre, 'male', sel_spot)
    write_to_tsv(female_reps, i2, out_pre, 'female', sel_spot)
