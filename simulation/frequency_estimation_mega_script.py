#! /usr/bin/env python3

import random
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks_cwt


LENGTH = 50000000
N_TRIALS = 350
selection_strength = 1

replicates = 3000
plot_all = False
window = 1000000
step = 10000

out_pre = '/Volumes/Jacob_2TB_storage/sim_sec_recombination_mapping/genomics_scripts/analysis/'


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
    """Compute the gamma distrubtion model for recombination.
    Return a dictionary with 100 steps (should never need that many)
    """
    model = {}
    prev = 0
    for i in range(1,101):
        # the computation:
            # (5*1morgan) (sum_to_i)(
            #   1 / (5i-1) * 1M^(5i) * dist^(5i-1) * e^(-(1M)*dist)  
            # )
            # where distance is in morgans
        prev += ((1**(5*i)) * ((0.01*i)**((5*i)-1)) * (np.e ** (-(0.01*i)))) / (5*(i) - 1)
        model[i] = 5 * prev
    return model

def recombine(seq1, seq2, model):
    """Make a recombinant sequence from sequence 1 and 2"""

    idxs = [idx for idx, x in enumerate(seq1)]

    co1 = random.choice(idxs)

    half = int(len(idxs) / 2)

    if co1 < half:
        co2 = random.choice(idxs[half:])
    else:
        co2 = random.choice(idxs[:half])

    co1, co2 = reorder(co1, co2)

    if bool_chance(0.5):
        rec_seq = seq1[:co1] + seq2[co1:co2] + seq1[co2:]
    else:
        rec_seq = seq2[:co1] + seq1[co1:co2] + seq2[co2:]

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

def sex_difference(male_reps, female_reps):
    reps = {}

    for i in male_reps:
        reps[i] = []
        for idx, entry in enumerate(male_reps[i]):
            male = male_reps[i][idx]
            female = female_reps[i][idx]
            freq = [(male[0]), (male[1] - female[1]), (male[2] - female[2]), (male[3] - female[3])]
            reps[i].append(freq)
    
    return reps

def average_replicates(reps):
    table = []
    for pos, lis in enumerate(reps[0]):
        melav = (reps[0][pos][1] + reps[1][pos][1] + reps[2][pos][1]) / 3
        simav = (reps[0][pos][2] + reps[1][pos][2] + reps[2][pos][2]) / 3
        secav = (reps[0][pos][3] + reps[1][pos][3] + reps[2][pos][3]) / 3
        table.append([reps[0][pos][0], melav, simav, secav])
    return table

def rolling_average(table, window, step):
    
    win2 = window/2
    pos = window/2

    window_table = []
    posits = [x[0] for x in table]
    while pos < max(posits):
        melav = []
        simav = []
        secav = []
        for x in table:
            if (x[0] > pos-win2):
                if (x[0] < pos+win2):
                    melav.append(x[1])
                    simav.append(x[2])
                    secav.append(x[3])
                else:
                    break
        if melav and simav and secav:
            window_table.append([pos, np.mean(melav), np.mean(simav), np.mean(secav)])
        pos += step

    
    newtable = [['position','mel','sim','sec']]
    for x in window_table:
        newtable.append(x)

    return newtable

def estimate_max(table):
    table = table[1:]
    sim_freqs = [(x[0],(x[2]-x[3]))for x in table]
    sim_freqs = sorted(sim_freqs, key=lambda x: x[1])

    # find the peaks
    xs = [x[1] for x in sim_freqs]
    peaks = list(find_peaks_cwt(xs, np.arange(50, 200)))

    # this produces a list. Find the biggest one in the list
    big = (0,0)
    for peak in peaks:
        if sim_freqs[peak][1] > big[1]:
            big = (sim_freqs[peak][0], sim_freqs[peak][1])

    return big[0]

def write_to_csv(table, filename):
    with open(filename, 'w') as f:
        for line in table:
            line = [str(x) for x in line]
            f.write(','.join(line))
            f.write('\n')

def plot_frequencies(filename, sel_spot, esitmate, num):
    df = pd.read_csv(filename)
    fig = plt.figure(figsize=(6, 8))
    plt.plot(df['position'], df['mel'], color = 'blue', label = 'D.mel')
    plt.plot(df['position'], df['sim'], color = 'orange', label = 'D.sim')
    plt.plot(df['position'], df['sec'], color = 'red', label = 'D.sec')
    plt.axvline(x=sel_spot, color='black', label = 'actual site')
    plt.axvline(x=esitmate, color='green', label = 'estimated site')
    plt.ylim(-0.4,0.4)
    plt.legend()
    plt.ylabel('Allele Frequency (Male - Female)')
    plt.xlabel('Genomic position')
    plotname = filename.split('.csv')[0] + str(num) + '.pdf'
    plt.savefig(plotname)


parents, id_spots = generate_parents()
model = gamma_model()
differences = []
for i2 in range(replicates):
    sel_idx, sel_spot = pick_selected_index(id_spots)
    male_reps = {}
    female_reps = {}
    reps = {}
    # replicate 3 times
    for i in range(3):
        male_pool = make_pool(parents, sel_idx, 'male', selection_strength, model)
        female_pool = make_pool(parents, sel_idx, 'female', selection_strength, model)
        male_reps[i] = simulate_sequencing(male_pool, id_spots)
        female_reps[i] = simulate_sequencing(female_pool, id_spots)
    reps = sex_difference(male_reps, female_reps)
    table = average_replicates(reps)
    table = rolling_average(table, window, step)
    estimated_site = estimate_max(table)
    differences.append(estimated_site-sel_spot)
    if plot_all:
        outtable = out_pre + str(i2) + 'simulation_out.csv'
        write_to_csv(table, outtable)
        plot_frequencies(outtable, sel_spot, estimated_site, i2)


print(differences)
print("2x std:", np.std(differences))
print("Average:", np.average(differences))
plt.hist(differences, bins=20)
plt.savefig(out_pre+'confidence_hist.pdf')

outtable = out_pre + 'final_simulation_out.csv'
write_to_csv(table, outtable)
plot_frequencies(outtable, sel_spot, estimated_site, 'final')
