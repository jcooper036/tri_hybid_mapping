#! /usr/bin/env python3

import random
import copy
import numpy as np
from scipy.signal import find_peaks_cwt
import pickle


LENGTH = 50000000
N_TRIALS = 350
selection_strength = 1

replicates = 10
window = 1000000
step = 20000

out_pre = '/Volumes/Jacob_2TB_storage/sim_sec_recombination_mapping/simulation/'

def load_tsv(idx, out_pre):
    table = {}

    # load files
    male_reps = {}
    female_reps = {}
    for idx2 in range(3):
        male_reps[idx2] = []
        filename = out_pre + 'data/simulated_data/' + str(idx) + '_' + str(idx2) + '_male.tsv'
        with open(filename, 'r') as f:
            for line in f:
                line = line.rstrip()
                if '@@@' in line:
                    sel_spot = int(line.split(':')[1])
                elif 'CHROM' not in line:
                    line = line.split('\t')
                    male_reps[idx2].append(line)
        female_reps[idx2] = []
        filename = out_pre + 'data/simulated_data/' + str(idx) + '_' + str(idx2) + '_female.tsv'
        with open(filename, 'r') as f:
            for line in f:
                line = line.rstrip()
                if '@@@' in line:
                    sel_spot = int(line.split(':')[1])
                elif 'CHROM' not in line:
                    line = line.split('\t')
                    female_reps[idx2].append(line)

    return male_reps, female_reps, sel_spot

def window_average(reps, window, step):

    win_reps = {}

    for rep in reps:
        win2 = window/2
        pos = window/2
        winds = []
        posits = [int(x[1]) for x in reps[rep]]
        while pos < max(posits):
            melav = []
            simav = []
            secav = []
            start = 0
            for idx in range(start, len(reps[rep])):
                x = reps[rep][idx]
                if (int(x[1]) > pos-win2):
                    start = idx
                    if (int(x[1]) < pos+win2):
                        melav.append(float(x[2]))
                        simav.append(float(x[3]))
                        secav.append(float(x[4]))
                    else:
                        break                

            if melav and simav and secav:
                winds.append([pos, np.mean(melav), np.mean(simav), np.mean(secav)])
            pos += step
        win_reps[rep] = winds

    return win_reps

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

differences = []
for i2 in range(replicates):
    male_reps, female_reps, sel_spot = load_tsv(i2, out_pre)
    male_reps = window_average(male_reps, window, step)
    female_reps = window_average(female_reps, window, step)
    table = sex_difference(male_reps, female_reps)
    table = average_replicates(table)
    estimated_site = estimate_max(table)
    out = {
        'est_site' : estimated_site,
        'difference' : estimated_site-sel_spot,
        'sel_site' : sel_spot,
        'table' : table
    }
    pickle_file = out_pre + 'data/parsed_data/' + str(i2) + '.pkl'
    with open(pickle_file, 'wb') as f:
        pickle.dump(out, f, pickle.HIGHEST_PROTOCOL)


# print(differences)
# print("2x std:", np.std(differences))
# print("Average:", np.average(differences))
# plt.hist(differences, bins=20)
# plt.savefig(out_pre+'confidence_hist.pdf')

# plot_frequencies(table, sel_spot, estimated_site, 'final')
