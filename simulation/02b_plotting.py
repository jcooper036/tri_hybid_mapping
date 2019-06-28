#! /usr/bin/env python3

import random
import copy
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks_cwt


LENGTH = 50000000
N_TRIALS = 350
selection_strength = 1

replicates = 10
window = 1000000
step = 20000

out_pre = '/Volumes/Jacob_2TB_storage/sim_sec_recombination_mapping/simulation/'

def load_table(filename):
    with open(filename, 'rb') as f:
        data = pickle.load(f)
        return data['table'], data['sel_site'], data['est_site']

def plot_frequencies(table, sel_spot, esitmate, num, out_pre):
    fig = plt.figure(figsize=(6, 8))
    table = np.array(table)
    plt.plot(table[::,0], table[::,1], color = 'blue', label = 'D.mel')
    plt.plot(table[::,0], table[::,2], color = 'orange', label = 'D.sim')
    plt.plot(table[::,0], table[::,3], color = 'red', label = 'D.sec')
    plt.axvline(x=sel_spot, color='black', label = 'actual site')
    plt.axvline(x=esitmate, color='green', label = 'estimated site')
    plt.ylim(-0.4,0.4)
    plt.legend()
    plt.ylabel('Allele Frequency (Male - Female)')
    plt.xlabel('Genomic position')
    plotname = out_pre + 'simulated_graph.pdf'
    plt.savefig(plotname)


differences = []
for i2 in range(replicates):
    filename = out_pre + 'data/parsed_data/' + str(i2) + '.pkl'
    table, sel_spot, estimated_site = load_table(filename)
    differences.append(estimated_site-sel_spot)

print(differences)
print("2x std:", np.std(differences))
print("Average:", np.average(differences))
plt.hist(differences, bins=20)
plt.savefig(out_pre+'confidence_hist.pdf')

plot_frequencies(table, sel_spot, estimated_site, 'final', out_pre)
