#!/usr/bin/env python3
"""
Looks for common windows between two window average tsv files, and normalizes the second file by the first file

Input: Two window average tsvs, the first is the control and the second the experiment. Name of tsv output

"""

import time
import sys


controlFile = sys.argv[1]
experFile = sys.argv[2]
outName = sys.argv[3]

def readTsv(file):
    """Reads a TSV file"""
    file = open(file, 'r')
    outdict = {}
    ## skip the first line
    i = 0
    for line in file.readlines():
        if i > 0:
            line = line.strip()
            line = line.split('\t')
            if line[0] not in outdict:
                outdict[line[0]] = {}
            outdict[line[0]][line[1]] = [line[2], line[3], line[4]]
        i += 1     
    file.close()
    return outdict

def normalizeTsv(controlFile, experFile):
    """Gets dictionaries, finds difference between values, returns dictionary"""
    controlDict = readTsv(controlFile)
    experDict = readTsv(experFile)
    diffDict = {}
    for scaff in controlDict:
        if scaff in experDict:
            if scaff not in diffDict:
                diffDict[scaff] = {}
            for pos in controlDict[scaff]:
                if pos in experDict[scaff]:
                    diffDict[scaff][pos] = [float(experDict[scaff][pos][idx]) - float(controlDict[scaff][pos][idx]) for idx in range(0, len(experDict[scaff][pos]))]
                    
    return diffDict

def outWrite(outFile, diffDict):
    """Writes the results to an output file"""
    with open(outFile, 'w') as f:
        f.write('CHROM\tPOS\tmel\tsim\tsec\n')
        for scaff in diffDict:
            for pos in diffDict[scaff]:
                f.write(str(scaff) + '\t' +
                        str(pos) + '\t' +
                        str(diffDict[scaff][pos][0]) + '\t' +
                        str(diffDict[scaff][pos][1]) + '\t' +
                        str(diffDict[scaff][pos][2]) + '\n')

diffDict = normalizeTsv(controlFile, experFile)
outWrite(outName, diffDict)