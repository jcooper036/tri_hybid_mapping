#!/usr/bin/env python3

"""
Input is VCF.tsv, .bam file
Output is a pickle file that is the dictionary of all genome postions
"""

import sys
import pickle
import time

## variables
scaffolds = [
    'X',
    'Y',
    '2L',
    '2R',
    '3L',
    '3R',
    '4'
]

bamFile = sys.argv[1]
vcfTable = sys.argv[2]
outFile = sys.argv[3]

timeStart = time.clock()

## make a dictionary of all postions by scaffold
countTable = {}
for scaf in scaffolds:
    countTable[scaf] = {}

## load in the SNP table
snpDict = {}
for scaf in scaffolds:
    snpDict[scaf] = {}
with open (vcfTable) as vcf:
    for line in vcf:
        line = line.strip()
        if 'CHRM' not in line:
            line = line.split('\t')
            if line[0] not in snpDict.keys():
                snpDict[line[0]] = {}
            snpDict[line[0]][int(line[1])] = [line[2], line[3], line[4]]

## make the header for the out file
printOut = open(outFile, 'w')
printOut.write('CHRM\tStartPOS\tEndPOS\ttype')

## load in the bam file, and process each read
with open (bamFile, 'r') as bam:
    for line in bam:
        if ('@SQ' not in line) and ('@PG' not in line) and ('@HD' not in line) and ('@RG' not in line):
            # line.split('\t')[2] is the CHR, [3] is the starting position, [9] is the sequence
            print(line)
            chrom = line.split('\t')[2]
            posit = line.split('\t')[3]
            seq = line.split('\t')[9]
            mel = 0
            sim = 0
            sec = 0
            ukn = True
            if chrom not in scaffolds:
                continue
            
            ## sort through the postions in a read and find out if it matches anything
            for pos in range(0, len(seq)):
                genPOS = pos + int(posit)
                if genPOS in snpDict[chrom].keys():

                    # print(snpDict[chrom][genPOS])
                    # print(line)
                    ukn = False
                    if snpDict[chrom][genPOS][0] == seq[pos]:
                        mel += 1
                    if snpDict[chrom][genPOS][1] == seq[pos]:
                        sim += 1
                    if snpDict[chrom][genPOS][2] == seq[pos]:
                        sec += 1
            
            ## now classify the read
            if ukn:
                readType = 'ukn'
            if (mel > sim) and (mel > sec):
                readType = 'mel'
            if (mel < sim) and (mel < sec):
                readType = 'simORsec'
            if readType == 'simORsec':
                if sim > sec:
                    readType = 'sim'
                if sec > sim:
                    readType = 'sec'
            # if (mel == 0) and (sim == 0) and (sec > 0):
            #     readType = 'sec'
            # if (mel > 0) and (sim == 0) and (sec == 0):
            #     readType = 'mel'
            # if (mel == 0) and (sim > 0) and (sec == 0):
            #     readType = 'sim'
            # if (mel == 0) and (sim > 0) and (sec > 0) and (sim == sec):
            #     readType = 'simORsec'
            # if (mel == 0) and (sim > 0) and (sec > 0) and (sim > sec):
            #     readType = 'sim'
            # if (mel == 0) and (sim > 0) and (sec > 0) and (sim < sec):
            #     readType = 'sec'

            # write it to the output file
            if readType != 'ukn':
                printOut.write('\n' + str(chrom) + '\t' +
                            str(posit) + '\t' + 
                            str(int(posit) + len(seq)-1) + '\t' + 
                            str(readType))


            # now add the read's ID to the table
            # for pos in range(0, len(seq)):
            #     genPOS = pos + int(posit)
            #     if genPOS in countTable[chrom].keys():
            #         countTable[chrom][genPOS][readType] += 1
            #     else:
            #         countTable[chrom][genPOS] = {readType : 1}

## close the file handle
printOut.close()

# measure time
print(str(time.clock()-timeStart) + ' seconds to pickle bam file')