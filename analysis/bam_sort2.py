#!/usr/bin/env python3

"""
Input is VCF.tsv, .bam file
Output is a tsv that has all the SNP positions with 
"""

import sys
import pickle
import time

## variables
scaffolds = [
    # 'X',
    # 'Y',
    '2L',
    # '2R',
    # '3L',
    # '3R',
    # '4'
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

## make the counting dictionary
countDict = {}
for scaf in scaffolds:
    countDict[scaf] = {}
    for kpos in snpDict[scaf].keys():
        countDict[scaf][kpos] = {'mel' : 0, 
                                 'sim' : 0, 
                                 'sec' : 0,
                                 'melORsim' : 0,
                                 'melORsec' : 0, 
                                 'simORsec' : 0
        }

print("Started analysis for " + str(bamFile))

## load in the bam file, and process each read
with open (bamFile, 'r') as bam:
    i = 0
    for line in bam:
        if ('@SQ' not in line) and ('@PG' not in line) and ('@HD' not in line) and ('@RG' not in line):
            # line.split('\t')[2] is the CHR, [3] is the starting position, [9] is the sequence
            i += 1
            chrom = line.split('\t')[2]
            posit = line.split('\t')[3]
            seq = line.split('\t')[9]
            ## quick check
            if chrom not in scaffolds:
                continue
            
            ## sort through the postions in a read and find out if it matches anything
            for pos in range(0, len(seq)):
                genPOS = pos + int(posit)
                if genPOS in snpDict[chrom]:
                    print(snpDict[chrom][genPOS], seq[pos-1:pos+2], seq[pos])
                     
                    if (snpDict[chrom][genPOS][0] == seq[pos]) and (snpDict[chrom][genPOS][1] != seq[pos]) and (snpDict[chrom][genPOS][2] != seq[pos]):
                        countDict[chrom][genPOS]['mel'] += 1
                    if (snpDict[chrom][genPOS][0] != seq[pos]) and (snpDict[chrom][genPOS][1] == seq[pos]) and (snpDict[chrom][genPOS][2] != seq[pos]):
                        countDict[chrom][genPOS]['sim'] += 1
                    if (snpDict[chrom][genPOS][0] != seq[pos]) and (snpDict[chrom][genPOS][1] != seq[pos]) and (snpDict[chrom][genPOS][2] == seq[pos]):
                        countDict[chrom][genPOS]['sec'] += 1
                    if (snpDict[chrom][genPOS][0] == seq[pos]) and (snpDict[chrom][genPOS][1] == seq[pos]) and (snpDict[chrom][genPOS][2] != seq[pos]):
                        countDict[chrom][genPOS]['melORsim'] += 1
                    if (snpDict[chrom][genPOS][0] == seq[pos]) and (snpDict[chrom][genPOS][1] != seq[pos]) and (snpDict[chrom][genPOS][2] == seq[pos]):
                        countDict[chrom][genPOS]['melORsec'] += 1
                    if (snpDict[chrom][genPOS][0] != seq[pos]) and (snpDict[chrom][genPOS][1] == seq[pos]) and (snpDict[chrom][genPOS][2] == seq[pos]):
                        countDict[chrom][genPOS]['simORsec'] += 1
            
            ## generate some terminal output
            if not (i % 1000000):
                print(str(i/1000000) + ' MIL reads analyzed for ' + str(bamFile) + ' : ' + str(time.clock()-timeStart) + ' seconds')

with open(outFile, 'w') as printOut:
    printOut.write('CHRM\tsnpPOS\tmelREADs\tsimREADs\tsecREADs\tmelORsim\tmelORsec\tsimORsec')
    for scaf in countDict:
        for kpos in countDict[scaf]:
            printOut.write('\n' + str(scaf) + '\t' +
                str(kpos) + '\t' + 
                str(countDict[scaf][kpos]['mel']) + '\t' +
                str(countDict[scaf][kpos]['sim']) + '\t' +
                str(countDict[scaf][kpos]['sec']) + '\t' +
                str(countDict[scaf][kpos]['melORsim']) + '\t' +
                str(countDict[scaf][kpos]['melORsec']) + '\t' +
                str(countDict[scaf][kpos]['simORsec']))

# measure time
print(str(time.clock()-timeStart) + ' seconds to sort through bam file')