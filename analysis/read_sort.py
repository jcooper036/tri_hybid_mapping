#!/usr/bin/env python3

"""
Input is VCF.tsv, .bam file
Output is a pickle file that is the dictionary of all genome postions
"""

import sys
import pickle
import time
import re


def cigarParse(genPos, seq, cigarString):
    """returns a new sequence and/or modified start position based on the cigar string"""
    
    cigarString = re.split('(\d+)',cigarString)
    oper = []
    cigarString = cigarString[1:]
    for idx, item in enumerate(cigarString):
        if item.isdigit():
            ent = (int(item),cigarString[idx+1])
            oper.append(ent)
        
        if (len(oper) == 0):
            return genPos, seq, False

    readPos = 0
    for idx, ent in enumerate(oper):
        ## soft clip
        if idx == 0:
            # if ent[1] == 'S':
            #     seq = seq[ent[0]:] # clip out the sequence that soft clip actually clipped
            if ent[1] == 'S' or ent[1] == 'I' or ent[1] == 'D':
                return genPos, seq, False
        
        else:
            if ent[1] == 'S':
                seq = seq[:-ent[0]] # clip out the sequence that soft clip actually clipped
            # if ent[1] == 'S':
            #     return genPos, seq, False


            if ent[1] == 'I':
                seq = seq[:readPos] + seq[(readPos + ent[0]):]

            if ent[1] == 'P':
                seq = seq[:readPos] + seq[(readPos + ent[0]):]

            if ent[1] == 'D':
                add = '-' * ent[0]
                seq = seq[:readPos] + add + seq[readPos:]            
        
        ## add on the number of BP that have been covered
        readPos += ent[0]


    return genPos, seq, True

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
printOut.write('CHRM\tPOS1\tPOS2\tTYPE1\tTYPE2\tRead1\tRead2\tSAMPLE\n')
print('CHRM\tPOS1\tPOS2\tTYPE1\tTYPE2\tRead1\tRead2\tSAMPLE')
## load in the bam file, and process each read
with open (bamFile, 'r') as bam:
    readDict = {}
    saveDict = {}
    for line in bam:
        if ('@SQ' not in line) and ('@PG' not in line) and ('@HD' not in line) and ('@RG' not in line):
            # line.split('\t')[2] is the CHR, [3] is the starting position, [9] is the sequence
            read = line.split('\t')[0]
            chrom = line.split('\t')[2]
            posit = line.split('\t')[3]
            seq = line.split('\t')[9]
            qualCheck = False
            if int(line.split('\t')[4]) > 36: qualCheck = True
            sample = line.split('RG:Z:')[1].split('\t')[0]
            mateCheck = False
            if line.split('\t')[6] == '=': 
                mateCheck = True
                matePos = line.split('\t')[7]
            cigarString = line.split('\t')[5]
            mel = 0
            sim = 0
            sec = 0
            ukn = True
            readType = 'ukn'
            if sample not in readDict:
                readDict[sample] = {}
            if sample not in saveDict:
                saveDict[sample] = {}                
            
            ## modify the sequence and start position based of the cigar string
            posit, seq, cigFlag = cigarParse(posit, seq, cigarString)

            ## just don't do anything if these fail
            if chrom not in scaffolds or not qualCheck or not mateCheck or not cigFlag:
                continue

            ## sort through the postions in a read and find out if it matches anything
            for pos in range(0, len(seq)):
                genPOS = pos + int(posit)
                if genPOS in snpDict[chrom]:

                    # print(snpDict[chrom][genPOS], seq[pos-1:pos+2], seq[pos], sample)
                    ukn = False
                    if snpDict[chrom][genPOS][0] == seq[pos]:
                        mel += 1
                    if snpDict[chrom][genPOS][1] == seq[pos]:
                        sim += 1
                    if snpDict[chrom][genPOS][2] == seq[pos]:
                        sec += 1

            ## now classify the read

            # if (mel > 0) and (sim == 0) and (sec == 0):
            #     readType = 'mel'
            if (mel == 0) and (sim > 0) and (sec == 0):
                readType = 'sim'            
            elif (mel == 0) and (sim == 0) and (sec > 0):
                readType = 'sec'             
            else:
                ukn = True

            # if (mel < sim) and (mel < sec):
            #     if sim > sec:
            #         readType = 'sim'
            #     elif sec > sim:
            #         readType = 'sec'
            #     else:
            #         ukn = True 
            # else:
            #     ukn = True

            # add it to the dictionary of the read file
            if not ukn and mateCheck:
                if (read in readDict[sample]) and (matePos == readDict[sample][read]['pos1']):
                    if readDict[sample][read]['type1'] != readType:
                        saveDict[sample][read] = readDict[sample][read]
                        saveDict[sample][read]['pos2'] = posit
                        saveDict[sample][read]['type2'] = readType
                        saveDict[sample][read]['read2'] = read
                        printOut.write(str(saveDict[sample][read]['CHROM']) + '\t' +
                                       str(saveDict[sample][read]['pos1']) + '\t' +
                                       str(saveDict[sample][read]['pos2']) + '\t' +
                                       str(saveDict[sample][read]['type1']) + '\t' +
                                       str(saveDict[sample][read]['type2']) + '\t' +
                                       str(saveDict[sample][read]['read1']) + '\t' +
                                       str(saveDict[sample][read]['read2']) + '\t' +
                                       str(sample) + '\n')
                        print(str(saveDict[sample][read]['CHROM']) + '\t' +
                              str(saveDict[sample][read]['pos1']) + '\t' +
                              str(saveDict[sample][read]['pos2']) + '\t' +
                              str(saveDict[sample][read]['type1']) + '\t' +
                              str(saveDict[sample][read]['type2']) + '\t' +
                              str(saveDict[sample][read]['read1']) + '\t' +
                              str(saveDict[sample][read]['read2']) + '\t' +
                              str(sample))
                else:
                    readDict[sample][read] = {'pos1' : posit, 
                                              'pos2' : '', 
                                              'type1' : readType, 
                                              'type2' : 'ukn',
                                              'CHROM' : chrom,
                                              'read1' : read}

## close the file handle
printOut.close()

# measure time
print(str(time.clock()-timeStart) + ' seconds to process bam file')