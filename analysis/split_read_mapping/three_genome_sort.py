#!/usr/bin/env python3

"""
Input is 3 files (in this order)
    - sam file of all the reads aligned to all 3 genomes
    - sam file of all the reads aligned to just 1 genome
    - output file name
Output is a tsv that has all the reads that were split mapped to different scaffolds
"""

import sys
import pickle
import time
import re
from subprocess import call, DEVNULL

# samFile = call("samtools view -h " + sys.argv[1], shell=True)
# samFile = pysam.AlignmentFile(sys.argv[1], rb)
samFile1 = sys.argv[1]
samFile2 = sys.argv[2]
oneGenomeSam = sys.argv[3]
vcfFile = sys.argv[4]
outFile = sys.argv[5]

offsets = {'Dmel' : 0,
           'Dsim' : 222903,
           'Dsec' : 7736251,
           'flux' : 40000,
           'melflux' : 5000      
          }

mainChrom = '3L' ## for the VCF parsing

timeStart = time.clock()

def samSort(samFile, timeStart):
    """
    Takes a sam file handle (using bash, this can be a bam file from SAMTOOLS).
    Returns a dictionary of the split reads, with read names as keys
    """

    readDict = {}
    splits = {}
    count = 1
    for line in samFile:
        line = line.strip()

        ## take care of the header
        if (line[0] == '@'):
            continue
        
        ## split into read name, chromosome, and position
        line = line.split('\t')
        read = str(line[0])
        flag = int(line[1])
        chrom = line[2]
        pos = line[3]
        qual = int(line[4])
        
        ## check read quality
        if qual < 36: 
            continue
        
        ## check that it is not a secondary or a supp. read
        if len(bin(flag)) >= 11:
            ## secondary
            if bin(flag)[-9]:
                continue
            ## supplimentary
            if (len(bin(flag)) >= 14) and (bin(flag)[-12]):
                continue
        
        ## check to see if we have seen that read
        if read in readDict:
            ## check to see if it is a new scaffold
            if chrom not in readDict[read]['chroms']:
                ## add the new chrom and position, copy it to the return dictionary
                readDict[read]['chroms'].append(chrom)
                readDict[read]['positions'].append(pos)
                splits[read] = readDict[read]
        
        ## otherwise, enter it for the first time
        else:
            readDict[read] = {'chroms' : [chrom],
                              'positions' : [pos],
                              'melpos' : []}

        if not (count%1000000):
            print(str(count) + ' reads processed: ' + str(time.clock()-timeStart) + ' sec')
            # if count == 4000000:
            #     return splits
        count+=1
    
    return splits


def positionFilter(splits, offsets):
    """
    Use the offsets between the different chromosomes to determine if the
    reads are somewhat near each other
    """
    print("Filtering positions")
    popList = []

    for read in splits:
        dump = False
        pos = []
        for idx,spec in enumerate(splits[read]['chroms']):
            pair = (splits[read]['chroms'][idx], splits[read]['positions'][idx])
            pos.append(pair)
        for i in range(0, (len(pos)-1)):
            spec1 = pos[i][0].split('_')[0]
            spec2 = pos[i+1][0].split('_')[0]
            if len(pos) > 2:
                dump = True
            if (spec1 == 'Dmel') or (spec2 == 'Dmel'):
                dump = True
            upper = abs(offsets[spec1] - offsets[spec2] - offsets['flux'])
            lower = abs(offsets[spec1] - offsets[spec2] + offsets['flux'])
            if not upper > abs(int(pos[i][1]) - int(pos[i+1][1])) > lower:
                dump = True
        
        if dump:
            popList.append(read)

    for read in popList:
        splits.pop(read)
    return splits

def getReadInfo(splits, samFile):
    """Go back through the SAM file and find the read sequence of all the split reads
       I don't do this the first time because that data set would get way too big"""
    print("Re-reading for sequences")

    newSplits = {}
    for line in samFile:
        line = line.strip()
        ## take care of the header
        if (line[0] == '@'):
            continue
        
        ## split into read name, chromosome, and position
        line = line.split('\t')
        read = line[0]
        chrom = line[2]
        pos = line[3]
        seq = line[9]

        ## only work with the ones that we need
        if read in splits:


            ## if the new dictionary hasn't seen it
            if read not in newSplits:
                newSplits[read] = {}
            
            ## gather info
            readID = chrom + '_' + pos
            spec = chrom.split('_')[0]
            newSplits[read][readID] = [spec, seq]

    return newSplits


def checkAlignment(splits, oneGenomeSam, offsets):
    """Checks read positions against a single aligned file. Returns the splits dicitonay"""

    ## pop out anything that was too far
    print("Cross referencing with other file")
    print('Splits before cross reference: ' + str(len(splits.keys())))   
    
    ## go through the sam file, add positions
    for line in oneGenomeSam.readlines():
        line = line.strip()

        ## take care of the header
        if (line[0] == '@'):
            continue

        ## split into read name, chromosome, and position
        line = line.split('\t')
        read = line[0]
        qual = int(line[4])
        pos = int(line[3])       

        ## skip the loop if its not in there
        if read not in splits:
            continue
        if qual < 36:
            continue
        
        ## get the sequence
        seq = line[9]
        
        ## check to see what read this is
        for mate in splits[read]:
            if seq in splits[read][mate][1]:
                splits[read][mate].append(pos)



    ## see how far appart those positions are
    popList = []
    for read in splits:
        trash = False

        ## too many or too few enteries
        if len(splits[read].keys()) != 2:
            trash = True

        ## if they didn't all have mel positions
        for mate in splits[read]:
            if len(splits[read][mate]) != 3:
                trash = True

        ## too far appart
        if not trash:
            tally = []
            for mate in splits[read]:
                tally.append(int(splits[read][mate][2]))
            if abs(tally[0] - tally[1]) > offsets['melflux']:
                trash = True

        ## to the trash bin
        if trash:
            popList.append(read)
    

    for read in popList:
        splits.pop(read)
    print('Splits after cross reference: ' + str(len(splits.keys())))
    
    return splits

def isNear(vcfPos, pos):
    pos = int(pos)
    upper = pos + 150
    for idx in range(pos, upper):
        if idx in vcfPos:
            return True
    return False


def checkVCFpositions(splits, vcfFile):
    """Check that the reads are near a position where we could actually tell the difference"""
    ## make a list of all the positions where there could be a SNP
    vcfPos = {}
    
    print('Splits before SNP parsing: ' + str(len(splits.keys())))
    for line in vcfFile:
        line = line.strip()
        if 'CHRM' not in line and mainChrom in line:
            vcfPos[int(line.split('\t')[1])] = 0

    popList = []
    for read in splits:
        trash = False
        for mate in splits[read]:
            if not isNear(vcfPos, splits[read][mate][2]):
                trash = True
        if trash:
            popList.append(read)

    for read in popList:
        splits.pop(read)    
    print('Splits after SNP parsing: ' + str(len(splits.keys())))
    return splits


def printResults(splits, printOut, offsets):
    """
    Takes the output file handle
    Print the results
    Returns nothing
    """
    printOut.write('READname\tSPEC1\tPOS1\tSPEC2\tPOS2\n')  
    outList = []

    ## adjust all positions. However, this is not exact, only an estimate
    for read in splits:
        subL = []
        
        for mate in splits[read]:
            spec = splits[read][mate][0]
            pos = splits[read][mate][2]
            subL.append((spec,pos))
        
        # adjust
        # for i in range(0,2):
        #     spec = splits[read]['chroms'][i].split('_')[0]
        #     splits[read]['positions'][i] = int(splits[read]['positions'][i]) + offsets[spec]
        
        ## make a list of the two (spec,pos) pairs
        # subL.append((splits[read]['chroms'][0].split('_')[0], splits[read]['positions'][0]))
        # subL.append((splits[read]['chroms'][1].split('_')[0], splits[read]['positions'][1]))
        
        ## sort that list
        subL = sorted(subL, key=lambda x: x[1])
        
        ## add the read name at the end
        subL.append(read)
        outList.append(subL)

    outList = sorted(outList, key=lambda x: x[0][1], reverse = True)

    ## print to the file
    for read in outList:
        printOut.write(str(read[2]) + '\t' +
                       str(read[0][0]) + '\t' +
                       str(read[0][1]) + '\t' +
                       str(read[1][0]) + '\t' +
                       str(read[1][1]) + '\n')


## process the sam file
samFile1 = open(samFile1, 'r')
splits = samSort(samFile1, timeStart)
samFile1.close()

splits = positionFilter(splits, offsets)
samFile2 = open(samFile2, 'r')
splits = getReadInfo(splits, samFile2)
samFile2.close()

oneGenomeSam = open(oneGenomeSam, 'r')
splits = checkAlignment(splits, oneGenomeSam, offsets)
oneGenomeSam.close()

vcfFile = open(vcfFile, 'r')
splits = checkVCFpositions(splits, vcfFile)
vcfFile.close()


print('Writing results')
## print the results
printOut = open(outFile, 'w')
printResults(splits, printOut, offsets)
printOut.close()

# measure time
print(str(time.clock()-timeStart) + ' seconds to process sam file')