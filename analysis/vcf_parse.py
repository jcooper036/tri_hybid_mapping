#!/usr/bin/env python3
"""
The point of this script is to find the positions in a VCF file that are tripartite
between the different samples, and then make a new parsed table of those SNPs with
their identies. The downstream application is to comare these to BAM files and 
Needs the vcf as an input after the script is invoked.

Input: VCF file, output file

"""

import time
import sys

conThresh = 80 # confidence threshold for read calling
outfile = sys.argv[2]

x = time.clock()

## we are looking for postions that are fixed between the different species
## then we want to build a table that says
#CHR #POS #MEL #SIM #SEC
# 2L XXXX    A    T    G

## where MEL, SIM, and SEC columns will have the identity of the base at that
## postion.Then use the table to go through reads and look at positions within 
# those reads to determine what parent they came from

def findNT(line, varNum):
    """Takes a line from a GATK file, and the varient number, and returns the
       nucleotide that corresponds to that variant ID"""
    # nts is a list where 0 is the reference, and then each postion after that
    # is the postion of the variable
    nts = []
    nts.append(line.split('\t')[3])
    seconds = line.split('\t')[4].split(',')
    for var in seconds:
        nts.append(var)
    if varNum == '.':
        varNum = 0
    return nts[int(varNum)]

def outWrite(outfile, varList):
    """Writes the list to the output file. Each line should already be in a nicely
       writable format"""
    with open (outfile, 'w') as f:
        for line in varList:
            f.write(line + '\n')


def all_same(items):
    return all(x == items[0] for x in items)

### some counters to track things
tries = 0
madeit = 0

# open the file
with open (sys.argv[1], "r") as vcfFile:
    headerDone = False # make sure we don't output anything from the header
    varList = ['CHRM\tPOS\tMEL\tSIM\tSEC']
    for linex in vcfFile:
        tries += 1
        spec = {}
        if headerDone:
            line = linex.strip()
            info = linex.strip()
            info = info.split('\t')[9:]
            
            ## at this point, there is a list of all the info
            ## stats that goes [J10075, J10224, J20039, J40022, J40023, J40024]
            spec['mel'] = [info[0], info[1]]
            spec['sim'] = [info[2]]
            spec['sec'] = [info[3], info[4], info[5]]

            ## first thing to do is to check that everything passes the quality score
            ## next, make sure that all the samples are homozygous for the call

            qScore = True
            allHom = True
            for ty in spec:
                for sample in spec[ty]:
                    #qscore
                    if sample.split(':')[0] != "./.":
                        if float(sample.split(':')[3]) < conThresh:
                            qScore = False 
                        
                        #homozygous
                        y = sample.split(':')[0]
                        y = y.split('/')
                        if y[0] != y[1]:
                            allHom = False
                        
                        ## filter for anything where any other allele was read
                        counts = sample.split(':')[1].split(',')
                        nonZero = False
                        for csa in counts:
                            if nonZero and int(csa) != 0:
                                allHom = False
                            if int(csa) != 0:
                                nonZero = True
                            
            
            ## next we want to know if the species are segregated int types
            specSeg = False
            if qScore and allHom:
                melTypes = [n.split(':')[0].split('/')[0].replace('.', '0') for n in spec['mel']]
                simTypes = [n.split(':')[0].split('/')[0].replace('.', '0') for n in spec['sim']]
                secTypes = [n.split(':')[0].split('/')[0].replace('.', '0') for n in spec['sec']]

                ## make sure they are all the same in a species
                if all_same(melTypes):
                    melType = melTypes[0]
                else:
                    continue
                
                if all_same(simTypes):
                    simType = simTypes[0]
                else:
                    continue
                
                if all_same(secTypes):
                    secType = secTypes[0]
                else:
                    continue

                ## if any one of them is unique
                if (melType != simType) or (melType != secType) or (simType != secType):
                    specSeg = True
            
            ## if we should save them, figure out how to do so
            if specSeg:                
                madeit += 1
                melNT = ''
                simNT = ''
                secNT = ''
                
                ## might be the reference
                melNT = findNT(line, melType)
                simNT = findNT(line, simType)
                secNT = findNT(line, secType)

                ## add all the good variants to the list
                varList.append(str(line.split('\t')[0]) + '\t' +
                               str(line.split('\t')[1]) + '\t' +
                               melNT + '\t' +
                               simNT + '\t' +
                               secNT)

        if "#CHROM" in linex:
            # check to start reading. This comes after everything
            # else so we will start reading AFTER it procs
            headerDone = True

outWrite(outfile, varList)

print(str(tries) + ' variants checked, ' + 
      str(madeit) + ' informative SNPs (' +
      str(100*madeit/tries) + ' percent)')
print(str(time.clock()-x) + ' seconds to parse VCF')


## notes: GT:AD:DP:GQ:PL
    # GT - genotype at that position for the two haplotypes 
        # 0/1 = het for the alt allele
        # 1/1 = homo for alt allele 1
        # 2/2 = homo for alt allele 2
    # AD - total unfiltered allele depth for each allele
        # but does not include uniformative reads
    # DP - filted allele depth, total for the postion
        # does include uniformative reads
    # GQ - quality of the assigned genotype
        # high = large difference between likelihood of this and 2nd most likely
        # low = little difference, less reason to pick on over the other
    # PL - normalized Phred score for all genotypes
