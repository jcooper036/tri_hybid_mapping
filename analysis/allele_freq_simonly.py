#!/usr/bin/env python3
"""
The point of this script is to find the positions in a VCF file that are tripartite
between the different samples, and then make a new parsed table of those SNPs with
their identies. The output is a tsv of the SNP position and all the allele freq
for each genotype at that postions.

Warning: This script requires that 

Input: VCF file, output file tag (will make indivdual files for all sequences using this header)

"""

import time
import sys

conThresh = 80 # confidence threshold for read calling
outfile = sys.argv[2]

## has to be samples in order that they will appear in the vcf
allSamples = ['J10075', 'J10224','J20039',
              'J40022','J40023','J40024',
              'rep1_female','rep1_male','rep2_female',
              'rep2_male','rep3_female','rep3_male']

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

def listWrite(lt):
    """Takes a list and returns the first element as a string"""
    return str(lt[0])


def outWrite(outfile, alleleTable, sample):
    """Writes the list to the output file. Each line should already be in a nicely
       writable format"""
    with open ((outfile + '_' + sample + '.tsv'), 'w') as f:
        f.write('CHROM\tPOS\tmel\tsim\tsec\n')
        for line in alleleTable:
            f.write('\t'.join(map(str,line)))
            f.write('\n')

### some counters to track things
tries = 0
madeit = 0

## make a dictionary that has a lists for all the sample names
alleleTables = {}
for smpName in allSamples:
    alleleTables[smpName] = []

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
            ## stats that goes [J10075, J10224, J20039, J40022, J40023, J40024, rep1_female, rep1_male, rep2_female, rep2_male, rep3_female, rep3_male]
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
            
            ## next we want to know if the species are segregated int types
            specSeg = False
            if qScore and allHom:
                melTypes = [n.split(':')[0].split('/')[0].replace('.', '0') for n in spec['mel']]
                simTypes = [n.split(':')[0].split('/')[0].replace('.', '0') for n in spec['sim']]
                secTypes = [n.split(':')[0].split('/')[0].replace('.', '0') for n in spec['sec']]
                
                ## if there are 3 unique SNPs at that postition
                for var in melTypes:
                    if var not in simTypes and var not in secTypes:
                        for var in simTypes:
                            if var not in secTypes:
                                specSeg = True
            
            ## if we should save them, figure out how to do so
            if specSeg:                
                madeit += 1
                numlist = [int(melTypes[0]), int(simTypes[0]), int(secTypes[0])]  #for sorting, since the allele numbers will be given in the order 0, 1, 2

                ## now, take the sample info from each line, add it to the allele table
                for idx, alChart in enumerate(info):
                    sample = alChart
                    if sample != './.': # just don't use these
                        sample = sample.split(':')[1].split(',')
                        
                        ## turn the counts into frequencies
                        sample = [int(u) for u in sample]

                        tots = sum(sample)
                        sample = [(u/tots) for u in sample]
                        
                        ## assign the variables to all the things
                        chrm = str(line.split('\t')[0])
                        bpidx = str(line.split('\t')[1])
                        melF = sample[numlist[0]]
                        simF = sample[numlist[1]]
                        secF = sample[numlist[2]]
                        if simF > 0 or secF > 0:
                            tot = simF + secF
                            simF /= tot
                            secF /= tot
                            melF = 0
                            alleleTables[allSamples[idx]].append([chrm, bpidx, melF, simF, secF])

        if "#CHROM" in linex:
            # check to start reading. This comes after everything
            # else so we will start reading AFTER it procs
            headerDone = True

for sample in alleleTables:
    outWrite(outfile, alleleTables[sample], sample)

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
