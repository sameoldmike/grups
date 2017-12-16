"""
usage: python PWD_from_stdin.py min_qual=int
                                [pileup=2sample_pileup_file.pileup] 
                                [verbose] 
                                [quiet] 
                                [min_depth=1,2] 
                                [known_variants]            # only count/allow bases that match a list of known REF/ALT alleles for each position
                                [samples=0,1] 
                                [targets=file.txt] 
                                [chr=int-int,int ...] 
                                [lines=int] 
                                [self_comparison] 
                                [filter_sites] 
                                [ignore_dels]
                                [block_JK]                  # turn on block jack-knife mode
                                [blocksize=int]             # block size for block jack-knife mode
"""




# import prerequisite packages
import sys
import numpy as np
import pylab as pl
import matplotlib
import gzip
from random import choice
import csv
import random
import argparse



# declare constants and options
mismatch_codes = {'A': 'A', \
                  'T': 'T', \
                  'C': 'C', \
                  'G': 'G', \
                  'N': 'N', \
                  'a': 'A', \
                  't': 'T', \
                  'c': 'C', \
                  'g': 'G', \
                  '*': '*', \
                  'n': 'N'}
transversions = [['C','A'], ['A','C'], ['G','T'], ['T','G']]
mismatches_1dir = [ set(['A','C']), set(['A','T']), set(['A','G']), \
                                    set(['C','T']), set(['C','G']), \
                                                    set(['T','G']), \
                    set(['A','N']), set(['C','N']), set(['T','N']), set(['G','N']), \
                    set(['A','*']), set(['C','*']), set(['T','*']), set(['G','*']) ]
mismatches_1dir_labels = ['A<>C', 'A<>T', 'A<>G', 'C<>T', 'C<>G', 'T<>G', 'A<>N', 'C<>N', 'T<>N', 'G<>N', 'A<>*', 'C<>*', 'T<>*', 'G<>*']

# need to update this to include indels
mismatches_2dir = [ ['A','C'], ['A','T'], ['A','G'], \
                    ['C','A'], ['C','T'], ['C','G'], \
                    ['T','A'], ['T','C'], ['T','G'], \
                    ['G','A'], ['G','C'], ['G','T'] ]   
mismatches_2dir_labels = ['A->C', 'A->T', 'A->G', 'C->A', 'C->T', 'C->G', 'T->A', 'T->C', 'T->G', 'G->A', 'G->C', 'G->T'] 
chrlist=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
chrlengths = [249250621, 243199373,  198022430,  191154276,  180915260,  171115067,  159138663,  146364022,  141213431,  135534747,  135006516,  133851895,  115169878,  107349540,  102531392,  90354753,   81195210,   78077248,   59128983,   63025520,   48129895,   51304566]



# parse command-line arguments
filename = ''
verbose = 0
quiet = 0
self_comparison = 0
min_depth = [0, 1]
known_variants = 0
min_qual = 0
samples = [0, 1]        # which samples to compare. default to pairwise comparison of both samples, but could be either sample against itself
targets_file = ''         
lines = 1000000
filter_sites = 0
ignore_dels = 0           
block_JK = 0            # block jack-knife mode OFF by default
blocksize = 10000000

parser = argparse.ArgumentParser()
parser.add_argument("--verbose", type=int, help="")
parser.add_argument("--quiet", type=int, help="")
parser.add_argument("--self_comparison", type=int, help="")
parser.add_argument("--known_variants", type=int, help="")
parser.add_argument("--filter_sites", type=int, help="")
parser.add_argument("--ignore_dels", type=int, help="")
parser.add_argument("--min_depth", type=str, help="")
parser.add_argument("--min_qual", type=str, help="")
parser.add_argument("--samples", type=str, help="")
parser.add_argument("--lines", type=str, help="")
parser.add_argument("--targets", type=str, help="")
parser.add_argument("--pileup", type=str, help="")
parser.add_argument("--chr", type=str, help="")
args = parser.parse_args()

if args.verbose:                                          
    verbose = 1
if args.quiet:   
    quiet = 1                                       
if args.self_comparison:
    self_comparison = 1
if args.known_variants:
    known_variants = 1
if args.filter_sites:
    filter_sites = 1
if args.ignore_dels:
    ignore_dels = 1
if args.min_depth:
    min_depth = args.min_depth.split(',')
    min_depth[0] = int(min_depth[0])
    min_depth[1] = int(min_depth[1])
if args.min_qual:
    min_qual = int(args.min_qual)
if args.samples:
    samples = args.samples.split(',')
    samples[0] = int(samples[0])
    samples[1] = int(samples[1])
if args.lines:
    lines = int(args.lines)
if args.targets:
    targets_file = args.targets
if args.pileup:
    filename = args.pileup
if args.chr:
    chrlist = args.chr
    chrlist = chrlist.split(',')
    new_chr = []
    for x in chrlist:
        x = x.split('-')
        if len(x) == 1:
            new_chr.append(int(x[0]))
        elif len(x) == 2:
            for z in range(int(x[0]), int(x[1])+1):
                new_chr.append(z)
        else:
            print 'error parsing chr= statement'
    chrlist = new_chr


# for this_arg in sys.argv[1:]:
#     #print this_arg
#     if this_arg.find('verbose') != -1:
#         verbose = 1
#         #print 'verbose = 1'
#     elif this_arg.find('quiet') != -1:
#         quiet = 1   
#         #print 'quiet = 1'
#     elif this_arg.find('self_comparison') != -1:
#         self_comparison = 1
#         #print 'self_comparison = 1'
#     elif this_arg.find('known_variants') != -1:
#         known_variants = 1
#     elif this_arg.find('filter_sites') != -1:
#         filter_sites = 1
#     elif this_arg.find('ignore_dels') != -1:
#         ignore_dels = 1
#     #elif this_arg.find('block_JK') != -1:
#     #    block_JK = 1        
#     #elif this_arg.find('blocksize=') != -1:
#     #    cmd_arg = 'blocksize='  
#     #    blocksize = int(this_arg[(this_arg.find(cmd_arg)+len(cmd_arg)):len(this_arg)])
#     elif this_arg.find('min_depth=') != -1:
#         cmd_arg = 'min_depth='  
#         min_depth = this_arg[(this_arg.find(cmd_arg)+len(cmd_arg)):len(this_arg)].split(',')
#         min_depth[0] = int(min_depth[0])
#         min_depth[1] = int(min_depth[1])
#     elif this_arg.find('min_qual=') != -1:
#         cmd_arg = 'min_qual='  
#         min_qual = int(this_arg[(this_arg.find(cmd_arg)+len(cmd_arg)):len(this_arg)])
#     elif this_arg.find('samples=') != -1:
#         cmd_arg = 'samples='  
#         samples = this_arg[(this_arg.find(cmd_arg)+len(cmd_arg)):len(this_arg)].split(',')
#         samples[0] = int(samples[0])
#         samples[1] = int(samples[1])
#     elif this_arg.find('lines=') != -1:
#         cmd_arg = 'lines='  
#         lines = int(this_arg[(this_arg.find(cmd_arg)+len(cmd_arg)):len(this_arg)])
#     elif this_arg.find('targets=') != -1:
#         cmd_arg = 'targets='  
#         targets_file = this_arg[(this_arg.find(cmd_arg)+len(cmd_arg)):len(this_arg)]
#     elif this_arg.find('pileup=') != -1:
#         cmd_arg = 'pileup='  
#         filename = this_arg[(this_arg.find(cmd_arg)+len(cmd_arg)):len(this_arg)]
#     elif this_arg.find('chr=') != -1:                                                   
#         cmd_arg = 'chr='
#         chrlist = this_arg[(this_arg.find(cmd_arg)+len(cmd_arg)):len(this_arg)]
#         chrlist = chrlist.split(',')
#         new_chr = []
#         for x in chrlist:
#             x = x.split('-')
#             if len(x) == 1:
#                 new_chr.append(int(x[0]))
#             elif len(x) == 2:
#                 for z in range(int(x[0]), int(x[1])+1):
#                     new_chr.append(z)
#             else:
#                 print 'error parsing chr= statement'
#         #print cmd_arg+this_arg[(this_arg.find(cmd_arg)+len(cmd_arg)):len(this_arg)]              
#         #print 'target contigs: ', new_chr
#         chrlist = new_chr
#         #chrlist = [x-1 for x in new_chr]
#         #random.shuffle(chrlist)
#     else:
#         print 'ERROR parsing command-line argument', this_arg




# ==========================================================================================================================================
# HASH IN TARGET SNP POSITIONS FILE (W/ ALT/REF ALLELE INFO IF INCLUDED)
# ==========================================================================================================================================
if targets_file != '':
    if verbose == 1:
        print 'hashing list of target positions'

    target_positions = {}
    input_file = open(targets_file, "rU")
    read_data = csv.reader(input_file, dialect=csv.excel_tab, lineterminator='\n', quoting=csv.QUOTE_NONE)
    for row in read_data:
        if len(row) == 2:
            target_positions[row[0]+'_'+row[1]] = 1
        elif len(row) == 4:
            target_positions[row[0]+'_'+row[1]] = [row[2], row[3]] # [alt, ref]
            
        else: 
            print 'ERROR parsing target SNP positions file'
    if verbose == 1:
        print '    # target positions hashed:\t'+str(len(target_positions.keys()))
    input_file.close()




# Calculate block edges using blocksize for block_JK mode
if block_JK == 1:
    block_edges = [None]*len(chrlengths)
    block_sitecounts = [None]*len(chrlengths)
    block_PWDcounts = [None]*len(chrlengths)
    for chromosome in range(len(chrlengths)):
        block_edges[chromosome] = []
        block_sitecounts[chromosome] = []
        block_PWDcounts[chromosome] = []

        for X in range(1, chrlengths[chromosome]+1, blocksize):
            if X != 1:
                block_edges[chromosome].append(X)
                block_sitecounts[chromosome].append(0)
                block_PWDcounts[chromosome].append(0)
        if block_edges[chromosome][-1] < chrlengths[chromosome]:
            block_edges[chromosome].append(chrlengths[chromosome])
            block_sitecounts[chromosome].append(0)
            block_PWDcounts[chromosome].append(0)
        if verbose == 1:
            print 'block_edges[chromosome] =', block_edges[chromosome]



# ==========================================================================================================================================
# READ IN PILEUP FILE
# ==========================================================================================================================================

    # pileup file format
    # --------------------------------------------------------------------------------------------------------------------------------------------------------
    # 0                                                 1                   2           3           4           5           6           7           8
    # chromosome            	                        1-based position    reference   bam1_depth  bam1_bases  bam1_quals  bam2_depth  bam2_bases  bam2_quals
    # supercont1.4920_of_Phytophthora_infestans_T30-4	297	                C	        1	        ,	        e           1           ,           e           
    # supercont1.4920_of_Phytophthora_infestans_T30-4	298	                A	        1	        ,	        e           1           ,           e
    # supercont1.4920_of_Phytophthora_infestans_T30-4	299	                C	        1	        ,	        e           1           ,           e
    # supercont1.4920_of_Phytophthora_infestans_T30-4	300	                A	        1	        ,	        e           1           ,           e
    # --------------------------------------------------------------------------------------------------------------------------------------------------------



# list of positions to output from unordered hash table to ordered FASTA
#hap_hash = {}
#list_of_positions = []

pileup_obs = {}

# open pileup
#input_file = gzip.open(filename)
#pileup_data = csv.reader(input_file, dialect=csv.excel_tab, lineterminator='\n', quoting=csv.QUOTE_NONE)

lines_read = 0
overlap_sites = 0
PWD = 0
all_quals = 0
mismatches_1dir_counts = [0]*len(mismatches_1dir)
mismatches_2dir_counts = [0]*len(mismatches_2dir)
filtered_triallelic_sites = 0
chr_labels = []
chr_sitecounts = []
chr_PWDcounts = []
#for row in pileup_data:
for line in sys.stdin:
    if verbose == 1:
        print line
    line = line.rstrip('\n')
    row = line.split('\t')
    #print row
    if quiet == 0:
        if lines_read == 0:
            print 'ex:\t', row
        if lines_read % lines == 0 and lines_read != 0:
            if overlap_sites != 0:
                print 'pileup lines read:', lines_read, '\t# PWDs:', PWD, '\tpost-filter sites:', overlap_sites, '\tPWD/site:', str(float(PWD)/overlap_sites)
            else:
                print 'pileup lines read:', lines_read, '\tpost-filter sites: 0'

    # save the read depth for this SNP
    #all_depths = all_depths + (int(row[3]),)
    #all_depths.append(int(row[3]))

    # decide if this genome site should be skipped
    skip_site = 1
    if row[0].isdigit():
        # if the contig ID is an autosomal chromosome
        if int(row[0]) >= 1 and int(row[0]) <= 22:        
            # if site is in a target chromosome
            if int(row[0]) in chrlist:
                skip_site = 0
                # if a target positions list was provided
                if targets_file != '':
                    # check if site is in target positions list
                    if row[0]+'_'+row[1] in target_positions:
                        skip_site = 0
                        ## if filtering transitions
                        #if only_tv == 1:
                        #    # if this position's hashed ALT/REF alleles represent a transition
                        #    if target_positions[row[0]+'_'+row[1]] not in transversions:
                        #        # skip this site
                        #        skip_site = 1
                    else:
                        skip_site = 1
            # if the pileup contig is after any of the target contigs, stop parsing the pileup data
            elif int(row[0]) > max(chrlist)+1:
                break
    
    # if neither sample is completely missing data at this site AND this site should not be ignored
    if row[4] != '*' and row[7] != '*' and skip_site == 0:
        #if verbose == 1:
        #    print line

        min_depth1 = row[3]
        min_depth2 = row[6]
        found_indel = 0
                
        if min_depth1 >= min_depth[0] and min_depth2 >= min_depth[1]:
            # read in the pileup data one sample at a time
            pileup_bases = [[], []]
            pileup_quals = [[], []]
            
            for K in [0, 1]:
        
                # read the base pileup field and base qual fields
                if samples[K] == 0:
                    pileup_line = row[4]
                    quals_line = row[5]
                elif samples[K] == 1:
                    pileup_line = row[7]
                    quals_line = row[8]
                
                line_pos = 0
                qual_pos = 0
                while line_pos < len(pileup_line):
                    # if base matches reference
                    if pileup_line[line_pos] == '.' or pileup_line[line_pos] == ',':
                        # if the reference base is not available in the pileup file
                        if row[2] == 'N':
                            # use the str 'r' to represent the reference base
                            pileup_bases[K].append('r')
                        else:
                            pileup_bases[K].append(row[2])
                        pileup_quals[K].append(int(ord(quals_line[qual_pos])-33))
                        #print int(ord(quals_line[qual_pos])-33)
                        line_pos += 1
                        qual_pos += 1

                    # if base is a mismatch    
                    elif int(pileup_line[line_pos] in mismatch_codes) == 1:
                        # if this is a deleted position and ignore_dels option is ON
                        if ignore_dels == 1 and pileup_line[line_pos] == '*':
                            pass
                        else:
                            pileup_bases[K].append(mismatch_codes[pileup_line[line_pos]])
                            pileup_quals[K].append(int(ord(quals_line[qual_pos])-33))
                            #print int(ord(quals_line[qual_pos])-33)
                        line_pos += 1
                        qual_pos += 1                

                    # if a read starts here
                    elif pileup_line[line_pos] == '^':
                        line_pos += 2
                        qual_pos += 0
 
                    # if a read ends here
                    elif pileup_line[line_pos] == '$':
                        line_pos += 1
                        qual_pos += 0

                    # if an INDEL starts here    
                    elif pileup_line[line_pos] == '+' or pileup_line[line_pos] == '-':
                        found_indel = 1
                    
                        # find the number of characters encoding the indel size integer
                        looker = line_pos+1
                        while int(pileup_line[looker].isdigit()) == 1:
                            looker += 1
                        indel_encode_digits = looker - line_pos - 1

                        # count this indel as an uncalled base
                        #pileup_bases[K].append('N')
                        indel_size = int(pileup_line[(line_pos+1):(line_pos+1+indel_encode_digits)])                
                        line_pos += (1 + indel_encode_digits + indel_size)
                        qual_pos += 0

                    elif pileup_line[line_pos] == '>' or pileup_line[line_pos] == '<':
                        print 'ERROR. must parse > or < in pileup bases column.'
                    else:
                        print 'ERROR!'
            if len(pileup_bases[0]) != len(pileup_quals[0]):
                print 'ERROR!'
            if len(pileup_bases[1]) != len(pileup_quals[1]):
                print 'ERROR!'
            
            #if verbose == 1 and found_indel == 1:
            #    print '---INDEL---'
            #    print row
            #    print ' col0_pileup:\t', pileup_bases[0], pileup_quals[0]
            #    print ' col1_pileup:\t', pileup_bases[1], pileup_quals[1]
            #    print '-----------'

            #if verbose == 1:
            #    print pileup_bases[0], pileup_quals[0], pileup_bases[1], pileup_quals[1]
        
            # filter bases by quality score (< min_qual) and, if known_variants == 1, if they are known variants
            base_choices1 = []
            base_choices2 = []
            qual_choices1 = []
            qual_choices2 = []
            for t in range(0, len(pileup_bases[0])):
                exclude_base = 0
                if pileup_quals[0][t] >= min_qual:
                    if known_variants == 1:
                        if pileup_bases[0][t] not in target_positions[row[0]+'_'+row[1]]:
                            exclude_base = 1
                    if exclude_base == 0:
                        base_choices1.append(pileup_bases[0][t])
                        qual_choices1.append(pileup_quals[0][t])
            for t in range(0, len(pileup_bases[1])):
                exclude_base = 0
                if pileup_quals[1][t] >= min_qual:
                    if known_variants == 1:
                        if pileup_bases[1][t] not in target_positions[row[0]+'_'+row[1]]:
                            exclude_base = 1
                    if exclude_base == 0:
                        base_choices2.append(pileup_bases[1][t])
                        qual_choices2.append(pileup_quals[1][t])
                                    
            # if number of bases to choose from is >= min_depth in both samples
            if len(base_choices1) >= min_depth[0] and len(base_choices2) >= min_depth[1]:

                if verbose == 1:
                    # print some quality control and debugging data
                    for cc in range(0, len(pileup_bases[0])):
                        if (pileup_bases[0][cc] == 'N' or pileup_bases[0][cc] == 'n') and pileup_quals[0][cc] >= 2:
                            print row, str(pileup_quals[0][cc])
                    for cc in range(0, len(pileup_bases[1])):
                        if (pileup_bases[1][cc] == 'N' or pileup_bases[1][cc] == 'n') and pileup_quals[1][cc] >= 2:
                            print row, str(pileup_quals[1][cc])
                    print pileup_bases[0], '...', pileup_bases[1]
                    print pileup_quals[0], '...', pileup_quals[1]
                    print base_choices1, '...', base_choices2
                    print '---'

            
                # if in self_comparison mode
                if self_comparison == 1:
                    valid_base_choices = []
                    valid_base_quals = []
                    for xxx in range(0, len(base_choices1)-1):
                        for yyy in range(xxx+1, len(base_choices2)):
                            # iterate through and save all possible pairwise quality-filtered base comparisons
                            choice1 = base_choices1[xxx]
                            choice2 = base_choices2[yyy]
                            qual1 = qual_choices1[xxx]
                            qual2 = qual_choices2[yyy]
                        
                            valid_base_choices.append([choice1, choice2])
                            valid_base_quals.append([qual1, qual2])
                                        
                    overlap_sites += 1
                
                    # if filter_sites mode is on, then simply print the contig and position
                    if filter_sites == 1:
                        print line
                        
                    # if filter_sites mode is off, then perform the PWD comparison
                    else:                    
                        # count the number of sites per chr
                        if row[0] in chr_labels:
                            chr_sitecounts[chr_labels.index(row[0])] += 1 
                        else:
                            chr_labels.append(row[0])
                            chr_sitecounts.append(1)
                            chr_PWDcounts.append(0)
        
                        # update block jack-knife stats if necessary
                        if block_JK == 1:
                            # find current chr and pos
                            current_chr = int(row[0])-1
                            current_pos = int(row[1])
                            # find current block                        
                            for current_block in range(len(block_edges[current_chr])):
                                if current_pos <= block_edges[current_chr]:
                                    break                              
                            # update count of # overlap sites within block  
                            block_sitecounts[current_chr][current_block] += 1
                            
                        # this chooses one of the valid base comparisons at this position and checks if it is a PWD
                        h = random.choice(range(0, len(valid_base_choices)))
                        all_quals += (float(valid_base_quals[h][0] + valid_base_quals[h][1])/2)
                        if valid_base_choices[h][0] != valid_base_choices[h][1]:
                            PWD += 1
                            chr_PWDcounts[chr_labels.index(row[0])] += 1
                            # count the mismatch type
                            mismatches_1dir_counts[mismatches_1dir.index(set([valid_base_choices[h][0], valid_base_choices[h][1]]))] += 1.0
                            
                            if block_JK == 1:
                                # update block-specific mismatch count
                                block_PWDcounts[current_chr][current_block] += 1
                            
                            if verbose == 1:
                                print '---PWD---'
                        elif verbose == 1:
                                print '---No PWD---'

                        if verbose == 1:
                            print line
                            print base_choices1, base_choices2
                            print valid_base_choices
                            print row[0], row[1], 'ref:'+row[2], valid_base_choices[h][0], valid_base_choices[h][1]
                            print valid_base_quals[h][0], '\t', valid_base_quals[h][1]

                        
                # otherwise, if not performing a self-comparison
                else:                    
                    # find if a mismatch or not by randomly comparing a single base from each sample    
                    choice1 = random.choice(range(0, len(base_choices1)))
                    choice2 = random.choice(range(0, len(base_choices2)))
                    qual1 = qual_choices1[choice1]
                    qual2 = qual_choices2[choice2]
                    choice1 = base_choices1[choice1]
                    choice2 = base_choices2[choice2]
            
                    count_site = 1
            
                    # if this site should be included in the analysis (after filtering out tri-allelic observations and transitions, if applicable)
                    if count_site == 1:                            
                        overlap_sites += 1

                        # if filter_sites mode is on, then simply print the contig and position
                        if filter_sites == 1:
                            print line
                        
                        # if filter_sites mode is off, then perform the PWD comparison
                        else:
                            # count the number of sites per chr
                            if row[0] in chr_labels:
                                chr_sitecounts[chr_labels.index(row[0])] += 1 
                            else:
                                chr_labels.append(row[0])
                                chr_sitecounts.append(1)
                                chr_PWDcounts.append(0)
                
                            # update block jack-knife stats if necessary
                            if block_JK == 1:
                                # find current chr and pos
                                current_chr = int(row[0])-1
                                current_pos = int(row[1])
                                # find current block                        
                                for current_block in range(len(block_edges[current_chr])):
                                    if current_pos <= block_edges[current_chr][current_block]:
                                        break                              
                                # update count of # overlap sites within block  
                                block_sitecounts[current_chr][current_block] += 1
                               
                                if verbose == 1:
                                    print 'current_chr =', current_chr
                                    print 'current_pos =', current_pos
                                    print 'current_block =', current_block
                                    print 'block_edges[current_chr] =', block_edges[current_chr]
                                    print 'block_sitecounts[current_chr][current_block] =', block_sitecounts[current_chr][current_block]

                            # calculate and save statistics of pairwise differences
                            all_quals = all_quals + (float(qual1 + qual2)/2.0)
                            if choice1 != choice2:
                                PWD += 1
                                chr_PWDcounts[chr_labels.index(row[0])] += 1 
                                # count the mismatch type
                                mismatches_1dir_counts[mismatches_1dir.index(set([choice1, choice2]))] += 1
                                
                                if block_JK == 1:
                                    # update block-specific mismatch count
                                    block_PWDcounts[current_chr][current_block] += 1

                                if verbose == 1:
                                    print '---PWD---'
                            elif verbose == 1:
                                    print '---No PWD---'

                            if verbose ==1:
                                print pileup_bases[0], '\t', pileup_bases[1], '\tref: ', row[2]+':' 
                                print pileup_quals[0], '\t', pileup_quals[1]
                                print base_choices1, '\t', base_choices2
                                print str(row[0])+'_'+str(row[1])+'\tref:'+row[2]+'\t'+choice1+'\t'+choice2+'\t'+str(qual1)+'\t'+str(qual2)
                                    #print '---------'

    lines_read += 1



if filter_sites == 0:
    # calculate PWD / site 
    print ''
    if block_JK == 0:    
        if quiet == 0:
            print '#PWD/site\t#PWD\t#overlap_sites\tmean_basequal\tfiltered_triallelic_sites'
        if overlap_sites != 0:
            print str(float(PWD)/overlap_sites)+'\t'+str(PWD)+'\t'+str(overlap_sites)+'\t'+str(float(all_quals)/overlap_sites)+'\t'+str(filtered_triallelic_sites)
        else:
            print 'NA'+'\t'+str(PWD)+'\t'+str(overlap_sites)+'\t'+'NA'
        if quiet == 0:
            # print mean PWD/site stats
            if block_JK == 0:
                # calculate per-chr PWD/site
                print 'contig\tPWD/site\tPWD\toverlap_sites'
                for t in range(0, len(chr_labels)):
                    print chr_labels[t]+'\t'+str(float(chr_PWDcounts[t])/chr_sitecounts[t])+'\t'+str(chr_PWDcounts[t])+'\t'+str(chr_sitecounts[t])

                # print histogram of mismatch types 
                for t in mismatches_1dir_labels:
                    print t+'\t',
                print ''
                mismatches_1dir_fracs = []
                for t in mismatches_1dir_counts:
                    mismatches_1dir_fracs.append(t/float(PWD))
                for t in mismatches_1dir_fracs:
                    print str(t)+'\t',
                print ''    
    
    # print block jack-knife mean PWD/site stats
    elif block_JK == 1:
        # calculate per-block PWD/site
        if quiet == 0:
            print 'contig\tPWD/site\tPWD\toverlap_sites'
        for t in chr_labels:
            this_chrlabel = int(t)-1
            for block in range(len(block_sitecounts[this_chrlabel])):
                if block_sitecounts[this_chrlabel][block] != 0:
                    print str(t)+'\t'+str(float(block_PWDcounts[this_chrlabel][block])/block_sitecounts[this_chrlabel][block])+'\t'+str(block_PWDcounts[this_chrlabel][block])+'\t'+str(block_sitecounts[this_chrlabel][block])
                else:
                    print str(t)+'\t'+'NA'+'\t'+str(block_PWDcounts[this_chrlabel][block])+'\t'+str(block_sitecounts[this_chrlabel][block])
