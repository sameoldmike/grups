# import prereq packages
import numpy as np
import random
import gzip
import pysam
#from bisect import bisect_left
import code #code.interact(local=locals())
import csv
#import bx
from bx.intervals.intersection import Intersecter, Interval


# define constants
CHROM,POS = 0,1
#seqerror_choices = {'A': ['C','T','G'], 'C': ['A','T','G'], 'T': ['A','C','G'], 'G': ['A','C','T']}
seqerror_choices = {0: [1, 2, 3], 1: [0, 2, 3], 2: [0, 1, 3], 3: [0, 1, 2]}
mismatch_codes = {'A': 'A', \
                  'T': 'T', \
                  'C': 'C', \
                  'G': 'G', \
                  'N': 'N', \
                  '*': '*', \
                  'a': 'A', \
                  't': 'T', \
                  'c': 'C', \
                  'g': 'G', \
                  'n': 'N'}

# define class
class Individual:
    def __init__(self, ID):
        self.chromosomeLength = 0
        self.ID = ID
        self.parent1 = ''
        self.parent2 = ''
        self.hap1 = []
        self.hap2 = []
        
    def __del__(self):
        del self.chromosomeLength
        del self.ID
        del self.parent1
        del self.parent2
        del self.hap1
        del self.hap2

    def update_default(self, size, brin1, brin2):
        self.hap1=np.array([brin1]*size, dtype=np.int8)
        self.hap2=np.array([brin2]*size, dtype=np.int8)

    def _update_1site_from1000g(self, d, p):
        self.hap1[p]=d[0]
        self.hap2[p]=d[2]

    def update(self, haplo1, haplo2, chromsize):
        self.chromosomeLength = chromsize
        self.hap1 = haplo1.copy()
        self.hap2 = haplo2.copy()

    def show(self):
        """
        Display both haplotypes
        """
        print self.hap1
        print self.hap2

    def meiosis(self, theta, chromosome_map, recomb_rate_regions, recomb_probs):
        """
        Simulate recombination at rate theta
        """
        
       # define recombination positions variables
        recombPos = []
        recomb_occurred = np.zeros(len(self.hap1), dtype=np.int)

        #  if using uniform recombination rate
        if len(recomb_rate_regions) == 0:   
            # define uniform per-base probability of recombination using Haldane's model of recombination (1919?)
            recomb_prob = (1 - np.exp(-2*theta))/2.0

            for y in range(1, len(chromosome_map)):
                interval_prob_recomb = recomb_prob*(chromosome_map[y] - chromosome_map[y-1])    # CHECK IF THIS PROBABILITY IS ALWAYS POSITIVE !!
                if np.random.random() < interval_prob_recomb:
                    recombPos.append(y)
        
        else:    
            # define an interval tree from recomb_rate_regions
            tree = Intersecter()
            for start, end in recomb_rate_regions:
                tree.add_interval(Interval(start, end))
                
            # define query intervals from chromosome_map 
            chromosome_map_pairs = [None]*(len(chromosome_map)-1)
            for y in range(1, len(chromosome_map)):
                chromosome_map_pairs[y-1] = (chromosome_map[y-1], chromosome_map[y])
            
            # find recomb_rate_regions intervals overlapping chromosome_map_pairs intervals
            overlaps = [None]*len(chromosome_map_pairs)
            pair_counter = 1
            total_prob = 0
            for start, end in chromosome_map_pairs:
                overlaps[pair_counter-1] = tree.find(start, end)
                #print '            chromosome_map_pairs interval matches: (%s, %s) -> %s' % (start, end, overlaps[-1])
                interval_prob_recomb = 0
                if len(overlaps[pair_counter-1]) > 0:
                    this_index = 0
                    for this_overlap in overlaps[pair_counter-1]:
                        #print 'this_overlap', this_overlap
                        # calculate total amount overlap with this region
                        if start < this_overlap.start:
                            real_start = this_overlap.start
                        else: 
                            real_start = start
                        if end > this_overlap.end:
                            real_end = this_overlap.end
                        else:
                            real_end = end
                        # calculate contribution to total probability of recombination in this region
                        this_interval_recomb_prob = recomb_probs[this_index]*(real_end-real_start+1)
                        interval_prob_recomb += this_interval_recomb_prob
                        total_prob += this_interval_recomb_prob
                        this_index += 1
                    
                if np.random.random() < interval_prob_recomb:
                    recombPos.append(pair_counter)
                pair_counter += 1 
            print '            total prob of recombination on this chr =', total_prob
        
        # print debugging stuff
        # for g in  range(0, len(recombPos)):
        #     print 'g =', g
        #     print 'recombPos[g] =', recombPos[g]
        #     print 'bisect_left(chromosome_map, recombPos[g]) =', bisect_left(chromosome_map, recombPos[g])
        #     print 'chromosome_map[bisect_left(chromosome_map, recombPos[g])] =', chromosome_map[bisect_left(chromosome_map, recombPos[g])]
        #     print 'recomb_occurred[bisect_left(chromosome_map, recombPos[g])] =', recomb_occurred[bisect_left(chromosome_map, recombPos[g])]
        #     print '---'


        if len(recombPos) != 0:
            for g in recombPos:
                recomb_occurred[g] = 1
                print '            Recombination at '+str(chromosome_map[g])+'. chromosome_map position: ['+str(g)+'/'+str(len(chromosome_map))+']'
        else:
            print '            No recombination to simulate'
        
        recombinants = recomb_occurred.cumsum() % 2
        
        original_hap1   = self.hap1.copy() 
        original_hap2   = self.hap2.copy() 
        recombined_hap1 = self.hap1.copy()
        recombined_hap2 = self.hap2.copy()
        
        if len(recombPos) != 0:
            recombined_hap1[recombinants == 1] = original_hap2[recombinants == 1]
            recombined_hap2[recombinants == 1] = original_hap1[recombinants == 1]
        
        return recombined_hap1, recombined_hap2


# END class Individual



# define functions
def update_individuals_from1000g(chr, name, sizeChrom, numvars, colInDatafile, listindiv, downsample_rate_variation, downsample_rate_numSNPs, pedigree_pop_ID_index, popAF_polymorph_filter, pileup_positions, use_pileup_positions, target_positions, use_target_positions, pedigree_pop, contam_pop, contam_ind_IDs, min_AF):
    """
    Update a list of objects of class Individual with genotype data from datafile 'name' (format type: vcf, eg 1000genomes)
    """
    
    nind = len(listindiv)
    if len(colInDatafile) != nind:
        print "Error in create_individuals_from1000g"; return

    pos_lookup = {}
    pos_rev_lookup = {}
    pos_vector_length = 0

    # create empty haplotype vectors
    haplo1 = [np.zeros(numvars, dtype=np.int8) for x in xrange(nind)]  # independant np arrays
    haplo2 = [np.zeros(numvars, dtype=np.int8) for x in xrange(nind)]  # independant np arrays

    # create empty ordered lists of allele frequencies in pedigree superpop and contaminating superpop
    pedigree_pop_AF = []
    contam_pop_AF   = [[], []]
    
    lines_read = 0
    fileop = pysam.Tabixfile(name)
    for line in fileop.fetch(str(chr+1)):       # compatible with pysam 0.8, but clunky as chr MUST be correct
    #for line in fileop.fetch():                # compatible with pysam 0.6
        lines_read += 1   
        #if (lines_read % 1000000 == 0) and (lines_read != 0):
        #    print 'read '+str(lines_read/1000000)+'M lines ...'
        
        #print line
        elem = line.strip('\n').split('\t')

        # perform initial filtering of variants
        include_SNP = 1
        
        # if a pileup and/or targets file with target positions has been provided [target positions file may contain a subset of pileup positions, but a position must be in both lists to be included]
        if use_pileup_positions == 1:
            # check if this BAM position is in the list of pileup targets
            if int(elem[1])-1 not in pileup_positions:
                # immediately filter out this position
                include_SNP = 0
        if use_target_positions == 1:
            # check if this BAM position is in the list of targets positions 
            if int(elem[1])-1 not in target_positions:
                include_SNP = 0

        # perform normal filtering on this position if it's not already filtered out
        if include_SNP == 1:
             
            # if this variant is not biallelic
            if len(elem[4].split(',')) != 1:
                # ignore variant
                include_SNP = 0
            
            # if this variant is not a SNP
            elif len(elem[3]) != 1 or len(elem[4]) != 1:
                # ignore variant
                include_SNP = 0

            else:
                #include_SNP = 1
                pos_on_chromosome = int(elem[POS])-1
        
                # count the number of different alleles present in the EUR superpop
                alleles = []
                for ind in pedigree_pop_ID_index:
                    alleles.append(elem[ind][0])
                    alleles.append(elem[ind][2])
            
                # filter (or not) sites that are monomorphic/polymorphic in EUR population
                if popAF_polymorph_filter == 'all':
                    if len(set(alleles)) == 0 or len(set(alleles)) >= 3:
                        include_SNP = 0
                else:
                    if popAF_polymorph_filter == 'mono':
                        if len(set(alleles)) != 1:
                            include_SNP = 0
                    if popAF_polymorph_filter == 'poly':
                        if len(set(alleles)) != 2:
                            include_SNP = 0
            
        # perform additional filtering if SNP is still not filtered out
        if include_SNP == 1:
            INFO = elem[7].split(';')
            
            # extract actual alternate allele frequency at this SNP for both pedigree superpopulation and contaminating superpopulation
            for t in INFO:
                if t.find(pedigree_pop+'_AF=') == 0:
                    this_pedigree_pop_AF = float(t[7:])
                if t.find(contam_pop[0]+'_AF=') == 0:
                    this_contam_pop0_AF = float(t[7:])
                if t.find(contam_pop[1]+'_AF=') == 0:
                    this_contam_pop1_AF = float(t[7:])

            # exclude SNPs below the threshold of pop frequency 
            if this_pedigree_pop_AF < min_AF:     
                include_SNP = 0
        
        # perform additional filtering
        if include_SNP == 1:
            if (random.random() > downsample_rate_numSNPs):
                include_SNP = 0
        
                
        # if necessary, replace contaminating allele frequency with frequency calculated from specific contaminating individuals
        if include_SNP == 1:
            #print 'contam_ind_IDs[0] =', contam_ind_IDs[0]
            if len(contam_ind_IDs[0]) > 0:
                # extract allele counts for contaminating individuals                     
                REF_allele_count = 0
                ALT_allele_count = 0
                for ind_index in contam_ind_IDs[0]:
                    ind_alleles = elem[ind_index].split('|')
                    #print 'ind_alleles =', ind_alleles
                    for a in [0, 1]:
                        if ind_alleles[a] == '0':
                            REF_allele_count += 1
                        elif ind_alleles[a] == '1':
                            ALT_allele_count += 1
                        else: 
                            print 'error parsing alleles from contaminating individual genome!'
                # calculate new value for this_contam_pop0_AF from allele counts of contaminating individual genomes
                #print 'changing this_contam_pop0_AF from', this_contam_pop0_AF, 'to', float(ALT_allele_count)/float(ALT_allele_count+REF_allele_count)
                this_contam_pop0_AF = float(ALT_allele_count)/float(ALT_allele_count+REF_allele_count)
            
            #print 'contam_ind_IDs[1] =', contam_ind_IDs[1]
            if len(contam_ind_IDs[1]) > 0:
                # extract allele counts for contaminating individuals                     
                REF_allele_count = 0
                ALT_allele_count = 0
                for ind_index in contam_ind_IDs[1]:
                    ind_alleles = elem[ind_index].split('|')
                    #print 'ind_alleles =', ind_alleles
                    for a in [0, 1]:
                        if ind_alleles[a] == '0':
                            REF_allele_count += 1
                        elif ind_alleles[a] == '1':
                            ALT_allele_count += 1
                        else: 
                            print 'error parsing alleles from contaminating individual genome!'
                # calculate new value for this_contam_pop1_AF from allele counts of contaminating individual genomes
                #print 'changing this_contam_pop1_AF from', this_contam_pop1_AF, 'to', float(ALT_allele_count)/float(ALT_allele_count+REF_allele_count)
                this_contam_pop1_AF = float(ALT_allele_count)/float(ALT_allele_count+REF_allele_count)

        if include_SNP == 1:
            # downsampling of variation (heterozygosity)
            # determine if freq of this variant should be changed to 0 in simulations
            dsSNP = 0
            if (random.random() > downsample_rate_variation):
                # do change freq of this SNP to 0 in the EUR pedigree pop
                dsSNP = 1      
            else:
                # do not change freq of this SNP (to 0) in the pedigree pop
                dsSNP = 0

            # update the chromosome position lookup and reverse-lookup dictionaries
            pos_lookup[pos_vector_length] = pos_on_chromosome
            if (pos_on_chromosome in pos_rev_lookup) == False:
                pos_rev_lookup[pos_on_chromosome] = pos_vector_length
            else:
                # this variant position is listed more than once in the VCF for some reason
                print '        ERROR! chr', str(chr), 'pos', str(pos_on_chromosome+1), '(1-based) is already in pos_rev_lookup{}. ignoring all but first variant'
                #print '            pos_on_chromosome =', pos_on_chromosome
                #print '            len(pos_rev_lookup) =', len(pos_rev_lookup)
                #print '            pos_rev_lookup[pos_on_chromosome] =', pos_rev_lookup[pos_on_chromosome]
                #print '            pos_vector_length =', pos_vector_length
                # skip to next line in file
                continue
                
            # extract genotypes for each individual
            for ind in xrange(nind):
                #print 'ind =', ind
                #print 'pos_vector_length=', pos_vector_length
                #print 'pos_on_chromosome =', pos_on_chromosome
                #print 'colInDatafile[ind] =', colInDatafile[ind]
                #print 'elem[colInDatafile[ind]] =', elem[colInDatafile[ind]]
                #print 'elem[colInDatafile[ind]][0] =', elem[colInDatafile[ind]][0]
                #print 'elem[colInDatafile[ind]][2] =', elem[colInDatafile[ind]][2]
                #print '---'
                
                if dsSNP == 0:
                    # read this individuals' haps and save sequence data
                    haplo1[ind][pos_vector_length] = int(elem[colInDatafile[ind]][0])
                    haplo2[ind][pos_vector_length] = int(elem[colInDatafile[ind]][2])
                else:
                    haplo1[ind][pos_vector_length] = 0
                    haplo2[ind][pos_vector_length] = 0
                #print 'allele1, allele2 =', elem[colInDatafile[ind]][0], elem[colInDatafile[ind]][2]
                    
            # increment length of position vector
            pos_vector_length += 1
            
            # extract pop allele frequencies at this SNP
            INFO = elem[7].split(';')
            
            # save actual alternate allele frequency at this SNP for pedigree superpopulation            
            if dsSNP == 0:
                pedigree_pop_AF.append(this_pedigree_pop_AF)
            else:
                pedigree_pop_AF.append(float(0))             
            # save alternate allele frequency at this SNP for contaminating superpopulation(s)
            contam_pop_AF[0].append(this_contam_pop0_AF)
            contam_pop_AF[1].append(this_contam_pop1_AF)
                                        
        else:
            # ignore filtered variant
            continue   
    
    # re-assign the downsampled number of variants to numvars 
    numvars = pos_vector_length
    if numvars != len(pos_rev_lookup):
        print '        ERROR! pos_vector_length =', int(pos_vector_length), 'and len(pos_rev_lookup) =', len(pos_rev_lookup), 'at end of update_individuals_from1000g()\n'  
        
        
        
    # truncate haplo1 and haplo2 to downsampled size
    for n in range(0, nind):
        haplo1[n] = haplo1[n][0:pos_vector_length]
        haplo2[n] = haplo2[n][0:pos_vector_length]
        
    [listindiv[ind].update(haplo1[ind], haplo2[ind], sizeChrom)  for ind in xrange(nind)]
    
    
    return pos_rev_lookup, numvars, contam_pop_AF, pedigree_pop_AF



def repro(father, mother, theta, chromosome_map, chromsize, kidID, recomb_rate_regions, recomb_probs):
    """
    Outcome of reproduction
    """
    kid = Individual(kidID)
    kid.parent1 = father.ID
    kid.parent2 = mother.ID
    print '        performing meiosis in parent1 ('+str(kid.parent1)+') ...'
    gam1_father, gam2_father = father.meiosis(theta, chromosome_map, recomb_rate_regions, recomb_probs)
    print '        performing meiosis in parent2 ('+str(kid.parent2)+') ...'
    gam1_mother, gam2_mother = mother.meiosis(theta, chromosome_map, recomb_rate_regions, recomb_probs)
    
    if np.random.random() < 0.5:
        hap1_kid = gam1_father
    else:
        hap1_kid = gam2_father
    if np.random.random() < 0.5:
        hap2_kid = gam1_mother
    else:
        hap2_kid = gam2_mother    
    kid.update(hap1_kid, hap2_kid, chromsize)

    return kid




def pairwise_diff(indiv1, indiv2, mean_coverage1, mean_coverage2, contam_rate1, contam_rate2, error_rate1, error_rate2, contam_pop_AF, pedigree_pop_AF, chromosome_map):  
    total_SNPs = len(indiv1.hap1)
    
    # create list of all sites (containing /1 for pairwise difference) for block jackknife on data
    PWD_data_out = [['', -1, -1]]*total_SNPs
    
    # simulate drawing a single read from the true genotypes of each individual
    # FIRST, simulate a number of reads covering the SNP site from a poisson distribution with mean equal to mean_coverage
    SNPcov_indiv1 = np.random.poisson(mean_coverage1, total_SNPs)
    SNPcov_indiv2 = np.random.poisson(mean_coverage2, total_SNPs)
    
    # create empty vectors for read observations
    indiv1_obs = np.empty(total_SNPs, dtype=np.int8)
    indiv2_obs = np.empty(total_SNPs, dtype=np.int8)
        
    # SECOND, for individual 1, simulate from a binomial distribution a single read for those sites that have some coverage
    #print 'simulating cov of indiv1'
    for g in range(0, total_SNPs):
        obsread = -1
        if SNPcov_indiv1[g] != 0:
            # simulate the alleles observed in g reads at this site, assign the single observed allele to indiv1
            obsread = simObservedRead(SNPcov_indiv1[g], contam_rate1, contam_pop_AF[0][g], error_rate1, indiv1.hap1[g], indiv1.hap2[g])
            indiv1_obs[g] = obsread

        else:
            # mark this site for removal from the difference comparison
            indiv1_obs[g] = -1
            

    # THIRD, for individual 2, simulate from a binomial distribution a single read for those sites that have some coverage
    #print 'simulating cov of indiv2'
    for g in range(0, total_SNPs):
        obsread = -1
        if SNPcov_indiv2[g] != 0:
            # simulate the alleles observed in g reads at this site, assign the single observed allele to indiv1
            obsread = simObservedRead(SNPcov_indiv2[g], contam_rate2, contam_pop_AF[1][g], error_rate2, indiv2.hap1[g], indiv2.hap2[g])
            #print 'indiv2 obsread =', obsread
            indiv2_obs[g] = obsread        
        else:
            # mark this site for removal from the difference comparison
            indiv2_obs[g] = -1
    
    # identify overlapping sites
    overlap_sites = np.where(np.logical_and(indiv1_obs != -1, indiv2_obs != -1) == True)[0]
    
    # count overlapping sites from both lists of observed alleles
    num_SNPs_covered = len(overlap_sites)
    
    # calculate the number of pairwise diffs only examining overlapping sites
    pairwise_diff_positions = np.where(np.logical_and(np.logical_and(indiv1_obs != -1, indiv2_obs != -1), indiv1_obs != indiv2_obs) == True)[0]
    pairwise_diffs = len(pairwise_diff_positions)
    fraction = float(pairwise_diffs) / num_SNPs_covered
    
    # save PWD position data for export
    for T in range(0, total_SNPs):
        if T in pairwise_diff_positions:
            PWD_data_out[T] = [chromosome_map[T], 1]
        else:
            #print 'total_SNPs =', total_SNPs
            #print 'T =', T
            #print 'len(PWD_data_out) =', len(PWD_data_out)
            #print 'len(chromosome_map) =', len(chromosome_map)
            #print 'chromosome_map[T] =', chromosome_map[T]
            PWD_data_out[T] = [chromosome_map[T], 0]        
    
        
    # calculate theoretical expected values (including contam, assume contam_rate1=contam_rate2) for pairwise differences considering only SNPs covered by both individuals
    sum_unrelated   = 0.0   # unrelated relationship
    sum_sibs        = 0.0   # sibling relationship
    sum_gpgc        = 0.0   # grandparent/grandchild or uncle/nephew relationship
    sum_twins       = 0.0   # self relationship
    sum_fc          = 0.0   # first cousin relationship
    for X in overlap_sites:
        c1 = contam_rate1
        c2 = contam_rate2
        pN = float(pedigree_pop_AF[X])       
        pC = float(contam_pop_AF[0][X] + contam_pop_AF[1][X]) / 2.0     # use mean value until theoretical expectations with 2 contaminating allele freqs are implemented
        q = (error_rate1 + error_rate2) / 2.0                           # use mean value until theoretical expectations with 2 error rates are implemented
                                        
        term2 = c1*(1-c2)*  (pN + pC - (2*pN*pC) + (2*q) - ((8.0/3)*pN*q)    - ((8.0/3)*pC*q)    + ((16.0/3)*pN*pC*q)    - ((4.0/3)*q*q) + ((16.0/9)*pN*q*q) + ((16.0/9)*pC*q*q) - ((32.0/9)*pN*pC*q*q))
        term3 = c2*(1-c1)*  (pN + pC - (2*pN*pC) + (2*q) - ((8.0/3)*pN*q)    - ((8.0/3)*pC*q)    + ((16.0/3)*pN*pC*q)    - ((4.0/3)*q*q) + ((16.0/9)*pN*q*q) + ((16.0/9)*pC*q*q) - ((32.0/9)*pN*pC*q*q))
        term4 = c1*c2*      (pC + pC - (2*pC*pC) + (2*q) - ((8.0/3)*pC*q)    - ((8.0/3)*pC*q)    + ((16.0/3)*pC*pC*q)    - ((4.0/3)*q*q) + ((16.0/9)*pC*q*q) + ((16.0/9)*pC*q*q) - ((32.0/9)*pC*pC*q*q))
  
        sum_unrelated = sum_unrelated   + term2 + term3 + term4 +   (  (1-c1)*(1-c2)* ((2*pN)           - (2*pN*pN)         + (2*q) - ((16.0/3)*pN*q)   + ((16.0/3)*pN*pN*q)    - ((4.0/3)*q*q) + ((32.0/9)*pN*q*q) - ((32.0/9)*pN*pN*q*q)) )
        sum_sibs = sum_sibs             + term2 + term3 + term4 +   (  (1-c1)*(1-c2)* (((3.0/2)*pN)     - ((3.0/2)*pN*pN)   + (2*q) - (4*pN*q)          + (4*pN*pN*q)           - ((4.0/3)*q*q) + ((8.0/3)*pN*q*q)  - ((8.0/3)*pN*pN*q*q))  )                                                                                                                                                               
        sum_gpgc = sum_gpgc             + term2 + term3 + term4 +   (  (1-c1)*(1-c2)* (((7.0/4)*pN)     - ((7.0/4)*pN*pN)   + (2*q) - ((14.0/3)*pN*q)   + ((14.0/3)*pN*pN*q)    - ((4.0/3)*q*q) + ((28.0/9)*pN*q*q) - ((28.0/9)*pN*pN*q*q)) )                                                                                                           
        sum_fc = sum_fc                 + term2 + term3 + term4 +   (  (1-c1)*(1-c2)* (((15.0/8)*pN)    - ((15.0/8)*pN*pN)  + (2*q) - (5*pN*q)          + (5*pN*pN*q)           - ((4.0/3)*q*q) + ((10.0/3)*pN*q*q) - ((10.0/3)*pN*pN*q*q)) )
        sum_twins = sum_twins           + term2 + term3 + term4 +   (  (1-c1)*(1-c2)* (pN               - (pN*pN)           + (2*q) - ((8.0/3)*pN*q)    + ((8.0/3)*pN*pN*q)     - ((4.0/3)*q*q) + ((16.0/9)*pN*q*q) - ((16.0/9)*pN*pN*q*q)) )  

    return [fraction, pairwise_diffs, num_SNPs_covered, total_SNPs, sum_unrelated, sum_sibs, sum_gpgc, sum_twins, sum_fc, PWD_data_out]



def simObservedRead(num_observed_reads, contam_rate, contam_pop_AF, seq_errorrate, hap1, hap2):
    # simulate the alleles observed in g reads at this site
    alleles_observed = []
    for b in range(0, num_observed_reads):
        
        # simulate contamination
        if np.random.random() < contam_rate:
            # read is from the contaminating population
            if np.random.random() < contam_pop_AF:                          
                # contaminant read is alt allele
                chosen = 1
            else:
                # contaminant read is ref allele 
                chosen = 0
        else:
            # read is from the true source individual
            chosen = random.choice([hap1, hap2])
        
        # simulate sequencing error
        if np.random.random() < seq_errorrate:
            # base error at the SNP, observed base is not the chosen (alt|ref) allele
            alleles_observed.append(random.choice(seqerror_choices[chosen]))
        else:
            # no base error at the SNP, observed base is the chosen (alt|ref) allele
            alleles_observed.append(chosen)
    
    #print alleles_observed
    return random.choice(alleles_observed)            
    

        
def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

        
    
def pileup_PWD(indiv1, indiv2, contam_rate1, contam_rate2, error_rate1, error_rate2, contam_pop_AF, pedigree_pop_AF, pos_rev_lookup, pileup_positions, verbose):  
    # initialize counter variables
    overlap_sites = 0
    pairwise_diffs = 0
    
    # initialize sums for theoretical expected values
    sum_unrelated   = 0.0   # unrelated relationship
    sum_sibs        = 0.0   # sibling relationship
    sum_gpgc        = 0.0   # grandparent/grandchild or uncle/nephew relationship
    sum_twins       = 0.0   # self relationship
    sum_fc          = 0.0   # first cousin relationship
    
    # create list of all sites (containing /1 for pairwise difference) for block jackknife on data
    pileup_positions_keys = pileup_positions.keys()
    PWD_data_out = [['', -1, -1]]*len(pileup_positions_keys)
    
    for position in pileup_positions_keys:    
        # if this is a known variant site 
        if int(position in pos_rev_lookup) == 1:     

            # define chromosome_map vector position of this chromosomal position
            map_position = pos_rev_lookup[position]
                       
            # assign actual haplotypes
            actual_indiv1_hap1 = indiv1.hap1[map_position]
            actual_indiv1_hap2 = indiv1.hap2[map_position]
            actual_indiv2_hap1 = indiv2.hap1[map_position]
            actual_indiv2_hap2 = indiv2.hap2[map_position]

            # assign allele frequencies for theoretical expectations calcs
            pN = float(pedigree_pop_AF[map_position])       
            #pC = float(contam_pop_AF[map_position])  
            if contam_rate1+contam_rate2 != 0:
                pC = (contam_rate1*float(contam_pop_AF[0][map_position])/(contam_rate1+contam_rate2)) + (contam_rate2*float(contam_pop_AF[1][map_position])/(contam_rate1+contam_rate2))
                # this is just a weighted mean of the two contamination pops' allele frequencies because theory isn't worked out
            else:
                pC = float(contam_pop_AF[random.choice([0,1])][map_position])

        else:
            # if this is not a known variant site, assign only reference allele 0 to both individuals
            actual_indiv1_hap1 = 0
            actual_indiv1_hap2 = 0
            actual_indiv2_hap1 = 0
            actual_indiv2_hap2 = 0
            
            # assign allele frequencies for theoretical expectations calcs
            pN = 0.0     
            pC = 0.0

        # count this site as an overlapping base
        overlap_sites += 1
        
    
        # determine which error rates to use for this position
        if error_rate1 != 'pileup':
            error_rateA = error_rate1
            error_rateB = error_rate2
        else:
            error_rateA = (10.0 ** (random.choice(pileup_positions[position][0])*(-1.0/10)))
            error_rateB = (10.0 ** (random.choice(pileup_positions[position][1])*(-1.0/10)))

        # randomly choose a single base quality from each sample, and simulate their observed reads         
        chosen_base1 = simObservedRead(1, contam_rate1, pC, error_rateA, actual_indiv1_hap1, actual_indiv1_hap2)        
        chosen_base2 = simObservedRead(1, contam_rate2, pC, error_rateB, actual_indiv2_hap1, actual_indiv2_hap2)   
        
        if chosen_base1 != chosen_base2:
            pairwise_diffs += 1
            PWD_data_out[overlap_sites-1] = [position, 1]
        else:
            PWD_data_out[overlap_sites-1] = [position, 0]
        
        if verbose == 1:
            print '---'
            print 'position:', position
            print 'pileup:', pileup_positions[position][0], pileup_positions[position][1]
            print 'actual indiv1 haplotype:', actual_indiv1_hap1,'/', actual_indiv1_hap2 
            print 'actual indiv2 haplotype:', actual_indiv2_hap1,'/', actual_indiv2_hap2
            print 'sampled bases:', chosen_base1, chosen_base2 
           
        
        # use mean of the two assigned seq error rates for theoretical calculations
        q = (error_rateA + error_rateB) / 2.0
        
        # calculate this site's contribution to theoretical expected values (including contam, assume contam_rate1=contam_rate2) for pairwise differences considering only SNPs covered by both individuals
        c1 = contam_rate1
        c2 = contam_rate2

        term2 = c1*(1-c2)*  (pN + pC - (2*pN*pC) + (2*q) - ((8.0/3)*pN*q)    - ((8.0/3)*pC*q)    + ((16.0/3)*pN*pC*q)    - ((4.0/3)*q*q) + ((16.0/9)*pN*q*q) + ((16.0/9)*pC*q*q) - ((32.0/9)*pN*pC*q*q))
        term3 = c2*(1-c1)*  (pN + pC - (2*pN*pC) + (2*q) - ((8.0/3)*pN*q)    - ((8.0/3)*pC*q)    + ((16.0/3)*pN*pC*q)    - ((4.0/3)*q*q) + ((16.0/9)*pN*q*q) + ((16.0/9)*pC*q*q) - ((32.0/9)*pN*pC*q*q))
        term4 = c1*c2*      (pC + pC - (2*pC*pC) + (2*q) - ((8.0/3)*pC*q)    - ((8.0/3)*pC*q)    + ((16.0/3)*pC*pC*q)    - ((4.0/3)*q*q) + ((16.0/9)*pC*q*q) + ((16.0/9)*pC*q*q) - ((32.0/9)*pC*pC*q*q))

        sum_unrelated = sum_unrelated   + term2 + term3 + term4 +   (  (1-c1)*(1-c2)* ((2*pN)           - (2*pN*pN)         + (2*q) - ((16.0/3)*pN*q)   + ((16.0/3)*pN*pN*q)    - ((4.0/3)*q*q) + ((32.0/9)*pN*q*q) - ((32.0/9)*pN*pN*q*q)) )
        sum_sibs = sum_sibs             + term2 + term3 + term4 +   (  (1-c1)*(1-c2)* (((3.0/2)*pN)     - ((3.0/2)*pN*pN)   + (2*q) - (4*pN*q)          + (4*pN*pN*q)           - ((4.0/3)*q*q) + ((8.0/3)*pN*q*q)  - ((8.0/3)*pN*pN*q*q))  )                                                                                                                                                               
        sum_gpgc = sum_gpgc             + term2 + term3 + term4 +   (  (1-c1)*(1-c2)* (((7.0/4)*pN)     - ((7.0/4)*pN*pN)   + (2*q) - ((14.0/3)*pN*q)   + ((14.0/3)*pN*pN*q)    - ((4.0/3)*q*q) + ((28.0/9)*pN*q*q) - ((28.0/9)*pN*pN*q*q)) )                                                                                                           
        sum_fc = sum_fc                 + term2 + term3 + term4 +   (  (1-c1)*(1-c2)* (((15.0/8)*pN)    - ((15.0/8)*pN*pN)  + (2*q) - (5*pN*q)          + (5*pN*pN*q)           - ((4.0/3)*q*q) + ((10.0/3)*pN*q*q) - ((10.0/3)*pN*pN*q*q)) )
        sum_twins = sum_twins           + term2 + term3 + term4 +   (  (1-c1)*(1-c2)* (pN               - (pN*pN)           + (2*q) - ((8.0/3)*pN*q)    + ((8.0/3)*pN*pN*q)     - ((4.0/3)*q*q) + ((16.0/9)*pN*q*q) - ((16.0/9)*pN*pN*q*q)) )  

    # to prevent crashing when overlap_sites = 0
    if overlap_sites == 0:
        fraction = float(0)   
    else:
        fraction = float(pairwise_diffs) / overlap_sites    
    
    return [fraction, pairwise_diffs, overlap_sites, len(indiv1.hap1), sum_unrelated, sum_sibs, sum_gpgc, sum_twins, sum_fc, PWD_data_out]
    
    
    
def getRecombinationRates(recomb_rates_file, verbose):
    print '    loading region-specific recombination rates from:', recomb_rates_file
    openfile = open(recomb_rates_file, "rU")
    open_recomb_rates_file = csv.reader(openfile, dialect=csv.excel_tab, lineterminator='\n', quoting=csv.QUOTE_NONE)
    open_recomb_rates_file.next()
    recomb_rate_regions = []
    recomb_probs = [] # probability of an odd number of crossover events per base per generation
    reg_start = -1
    reg_stop = -1 
    for row in open_recomb_rates_file:
        if reg_start == -1:
            reg_start = int(row[1])
        elif reg_stop == -1:
            reg_stop = int(row[1])-1
            recomb_rate_regions.append((reg_start, reg_stop))
            recomb_probs.append(float(row[2])/100/1000000)
            reg_start = reg_stop+1
            reg_stop = -1
    openfile.close()
    
    # check the recombination probabilities for errors
    if verbose == 1:
        len_sum = 0
        prob_sum = 0
        for X in range(len(recomb_rate_regions)):
            len_sum += recomb_rate_regions[X][1]-recomb_rate_regions[X][0]+1
            prob_sum += recomb_probs[X]*(recomb_rate_regions[X][1]-recomb_rate_regions[X][0]+1)
        print '        sum of per-base recombination probs over all bases:', prob_sum
        print '        sum of region lengths where recomb rate != 0:', len_sum

    return recomb_rate_regions, recomb_probs                        
 
    
    
def getQualitiesFromPileup(pileup_positions, qualities_pileup_filelist, self_comparison, min_qual, verbose, chr):
    # read in pileup file, creating a dictionary of read depths and qualities for each position on this chromosome 
    print 'loading base qualities from pileup file '+qualities_pileup_filelist+' ...' 
    # open pileup file
    input_file = gzip.open(qualities_pileup_filelist)
    pileup_data = csv.reader(input_file, dialect=csv.excel_tab, lineterminator='\n', quoting=csv.QUOTE_NONE)

        # expected pileup file format
        # --------------------------------------------------------------------------------------------------------------------------------------------------------
        # 0                                                 1                   2           3           4           5           6           7           8
        # chromosome            	                        1-based position    reference   bam1_depth  bam1_bases  bam1_quals  bam2_depth  bam2_bases  bam2_quals
        # supercont1.4920                               	297	                C	        1	        ,	        e           1           ,           e           
        # --------------------------------------------------------------------------------------------------------------------------------------------------------

    # initialize counter variables and pileup_positions dictionary
    lines_read = 0
    #pileup_positions = {}   # pileup_positions[pos] = [[basequal1, basequal2, ...], [basequal1, basequal2, ...] ]

    # parse each line of the provided pileup file
    for row in pileup_data:
        if lines_read % 1000000 == 0 and lines_read != 0:
            print '    pileup lines read: '+str(lines_read/1000000)+'M'

        # if pileup file chromosome number is larger than the chromosome number being simulated
        if RepresentsInt(row[0]) == False:
            # stop reading through the file
            break            
        #elif int(row[0]) > chr:
        #    # stop reading through the file
        #    break
        elif chr+1 == int(row[0]):

            # extract read depth and associated base qualities for each sample
            pileup_quals = [[], []]
            pileup_bases = [[], []]
            
            for K in [0, 1]:
                # read the base qual field
                if K == 0:
                    pileup_line = row[4]
                    quals_line = row[5]
                elif K == 1:
                    pileup_line = row[7]
                    quals_line = row[8]

                line_pos = 0
                qual_pos = 0
                while line_pos < len(pileup_line):
                    # if base matches reference
                    if pileup_line[line_pos] == '.' or pileup_line[line_pos] == ',':
                        pileup_bases[K].append(row[2])
                        pileup_quals[K].append(int(ord(quals_line[qual_pos])-33))
                        #print int(ord(quals_line[qual_pos])-33)
                        line_pos += 1
                        qual_pos += 1

                    # if base is a mismatch to reference   
                    elif int(pileup_line[line_pos] in mismatch_codes) == 1:
                        pileup_bases[K].append(mismatch_codes[pileup_line[line_pos]])
                        pileup_quals[K].append(int(ord(quals_line[qual_pos])-33))
                        #print int(ord(quals_line[qual_pos])-33)

                        # if this is a deleted position and ignore_dels option is ON
                        #if ignore_dels == 1 and pileup_line[line_pos] == '*':
                        if pileup_line[line_pos] == '*':
                            pass
                        else:
                            pileup_bases[K].append(mismatch_codes[pileup_line[line_pos]])
                            pileup_quals[K].append(int(ord(quals_line[qual_pos])-33))
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

                        # find the number of characters encoding the indel size integer
                        looker = line_pos+1
                        while int(pileup_line[looker].isdigit()) == 1:
                            looker += 1
                        indel_encode_digits = looker - line_pos - 1

                        # count this indel as an uncalled base
                        indel_size = int(pileup_line[(line_pos+1):(line_pos+1+indel_encode_digits)])                
                        line_pos += (1 + indel_encode_digits + indel_size)
                        qual_pos += 0

                    else:
                        print 'ERROR!'

            # filter bases with quality score <= min_qual
            base_choices1 = []
            base_choices2 = []
            qual_choices1 = []
            qual_choices2 = []
            
            # if self_comparison mode is OFF
            if self_comparison == 0:
                for t in range(0, len(pileup_quals[0])):
                    if pileup_quals[0][t] >= min_qual:
                        qual_choices1.append(pileup_quals[0][t])
                for t in range(0, len(pileup_quals[1])):
                    if pileup_quals[1][t] >= min_qual:
                        qual_choices2.append(pileup_quals[1][t])

                for t in range(0, len(pileup_bases[0])):
                    exclude_base = 0
                    if pileup_quals[0][t] >= min_qual:
                        if exclude_base == 0:
                            base_choices1.append(pileup_bases[0][t])
                            qual_choices1.append(pileup_quals[0][t])
                for t in range(0, len(pileup_bases[1])):
                    exclude_base = 0
                    if pileup_quals[1][t] >= min_qual:
                        if exclude_base == 0:
                            base_choices2.append(pileup_bases[1][t])
                            qual_choices2.append(pileup_quals[1][t])                    
            
            
            # if self_comparison mode is ON
            else:
                # randomly alternate assignment of each base quality to qual_choices1 and qual_choices2
                qual_indices = range(0, len(pileup_quals[0]))
                random.shuffle(qual_indices)
                assigned = 0
                while len(qual_indices) > 0:
                    current_index = qual_indices.pop()
                    current_qual = pileup_quals[0][current_index]
                    current_base = pileup_bases[0][current_index]    
                    exclude_base = 0
                    if current_qual >= min_qual:
                        if exclude_base == 0:
                            if assigned % 2 == 0:
                                qual_choices1.append(current_qual) 
                            else:
                                qual_choices2.append(current_qual) 
                            assigned += 1
            
            # if after filtering, there is at least 1 observed base to simulate in both samples
            if len(qual_choices1) > 0 and len(qual_choices2) > 0:
                # this position passes filters and should be simulated 
                #pileup_positions[int(row[1])][0] = qual_choices1
                #pileup_positions[int(row[1])][1] = qual_choices2
                # pysam/BAM uses a 0-based coordinate system, so must subtract 1 from each coordinate
                pileup_positions[int(row[1])-1] = [qual_choices1, qual_choices2]
            elif verbose == 1:
                print 'this pileup position did not pass filters'
                print row
        lines_read += 1
    print '    done. hashed '+str(len(pileup_positions))+' genomic positions to pileup_positions.'  
    
    return pileup_positions
    