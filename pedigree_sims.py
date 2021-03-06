# define constants
recomb_rate = 0.000000013   # intra-locus uniform recombination rate [Prufer et al 2014, SI 16a]


# import packages
import sys
import os.path
import datetime
import pysam
import random
#import code #code.interact(local=locals())
import gzip
import numpy as np
import math
import pyximport; pyximport.install()
import grups_module as mfp
import argparse


# define functions
def writePWDtoFiles(relationship_index):
    outfile[param_set][rep].write(str(float(pw_sum))+'\t')
    outfile_numvars[param_set][rep].write(str(numvars_actual)+'\t')
    #outfile_pwdiffs[param_set][rep].write(str(pwdiffsresults[4])+'\t')
    #outfile_pwdiffs_sibs[param_set][rep].write(str(pwdiffsresults[5])+'\t')     #sum_sibs
    #outfile_pwdiffs_gpgc[param_set][rep].write(str(pwdiffsresults[6])+'\t')     #sum_gpargchild
    #outfile_pwdiffs_twins[param_set][rep].write(str(pwdiffsresults[7])+'\t')    #sum_IDtwins
    #outfile_pwdiffs_fcous[param_set][rep].write(str(pwdiffsresults[8])+'\t')    #sum_firstcousin
    outfile_SNPscovered[param_set][rep].write(str(pwdiffsresults[2])+'\t') 
    #if writeAllData == 1:
    #    all_PWD_data = sorted(pwdiffsresults[9], key=lambda element: element[0])
    #    for position in range(0, len(all_PWD_data)):
    #        #print 'position =', position
    #        #print 'all_PWD_data[position] =', all_PWD_data[position]
    #        outfile_PWDdata[param_set][relationship_index].write(str(chr+1)+'\t'+ \
    #                                                             str(all_PWD_data[position][0])+'\t'+ \
    #                                                             str(all_PWD_data[position][1])+'\n') 

def openOutputFiles():
    outfile[param_set].append(                          open(out_dir+scriptname+'.'+output_label[param_set]+'.paramRep'+str(rep)+'.'+timenow+'.out',                'a', 0))
    outfile_numvars[param_set].append(                  open(out_dir+scriptname+'.'+output_label[param_set]+'.paramRep'+str(rep)+'.'+timenow+'.numvars',            'a', 0))
    #outfile_pwdiffs[param_set].append(                  open(out_dir+scriptname+'.'+output_label[param_set]+'.paramRep'+str(rep)+'.'+timenow+'.pwexp',              'a', 0))    # expected pairwise difference for this chromosome for unrelated relationships 
    #outfile_pwdiffs_sibs[param_set].append(             open(out_dir+scriptname+'.'+output_label[param_set]+'.paramRep'+str(rep)+'.'+timenow+'.pwexp_sibs',         'a', 0))    # expected pairwise difference for this chromosome for siblings + parent-child relationships 
    #outfile_pwdiffs_gpgc[param_set].append(             open(out_dir+scriptname+'.'+output_label[param_set]+'.paramRep'+str(rep)+'.'+timenow+'.pwexp_gpgc',         'a', 0))    # expected pairwise difference for this chromosome for gparent/gchild + uncle/newphew + half-sibling relationships
    #outfile_pwdiffs_twins[param_set].append(            open(out_dir+scriptname+'.'+output_label[param_set]+'.paramRep'+str(rep)+'.'+timenow+'.pwexp_twins',        'a', 0))    # expected pairwise difference for this chromosome for identical twin relationships
    #outfile_pwdiffs_fcous[param_set].append(            open(out_dir+scriptname+'.'+output_label[param_set]+'.paramRep'+str(rep)+'.'+timenow+'.pwexp_fcous',        'a', 0))    # expected pairwise difference for this chromosome for first cousin relationships
    outfile_SNPscovered[param_set].append(              open(out_dir+scriptname+'.'+output_label[param_set]+'.paramRep'+str(rep)+'.'+timenow+'.numSNPscov',         'a', 0))
    #outfile_SNPfreqs[param_set].append(                open(out_dir+scriptname+'.'+output_label[param_set]+'.paramRep'+str(rep)+'.'+timenow+'.SNPfreqs',           'a', 0))    # row of freq(EUR_AF) and freq(pedigree population) for each included SNP
    outfile_labels[param_set].append(                   open(out_dir+scriptname+'.'+output_label[param_set]+'.paramRep'+str(rep)+'.'+timenow+'.labels',             'a', 0))    # column labels
    #if writeAllData == 1:
    #    outfile_PWDdata[param_set].append([])
    #    for rel in labels:
    #        outfile_PWDdata[param_set][rep].append(     open(out_dir+scriptname+'.'+output_label[param_set]+'.paramRep'+str(rep)+'.'+rel+'.'+timenow+'.PWDdata',    'a', 0))    # 0/1 simulated for pairwise difference at each site

def writeOutputFileEndlines():
    outfile[param_set][rep].write('\n')
    outfile_numvars[param_set][rep].write('\n')
    #outfile_pwdiffs[param_set][rep].write('\n')
    #outfile_pwdiffs_sibs[param_set][rep].write('\n')
    #outfile_pwdiffs_gpgc[param_set][rep].write('\n')
    #outfile_pwdiffs_twins[param_set][rep].write('\n')
    #outfile_pwdiffs_fcous[param_set][rep].write('\n')
    outfile_SNPscovered[param_set][rep].write('\n')
    #outfile_SNPfreqs[param_set][rep].write('\n')
    outfile_labels[param_set][rep].write('\n')
    #if writeAllData == 1:
    #    for R in range(len(outfile_PWDdata[param_set][rep])):
    #        outfile_PWDdata[param_set][rep][R].write('\n')

def closeOutputFiles():
    outfile[param_set][rep].close()
    outfile_numvars[param_set][rep].close()
    #outfile_pwdiffs[param_set][rep].close()
    #outfile_pwdiffs_sibs[param_set][rep].close()
    #outfile_pwdiffs_gpgc[param_set][rep].close()
    #outfile_pwdiffs_twins[param_set][rep].close()
    #outfile_pwdiffs_fcous[param_set][rep].close()
    outfile_SNPscovered[param_set][rep].close()
    #outfile_SNPfreqs[param_set][rep].close()
    outfile_labels[param_set][rep].close()
    #if writeAllData == 1:
    #    for rel in range(len(outfile_PWDdata[param_set][rep])):
    #        outfile_PWDdata[param_set][rep][rel].close()

def resetFileLists():
    outfile[param_set]                  = []
    outfile_numvars[param_set]          = []
    #outfile_pwdiffs[param_set]          = []
    #outfile_pwdiffs_sibs[param_set]     = []
    #outfile_pwdiffs_gpgc[param_set]     = []
    #outfile_pwdiffs_twins[param_set]    = []
    #outfile_pwdiffs_fcous[param_set]    = []
    outfile_SNPscovered[param_set]      = []
    #outfile_SNPfreqs[param_set]        = []
    outfile_labels[param_set]           = []
    #if writeAllData == 1:
    #    outfile_PWDdata[param_set]      = []

def parsePedigreeFile(ped_defintion):
    import csv
    
    input_file = open(ped_defintion, "rU")
    read_data = csv.reader(input_file, dialect=csv.excel_tab, lineterminator='\n', quoting=csv.QUOTE_NONE)
    relateds = []
    parsing_individuals = 0 
    parsing_relationships = 0
    parsing_comparisons = 0
    for row in read_data:
        if row != '' and len(row) != 0:
            if row[0][0] != "#":
                if row == ["INDIVIDUALS"]:
                    parsing_individuals = 1
                    parsing_relationships = 0
                    parsing_comparisons = 0
                elif row == ["RELATIONSHIPS"]:
                    parsing_individuals = 0
                    parsing_relationships = 1
                    parsing_comparisons = 0
                    print 'defined INDIVIDUALS:'
                    for x in individuals:
                        print '    ', x
                elif row == ["COMPARISONS"]:
                    parsing_individuals = 0
                    parsing_relationships = 0
                    parsing_comparisons = 1
                    print 'defined RELATIONSHIPS:'
                    for x in range(0, len(relationships)):
                        print '    ', relationships[x][0], 'offspring of', relationships[x][1], 'and', relationships[x][2]
                elif parsing_individuals == 1:
                    individuals.append(row[0])
                elif parsing_relationships == 1:
                    temp = row[0].split('=')
                    offspring = temp[0]
                    temp = temp[1].strip(')')
                    parents = temp.split('(')[1].split(',')
                    relationships.append([offspring, parents[0], parents[1]])
                    relateds.append(offspring)
                elif parsing_comparisons == 1:
                    temp = row[0].split('=')
                    printlabel = temp[0]
                    temp = temp[1].split(')')
                    pair = temp[0].split('(')[1].split(',')
                    comparisons.append([printlabel, pair[0], pair[1]])
    print 'defined COMPARISONS:'
    for x in range(0, len(comparisons)):
        print '    <', comparisons[x][0], '> label for comparison of', comparisons[x][1], 'and', comparisons[x][2]
    input_file.close()




# parse command-line constants and options 
scriptname = os.path.split(sys.argv[0])[1] # filename of the script being run
print '\n'
shuffle = 1
quiet = 0
popAF_polymorph_filter = 'all'
qualities_pileup_filename = ''
recomb_dir = ''
verbose = 0
min_qual = 0
self_comparison = 0
seq_errorrate = ['pileup', 'pileup']
mean_coverage = 'pileup'
targets_file = ''
min_AF = 0.0
pedigree_pop='EUR'
contam_pop=['EUR', 'EUR']
downsample_rate = 1
ds_rate_numSNPs = 1
contam_numind = 0,0
#writeAllData=0
paramNumReps=[1]
chr_type = 'whole'
data_dir='./'
out_dir='./'
individuals = []
relationships = []
comparisons = []


parser = argparse.ArgumentParser()
parser.add_argument("--ds_rate", type=str, help="FLOAT, proportion of SNPs to keep at true frequency (i.e. NOT change to 0% frequency) [default: 1.0]")
parser.add_argument('--ds_rate_numSNPs', type=str, help="FLOAT, proportion of filtered SNP positions include in the analysis [default: 1.0]")            
parser.add_argument('--c_rate', type=str, help="contamination rates (or rate ranges) for each desired scenario. format: FLOAT,FLOAT or FLOAT,FLOAT/FLOAT,FLOAT ... or FLOAT-FLOAT,FLOAT-FLOAT/FLOAT-FLOAT,FLOAT-FLOAT ...") 
parser.add_argument('--q_rate', type=str, help="sequencing error rates  (or rate ranges) for each desired scenario. format: FLOAT,FLOAT or FLOAT,FLOAT/FLOAT,FLOAT ... or FLOAT-FLOAT,FLOAT-FLOAT/FLOAT-FLOAT,FLOAT-FLOAT ...")                                       
parser.add_argument('--mean_cov', type=str, help="mean sequencing depths (or depth ranges) for each desired scenario. format: FLOAT,FLOAT or FLOAT,FLOAT/FLOAT,FLOAT ... or FLOAT-FLOAT,FLOAT-FLOAT/FLOAT-FLOAT,FLOAT-FLOAT ...") 
parser.add_argument('--paramNumReps', type=str, help="number of replicates to perform when randomly drawing values from specified ranges for parameters c_rate, mean_cov, q_rate. format: INT or INT/INT/INT ... ") 
parser.add_argument('--label', type=str, help="data output labels. format: STR/STR/STR ...") 
parser.add_argument('--data_dir', type=str, help="path to 1000genomes input VCFs [default: working directory]") 
parser.add_argument('--recomb_dir', type=str, help="path to recombination genetic map data") 
parser.add_argument('--reps', type=str, help="INT, number of pedigree replicates to perform [1]") 
parser.add_argument('--min_qual', type=str, help="INT, minimum base quality (phred-scale) of a pileup base to be included in simulations [default: 0]")                                       
parser.add_argument('--min_AF', type=str, help="FLOAT, minimum allele frequency in pedigree superpopulation for SNP to be included in simulations [default: 0]")                                      
#parser.add_argument('--chr_type', type=str, help="")                                           
parser.add_argument('--pedigree_pop', type=str, help="EAS|AMR|AFR|EUR|SAS, superpopulation with which to perform pedigree simulations [default: EUR]")                                         
parser.add_argument('--contam_pop', type=str, help="superpopulation with which to contaminate pedigree simulations. format: STR,STR [default: EUR,EUR]")                                           
parser.add_argument('--include', type=str, help="all|poly|mono")                                        
parser.add_argument('--pileup', type=str, help="path to input pileup file to be used in simulations")                                            
parser.add_argument('--chr', type=str, help="which chromosome(s) to include in the simulations. The desired chromosomes are selected automatically when provided with input pileup and/or targets files. format: INT or INT-INT or INT,INT or INT-INT,INT ..., etc.")                                               
parser.add_argument('--verbose', type=str, help="turn on verbose output mode, e.g. for debugging [default: off]")                                            
parser.add_argument('--self_comparison', type=str, help="turn on self-comparison mode (only used when an input pileup file is defined) [default: off]")                                           
parser.add_argument('--out', type=str, help="path to directory for output [default: working directory]") 
parser.add_argument('--targets', type=str, help="path to file defining target genomic positions")
parser.add_argument('--contam_numind', type=str, help="numbers of random individual genomes with which to contaminate pedigree simulations. format: INT,INT. [default: contaminate with population allele frequencies]")   
#parser.add_argument('--writeAllData', type=str, help="")   
parser.add_argument('--ped', type=str, help="path to input pedigree definition file")   

args = parser.parse_args()
if args.ds_rate:
    downsample_rate = float(args.ds_rate)
    print 'ds_rate =', downsample_rate
if args.ds_rate_numSNPs:
    ds_rate_numSNPs = float(args.ds_rate_numSNPs)
    print 'ds_rate_numSNPs =', ds_rate_numSNPs
if args.c_rate:
    contam_rate = args.c_rate
    param_sets = contam_rate.split('/')
    contam_rate = []
    for k in range(len(param_sets)):
        contam_rate.append([[], []])
        ind_params = param_sets[k].split(',')
        #print 'ind_params =', ind_params
        for ind in range(len(ind_params)):
            if ind_params[ind].find('-') != -1:
                ind_param_range = ind_params[ind].split('-')
            else:
                ind_param_range = [ind_params[ind], ind_params[ind]]
            contam_rate[k][ind] = [float(ind_param_range[0]), float(ind_param_range[1])]
    print 'c_rate =', contam_rate
if args.q_rate:                                     
    seq_errorrate = args.q_rate
    param_sets = seq_errorrate.split('/')
    seq_errorrate = []
    for k in range(len(param_sets)):
        seq_errorrate.append([[], []])
        ind_params = param_sets[k].split(',')
        #print 'ind_params =', ind_params
        for ind in range(len(ind_params)):
            if ind_params[ind].find('-') != -1:
                ind_param_range = ind_params[ind].split('-')
            else:
                ind_param_range = [ind_params[ind], ind_params[ind]]
            seq_errorrate[k][ind] = [float(ind_param_range[0]), float(ind_param_range[1])]
    print 'q_rate =', seq_errorrate
if args.mean_cov:
    mean_coverage = args.mean_cov  
    param_sets = mean_coverage.split('/')
    mean_coverage = []
    for k in range(len(param_sets)):
        mean_coverage.append([[], []])
        ind_params = param_sets[k].split(',')
        #print 'ind_params =', ind_params
        for ind in range(len(ind_params)):
            if ind_params[ind].find('-') != -1:
                ind_param_range = ind_params[ind].split('-')
            else:
                ind_param_range = [ind_params[ind], ind_params[ind]]
            mean_coverage[k][ind] = [float(ind_param_range[0]), float(ind_param_range[1])]
            if mean_coverage[k][ind] == 0:
                print 'mean_coverage[k][ind] == 0. Exiting now ...'
                exit()
    print 'mean_cov =', mean_coverage
if args.paramNumReps:
    paramNumReps = args.paramNumReps
    param_sets = paramNumReps.split('/')
    paramNumReps = []
    for k in range(len(param_sets)):
        paramNumReps.append(int(param_sets[k]))
    print 'paramNumReps =', paramNumReps
if args.label:
    output_label = args.label
    output_label = output_label.split('/')        
    print 'label =', output_label
if args.data_dir:
    data_dir = args.data_dir
    print 'data_dir =', data_dir
if args.recomb_dir:
    recomb_dir = args.recomb_dir
    print 'recomb_dir =', recomb_dir
if args.reps:
    num_reps = int(args.reps)
    print 'reps =', num_reps
if args.min_qual:                                    
    min_qual = float(args.min_qual)
    print 'min_qual =', min_qual
if args.min_AF:
    min_AF = float(args.min_AF)
    print 'min_AF =', min_AF                       
#if args.chr_type:
#    chr_type = args.chr_type
#    print 'chr_type =', chr_type
if args.pedigree_pop:
    pedigree_pop = args.pedigree_pop
    print 'pedigree_pop =', pedigree_pop         
if args.contam_pop:                                        
    contam_pop = args.contam_pop
    contam_pop = contam_pop.split(',')
    print 'contam_pop =', contam_pop
if args.include:
    popAF_polymorph_filter = args.include
    print 'include =', popAF_polymorph_filter 
if args.pileup:                                           
    qualities_pileup_filename = args.pileup
    qualities_pileup_filelist = []
    print 'pileup =', qualities_pileup_filename
    for yyy in range(1, 23):
        qualities_pileup_filelist.append(qualities_pileup_filename[0:len(qualities_pileup_filename)-3]+'.chr'+str(yyy)+'.gz')
        if os.path.exists(qualities_pileup_filelist[-1]) == False:
            print 'error. file '+qualities_pileup_filelist[-1]+' not found. exiting.'
            sys.exit()
        else:
            print '    found file '+qualities_pileup_filelist[-1]
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
    chrlist = [x-1 for x in new_chr]
    random.shuffle(chrlist)
    print 'chr =', chrlist
if args.verbose:                                          
    verbose = 1
    print 'verbose =', verbose    
if args.self_comparison:                                         
    self_comparison = 1
    print 'self_comparison =', self_comparison 
if args.out:
    out_dir = args.out
    if not out_dir.endswith('/'):
        out_dir = out_dir+'/'
    print 'out_dir =', out_dir
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
        print '    creating dir '+out_dir
if args.targets:
    cmd_arg = 'targets='  
    targets_file = args.targets
    print 'targets_file =', targets_file
if args.contam_numind:
    contam_numind = args.contam_numind
    contam_numind = contam_numind.split(',')
    contam_numind[0] = int(contam_numind[0])
    contam_numind[1] = int(contam_numind[1])
    print 'contam_numind =', contam_numind
#if args.writeAllData:   
#    writeAllData = int(args.writeAllData)
#    print 'writeAllData =', writeAllData
if args.ped:
    parsePedigreeFile(args.ped)
print ' '





# define individual IDs within specific superpopulations
ind_IDs = {}
ind_IDs['EUR'] = ['HG00096', 'HG00097', 'HG00099', 'HG00100', 'HG00101', 'HG00102', 'HG00103', 'HG00105', 'HG00106', 'HG00107', 'HG00108', 'HG00109', 'HG00110', 'HG00111', 'HG00112', 'HG00113', 'HG00114', 'HG00115', 'HG00116', 'HG00117', 'HG00118', 'HG00119', 'HG00120', 'HG00121', 'HG00122', 'HG00123', 'HG00125', 'HG00126', 'HG00127', 'HG00128', 'HG00129', 'HG00130', 'HG00131', 'HG00132', 'HG00133', 'HG00136', 'HG00137', 'HG00138', 'HG00139', 'HG00140', 'HG00141', 'HG00142', 'HG00143', 'HG00145', 'HG00146', 'HG00148', 'HG00149', 'HG00150', 'HG00151', 'HG00154', 'HG00155', 'HG00157', 'HG00158', 'HG00159', 'HG00160', 'HG00171', 'HG00173', 'HG00174', 'HG00176', 'HG00177', 'HG00178', 'HG00179', 'HG00180', 'HG00181', 'HG00182', 'HG00183', 'HG00185', 'HG00186', 'HG00187', 'HG00188', 'HG00189', 'HG00190', 'HG00231', 'HG00232', 'HG00233', 'HG00234', 'HG00235', 'HG00236', 'HG00237', 'HG00238', 'HG00239', 'HG00240', 'HG00242', 'HG00243', 'HG00244', 'HG00245', 'HG00246', 'HG00250', 'HG00251', 'HG00252', 'HG00253', 'HG00254', 'HG00255', 'HG00256', 'HG00257', 'HG00258', 'HG00259', 'HG00260', 'HG00261', 'HG00262', 'HG00263', 'HG00264', 'HG00265', 'HG00266', 'HG00267', 'HG00268', 'HG00269', 'HG00271', 'HG00272', 'HG00273', 'HG00274', 'HG00275', 'HG00276', 'HG00277', 'HG00278', 'HG00280', 'HG00281', 'HG00282', 'HG00284', 'HG00285', 'HG00288', 'HG00290', 'HG00304', 'HG00306', 'HG00308', 'HG00309', 'HG00310', 'HG00311', 'HG00313', 'HG00315', 'HG00318', 'HG00319', 'HG00320', 'HG00321', 'HG00323', 'HG00324', 'HG00325', 'HG00326', 'HG00327', 'HG00328', 'HG00329', 'HG00330', 'HG00331', 'HG00332', 'HG00334', 'HG00335', 'HG00336', 'HG00337', 'HG00338', 'HG00339', 'HG00341', 'HG00342', 'HG00343', 'HG00344', 'HG00345', 'HG00346', 'HG00349', 'HG00350', 'HG00351', 'HG00353', 'HG00355', 'HG00356', 'HG00357', 'HG00358', 'HG00360', 'HG00361', 'HG00362', 'HG00364', 'HG00365', 'HG00366', 'HG00367', 'HG00368', 'HG00369', 'HG00371', 'HG00372', 'HG00373', 'HG00375', 'HG00376', 'HG00378', 'HG00379', 'HG00380', 'HG00381', 'HG00382', 'HG00383', 'HG00384', 'HG01334', 'HG01500', 'HG01501', 'HG01503', 'HG01504', 'HG01506', 'HG01507', 'HG01509', 'HG01510', 'HG01512', 'HG01513', 'HG01515', 'HG01516', 'HG01518', 'HG01519', 'HG01521', 'HG01522', 'HG01524', 'HG01525', 'HG01527', 'HG01528', 'HG01530', 'HG01531', 'HG01536', 'HG01537', 'HG01602', 'HG01603', 'HG01605', 'HG01606', 'HG01607', 'HG01608', 'HG01610', 'HG01612', 'HG01613', 'HG01615', 'HG01617', 'HG01618', 'HG01619', 'HG01620', 'HG01623', 'HG01624', 'HG01625', 'HG01626', 'HG01628', 'HG01630', 'HG01631', 'HG01632', 'HG01668', 'HG01669', 'HG01670', 'HG01672', 'HG01673', 'HG01675', 'HG01676', 'HG01678', 'HG01679', 'HG01680', 'HG01682', 'HG01684', 'HG01685', 'HG01686', 'HG01694', 'HG01695', 'HG01697', 'HG01699', 'HG01700', 'HG01702', 'HG01704', 'HG01705', 'HG01707', 'HG01708', 'HG01709', 'HG01710', 'HG01746', 'HG01747', 'HG01756', 'HG01757', 'HG01761', 'HG01762', 'HG01765', 'HG01766', 'HG01767', 'HG01768', 'HG01770', 'HG01771', 'HG01773', 'HG01775', 'HG01776', 'HG01777', 'HG01779', 'HG01781', 'HG01783', 'HG01784', 'HG01785', 'HG01786', 'HG01789', 'HG01790', 'HG01791', 'HG02215', 'HG02219', 'HG02220', 'HG02221', 'HG02223', 'HG02224', 'HG02230', 'HG02231', 'HG02232', 'HG02233', 'HG02235', 'HG02236', 'HG02238', 'HG02239', 'NA06984', 'NA06985', 'NA06986', 'NA06989', 'NA06994', 'NA07000', 'NA07037', 'NA07048', 'NA07051', 'NA07056', 'NA07347', 'NA07357', 'NA10847', 'NA10851', 'NA11829', 'NA11830', 'NA11831', 'NA11832', 'NA11840', 'NA11843', 'NA11881', 'NA11892', 'NA11893', 'NA11894', 'NA11918', 'NA11919', 'NA11920', 'NA11930', 'NA11931', 'NA11932', 'NA11933', 'NA11992', 'NA11994', 'NA11995', 'NA12003', 'NA12004', 'NA12005', 'NA12006', 'NA12043', 'NA12044', 'NA12045', 'NA12046', 'NA12058', 'NA12144', 'NA12154', 'NA12155', 'NA12156', 'NA12234', 'NA12249', 'NA12272', 'NA12273', 'NA12275', 'NA12282', 'NA12283', 'NA12286', 'NA12287', 'NA12340', 'NA12341', 'NA12342', 'NA12347', 'NA12348', 'NA12383', 'NA12399', 'NA12400', 'NA12413', 'NA12414', 'NA12489', 'NA12546', 'NA12716', 'NA12717', 'NA12718', 'NA12748', 'NA12749', 'NA12750', 'NA12751', 'NA12760', 'NA12761', 'NA12762', 'NA12763', 'NA12775', 'NA12776', 'NA12777', 'NA12778', 'NA12812', 'NA12813', 'NA12814', 'NA12815', 'NA12827', 'NA12828', 'NA12829', 'NA12830', 'NA12842', 'NA12843', 'NA12872', 'NA12873', 'NA12874', 'NA12878', 'NA12889', 'NA12890', 'NA20502', 'NA20503', 'NA20504', 'NA20505', 'NA20506', 'NA20507', 'NA20508', 'NA20509', 'NA20510', 'NA20511', 'NA20512', 'NA20513', 'NA20514', 'NA20515', 'NA20516', 'NA20517', 'NA20518', 'NA20519', 'NA20520', 'NA20521', 'NA20522', 'NA20524', 'NA20525', 'NA20527', 'NA20528', 'NA20529', 'NA20530', 'NA20531', 'NA20532', 'NA20533', 'NA20534', 'NA20535', 'NA20536', 'NA20538', 'NA20539', 'NA20540', 'NA20541', 'NA20542', 'NA20543', 'NA20544', 'NA20581', 'NA20582', 'NA20585', 'NA20586', 'NA20587', 'NA20588', 'NA20589', 'NA20752', 'NA20753', 'NA20754', 'NA20755', 'NA20756', 'NA20757', 'NA20758', 'NA20759', 'NA20760', 'NA20761', 'NA20762', 'NA20763', 'NA20764', 'NA20765', 'NA20766', 'NA20767', 'NA20768', 'NA20769', 'NA20770', 'NA20771', 'NA20772', 'NA20773', 'NA20774', 'NA20775', 'NA20778', 'NA20783', 'NA20785', 'NA20786', 'NA20787', 'NA20790', 'NA20792', 'NA20795', 'NA20796', 'NA20797', 'NA20798', 'NA20799', 'NA20800', 'NA20801', 'NA20802', 'NA20803', 'NA20804', 'NA20805', 'NA20806', 'NA20807', 'NA20808', 'NA20809', 'NA20810', 'NA20811', 'NA20812', 'NA20813', 'NA20814', 'NA20815', 'NA20818', 'NA20819', 'NA20821', 'NA20822', 'NA20826', 'NA20827', 'NA20828', 'NA20832']
ind_IDs['EAS'] = ['HG00403', 'HG00404', 'HG00406', 'HG00407', 'HG00409', 'HG00410', 'HG00419', 'HG00421', 'HG00422', 'HG00428', 'HG00436', 'HG00437', 'HG00442', 'HG00443', 'HG00445', 'HG00446', 'HG00448', 'HG00449', 'HG00451', 'HG00452', 'HG00457', 'HG00458', 'HG00463', 'HG00464', 'HG00472', 'HG00473', 'HG00475', 'HG00476', 'HG00478', 'HG00479', 'HG00500', 'HG00513', 'HG00524', 'HG00525', 'HG00530', 'HG00531', 'HG00533', 'HG00534', 'HG00536', 'HG00537', 'HG00542', 'HG00543', 'HG00556', 'HG00557', 'HG00559', 'HG00560', 'HG00565', 'HG00566', 'HG00580', 'HG00581', 'HG00583', 'HG00584', 'HG00589', 'HG00590', 'HG00592', 'HG00593', 'HG00595', 'HG00596', 'HG00598', 'HG00599', 'HG00607', 'HG00608', 'HG00610', 'HG00611', 'HG00613', 'HG00614', 'HG00619', 'HG00620', 'HG00622', 'HG00623', 'HG00625', 'HG00626', 'HG00628', 'HG00629', 'HG00631', 'HG00632', 'HG00634', 'HG00650', 'HG00651', 'HG00653', 'HG00654', 'HG00656', 'HG00657', 'HG00662', 'HG00663', 'HG00671', 'HG00672', 'HG00674', 'HG00675', 'HG00683', 'HG00684', 'HG00689', 'HG00690', 'HG00692', 'HG00693', 'HG00698', 'HG00699', 'HG00701', 'HG00704', 'HG00705', 'HG00707', 'HG00708', 'HG00717', 'HG00728', 'HG00729', 'HG00759', 'HG00766', 'HG00844', 'HG00851', 'HG00864', 'HG00867', 'HG00879', 'HG00881', 'HG00956', 'HG00978', 'HG00982', 'HG01028', 'HG01029', 'HG01031', 'HG01046', 'HG01595', 'HG01596', 'HG01597', 'HG01598', 'HG01599', 'HG01600', 'HG01794', 'HG01795', 'HG01796', 'HG01797', 'HG01798', 'HG01799', 'HG01800', 'HG01801', 'HG01802', 'HG01804', 'HG01805', 'HG01806', 'HG01807', 'HG01808', 'HG01809', 'HG01810', 'HG01811', 'HG01812', 'HG01813', 'HG01815', 'HG01816', 'HG01817', 'HG01840', 'HG01841', 'HG01842', 'HG01843', 'HG01844', 'HG01845', 'HG01846', 'HG01847', 'HG01848', 'HG01849', 'HG01850', 'HG01851', 'HG01852', 'HG01853', 'HG01855', 'HG01857', 'HG01858', 'HG01859', 'HG01860', 'HG01861', 'HG01862', 'HG01863', 'HG01864', 'HG01865', 'HG01866', 'HG01867', 'HG01868', 'HG01869', 'HG01870', 'HG01871', 'HG01872', 'HG01873', 'HG01874', 'HG01878', 'HG02016', 'HG02017', 'HG02019', 'HG02020', 'HG02023', 'HG02025', 'HG02026', 'HG02028', 'HG02029', 'HG02031', 'HG02032', 'HG02035', 'HG02040', 'HG02047', 'HG02048', 'HG02049', 'HG02050', 'HG02057', 'HG02058', 'HG02060', 'HG02061', 'HG02064', 'HG02067', 'HG02069', 'HG02070', 'HG02072', 'HG02073', 'HG02075', 'HG02076', 'HG02078', 'HG02079', 'HG02081', 'HG02082', 'HG02084', 'HG02085', 'HG02086', 'HG02087', 'HG02088', 'HG02113', 'HG02116', 'HG02121', 'HG02122', 'HG02127', 'HG02128', 'HG02130', 'HG02131', 'HG02133', 'HG02134', 'HG02136', 'HG02137', 'HG02138', 'HG02139', 'HG02140', 'HG02141', 'HG02142', 'HG02151', 'HG02152', 'HG02153', 'HG02154', 'HG02155', 'HG02156', 'HG02164', 'HG02165', 'HG02166', 'HG02178', 'HG02179', 'HG02180', 'HG02181', 'HG02182', 'HG02184', 'HG02185', 'HG02186', 'HG02187', 'HG02188', 'HG02190', 'HG02250', 'HG02351', 'HG02353', 'HG02355', 'HG02356', 'HG02360', 'HG02364', 'HG02367', 'HG02371', 'HG02373', 'HG02374', 'HG02375', 'HG02379', 'HG02380', 'HG02382', 'HG02383', 'HG02384', 'HG02385', 'HG02386', 'HG02389', 'HG02390', 'HG02391', 'HG02392', 'HG02394', 'HG02395', 'HG02396', 'HG02397', 'HG02398', 'HG02399', 'HG02401', 'HG02402', 'HG02406', 'HG02407', 'HG02408', 'HG02409', 'HG02410', 'HG02512', 'HG02513', 'HG02521', 'HG02522', 'NA18525', 'NA18526', 'NA18528', 'NA18530', 'NA18531', 'NA18532', 'NA18533', 'NA18534', 'NA18535', 'NA18536', 'NA18537', 'NA18538', 'NA18539', 'NA18541', 'NA18542', 'NA18543', 'NA18544', 'NA18545', 'NA18546', 'NA18547', 'NA18548', 'NA18549', 'NA18550', 'NA18552', 'NA18553', 'NA18555', 'NA18557', 'NA18558', 'NA18559', 'NA18560', 'NA18561', 'NA18562', 'NA18563', 'NA18564', 'NA18565', 'NA18566', 'NA18567', 'NA18570', 'NA18571', 'NA18572', 'NA18573', 'NA18574', 'NA18577', 'NA18579', 'NA18582', 'NA18591', 'NA18592', 'NA18593', 'NA18595', 'NA18596', 'NA18597', 'NA18599', 'NA18602', 'NA18603', 'NA18605', 'NA18606', 'NA18608', 'NA18609', 'NA18610', 'NA18611', 'NA18612', 'NA18613', 'NA18614', 'NA18615', 'NA18616', 'NA18617', 'NA18618', 'NA18619', 'NA18620', 'NA18621', 'NA18622', 'NA18623', 'NA18624', 'NA18625', 'NA18626', 'NA18627', 'NA18628', 'NA18629', 'NA18630', 'NA18631', 'NA18632', 'NA18633', 'NA18634', 'NA18635', 'NA18636', 'NA18637', 'NA18638', 'NA18639', 'NA18640', 'NA18641', 'NA18642', 'NA18643', 'NA18644', 'NA18645', 'NA18646', 'NA18647', 'NA18648', 'NA18740', 'NA18745', 'NA18747', 'NA18748', 'NA18749', 'NA18757', 'NA18939', 'NA18940', 'NA18941', 'NA18942', 'NA18943', 'NA18944', 'NA18945', 'NA18946', 'NA18947', 'NA18948', 'NA18949', 'NA18950', 'NA18951', 'NA18952', 'NA18953', 'NA18954', 'NA18956', 'NA18957', 'NA18959', 'NA18960', 'NA18961', 'NA18962', 'NA18963', 'NA18964', 'NA18965', 'NA18966', 'NA18967', 'NA18968', 'NA18969', 'NA18970', 'NA18971', 'NA18972', 'NA18973', 'NA18974', 'NA18975', 'NA18976', 'NA18977', 'NA18978', 'NA18979', 'NA18980', 'NA18981', 'NA18982', 'NA18983', 'NA18984', 'NA18985', 'NA18986', 'NA18987', 'NA18988', 'NA18989', 'NA18990', 'NA18991', 'NA18992', 'NA18993', 'NA18994', 'NA18995', 'NA18997', 'NA18998', 'NA18999', 'NA19000', 'NA19001', 'NA19002', 'NA19003', 'NA19004', 'NA19005', 'NA19006', 'NA19007', 'NA19009', 'NA19010', 'NA19011', 'NA19012', 'NA19054', 'NA19055', 'NA19056', 'NA19057', 'NA19058', 'NA19059', 'NA19060', 'NA19062', 'NA19063', 'NA19064', 'NA19065', 'NA19066', 'NA19067', 'NA19068', 'NA19070', 'NA19072', 'NA19074', 'NA19075', 'NA19076', 'NA19077', 'NA19078', 'NA19079', 'NA19080', 'NA19081', 'NA19082', 'NA19083', 'NA19084', 'NA19085', 'NA19086', 'NA19087', 'NA19088', 'NA19089', 'NA19090', 'NA19091']
ind_IDs['AFR'] = ['HG01879', 'HG01880', 'HG01882', 'HG01883', 'HG01885', 'HG01886', 'HG01889', 'HG01890', 'HG01894', 'HG01896', 'HG01912', 'HG01914', 'HG01915', 'HG01956', 'HG01958', 'HG01985', 'HG01986', 'HG01988', 'HG01989', 'HG01990', 'HG02009', 'HG02010', 'HG02012', 'HG02013', 'HG02014', 'HG02051', 'HG02052', 'HG02053', 'HG02054', 'HG02095', 'HG02107', 'HG02108', 'HG02111', 'HG02143', 'HG02144', 'HG02255', 'HG02256', 'HG02281', 'HG02282', 'HG02283', 'HG02284', 'HG02307', 'HG02308', 'HG02309', 'HG02314', 'HG02315', 'HG02317', 'HG02318', 'HG02322', 'HG02323', 'HG02325', 'HG02330', 'HG02332', 'HG02334', 'HG02337', 'HG02339', 'HG02343', 'HG02419', 'HG02420', 'HG02427', 'HG02429', 'HG02433', 'HG02439', 'HG02442', 'HG02445', 'HG02449', 'HG02450', 'HG02455', 'HG02461', 'HG02462', 'HG02464', 'HG02465', 'HG02470', 'HG02471', 'HG02476', 'HG02477', 'HG02479', 'HG02481', 'HG02484', 'HG02485', 'HG02489', 'HG02496', 'HG02497', 'HG02501', 'HG02502', 'HG02505', 'HG02508', 'HG02511', 'HG02536', 'HG02537', 'HG02541', 'HG02545', 'HG02546', 'HG02549', 'HG02554', 'HG02555', 'HG02557', 'HG02558', 'HG02561', 'HG02562', 'HG02568', 'HG02570', 'HG02571', 'HG02573', 'HG02574', 'HG02577', 'HG02580', 'HG02582', 'HG02583', 'HG02585', 'HG02586', 'HG02588', 'HG02589', 'HG02594', 'HG02595', 'HG02610', 'HG02611', 'HG02613', 'HG02614', 'HG02620', 'HG02621', 'HG02623', 'HG02624', 'HG02628', 'HG02629', 'HG02634', 'HG02635', 'HG02642', 'HG02643', 'HG02645', 'HG02646', 'HG02666', 'HG02667', 'HG02675', 'HG02676', 'HG02678', 'HG02679', 'HG02702', 'HG02703', 'HG02715', 'HG02716', 'HG02721', 'HG02722', 'HG02756', 'HG02757', 'HG02759', 'HG02760', 'HG02763', 'HG02768', 'HG02769', 'HG02771', 'HG02772', 'HG02798', 'HG02799', 'HG02804', 'HG02805', 'HG02807', 'HG02808', 'HG02810', 'HG02811', 'HG02813', 'HG02814', 'HG02816', 'HG02817', 'HG02819', 'HG02820', 'HG02836', 'HG02837', 'HG02839', 'HG02840', 'HG02851', 'HG02852', 'HG02854', 'HG02855', 'HG02860', 'HG02861', 'HG02870', 'HG02878', 'HG02879', 'HG02881', 'HG02882', 'HG02884', 'HG02885', 'HG02887', 'HG02888', 'HG02890', 'HG02891', 'HG02895', 'HG02896', 'HG02922', 'HG02923', 'HG02938', 'HG02941', 'HG02943', 'HG02944', 'HG02946', 'HG02947', 'HG02952', 'HG02953', 'HG02968', 'HG02970', 'HG02971', 'HG02973', 'HG02974', 'HG02976', 'HG02977', 'HG02979', 'HG02981', 'HG02982', 'HG02983', 'HG03024', 'HG03025', 'HG03027', 'HG03028', 'HG03039', 'HG03040', 'HG03045', 'HG03046', 'HG03048', 'HG03049', 'HG03052', 'HG03054', 'HG03055', 'HG03057', 'HG03058', 'HG03060', 'HG03061', 'HG03063', 'HG03064', 'HG03066', 'HG03069', 'HG03072', 'HG03073', 'HG03074', 'HG03077', 'HG03078', 'HG03079', 'HG03081', 'HG03082', 'HG03084', 'HG03085', 'HG03086', 'HG03088', 'HG03091', 'HG03095', 'HG03096', 'HG03097', 'HG03099', 'HG03100', 'HG03103', 'HG03105', 'HG03108', 'HG03109', 'HG03111', 'HG03112', 'HG03114', 'HG03115', 'HG03117', 'HG03118', 'HG03120', 'HG03121', 'HG03123', 'HG03124', 'HG03126', 'HG03127', 'HG03129', 'HG03130', 'HG03132', 'HG03133', 'HG03135', 'HG03136', 'HG03139', 'HG03157', 'HG03159', 'HG03160', 'HG03162', 'HG03163', 'HG03166', 'HG03168', 'HG03169', 'HG03172', 'HG03175', 'HG03189', 'HG03190', 'HG03193', 'HG03195', 'HG03196', 'HG03198', 'HG03199', 'HG03202', 'HG03209', 'HG03212', 'HG03224', 'HG03225', 'HG03240', 'HG03241', 'HG03246', 'HG03247', 'HG03258', 'HG03259', 'HG03265', 'HG03267', 'HG03268', 'HG03270', 'HG03271', 'HG03279', 'HG03280', 'HG03291', 'HG03294', 'HG03295', 'HG03297', 'HG03298', 'HG03300', 'HG03301', 'HG03303', 'HG03304', 'HG03311', 'HG03313', 'HG03342', 'HG03343', 'HG03351', 'HG03352', 'HG03354', 'HG03363', 'HG03366', 'HG03367', 'HG03369', 'HG03370', 'HG03372', 'HG03376', 'HG03378', 'HG03380', 'HG03382', 'HG03385', 'HG03388', 'HG03391', 'HG03394', 'HG03397', 'HG03401', 'HG03410', 'HG03419', 'HG03428', 'HG03432', 'HG03433', 'HG03436', 'HG03437', 'HG03439', 'HG03442', 'HG03445', 'HG03446', 'HG03449', 'HG03451', 'HG03452', 'HG03455', 'HG03457', 'HG03458', 'HG03460', 'HG03461', 'HG03464', 'HG03469', 'HG03470', 'HG03472', 'HG03473', 'HG03476', 'HG03478', 'HG03479', 'HG03484', 'HG03485', 'HG03499', 'HG03511', 'HG03514', 'HG03515', 'HG03517', 'HG03518', 'HG03520', 'HG03521', 'HG03538', 'HG03539', 'HG03547', 'HG03548', 'HG03556', 'HG03557', 'HG03558', 'HG03559', 'HG03563', 'HG03565', 'HG03567', 'HG03571', 'HG03572', 'HG03575', 'HG03577', 'HG03578', 'HG03583', 'NA18486', 'NA18488', 'NA18489', 'NA18498', 'NA18499', 'NA18501', 'NA18502', 'NA18504', 'NA18505', 'NA18507', 'NA18508', 'NA18510', 'NA18511', 'NA18516', 'NA18517', 'NA18519', 'NA18520', 'NA18522', 'NA18523', 'NA18853', 'NA18856', 'NA18858', 'NA18861', 'NA18864', 'NA18865', 'NA18867', 'NA18868', 'NA18870', 'NA18871', 'NA18873', 'NA18874', 'NA18876', 'NA18877', 'NA18878', 'NA18879', 'NA18881', 'NA18907', 'NA18908', 'NA18909', 'NA18910', 'NA18912', 'NA18915', 'NA18916', 'NA18917', 'NA18923', 'NA18924', 'NA18933', 'NA18934', 'NA19017', 'NA19019', 'NA19020', 'NA19023', 'NA19024', 'NA19025', 'NA19026', 'NA19027', 'NA19028', 'NA19030', 'NA19031', 'NA19035', 'NA19036', 'NA19037', 'NA19038', 'NA19041', 'NA19042', 'NA19043', 'NA19092', 'NA19093', 'NA19095', 'NA19096', 'NA19098', 'NA19099', 'NA19102', 'NA19107', 'NA19108', 'NA19113', 'NA19114', 'NA19116', 'NA19117', 'NA19118', 'NA19119', 'NA19121', 'NA19129', 'NA19130', 'NA19131', 'NA19137', 'NA19138', 'NA19141', 'NA19143', 'NA19144', 'NA19146', 'NA19147', 'NA19149', 'NA19152', 'NA19153', 'NA19159', 'NA19160', 'NA19171', 'NA19172', 'NA19175', 'NA19184', 'NA19185', 'NA19189', 'NA19190', 'NA19197', 'NA19198', 'NA19200', 'NA19201', 'NA19204', 'NA19206', 'NA19207', 'NA19209', 'NA19210', 'NA19213', 'NA19214', 'NA19222', 'NA19223', 'NA19225', 'NA19235', 'NA19236', 'NA19238', 'NA19239', 'NA19247', 'NA19248', 'NA19256', 'NA19257', 'NA19307', 'NA19308', 'NA19309', 'NA19310', 'NA19312', 'NA19314', 'NA19315', 'NA19316', 'NA19317', 'NA19318', 'NA19319', 'NA19320', 'NA19321', 'NA19323', 'NA19324', 'NA19327', 'NA19328', 'NA19331', 'NA19332', 'NA19334', 'NA19338', 'NA19346', 'NA19347', 'NA19350', 'NA19351', 'NA19355', 'NA19360', 'NA19372', 'NA19374', 'NA19375', 'NA19376', 'NA19377', 'NA19378', 'NA19379', 'NA19380', 'NA19383', 'NA19384', 'NA19385', 'NA19390', 'NA19391', 'NA19393', 'NA19394', 'NA19395', 'NA19397', 'NA19399', 'NA19401', 'NA19403', 'NA19404', 'NA19428', 'NA19429', 'NA19430', 'NA19431', 'NA19434', 'NA19435', 'NA19436', 'NA19437', 'NA19438', 'NA19439', 'NA19440', 'NA19443', 'NA19445', 'NA19446', 'NA19448', 'NA19449', 'NA19451', 'NA19452', 'NA19454', 'NA19455', 'NA19456', 'NA19457', 'NA19461', 'NA19462', 'NA19463', 'NA19466', 'NA19467', 'NA19468', 'NA19471', 'NA19472', 'NA19473', 'NA19474', 'NA19475', 'NA19625', 'NA19700', 'NA19701', 'NA19703', 'NA19704', 'NA19707', 'NA19711', 'NA19712', 'NA19713', 'NA19818', 'NA19819', 'NA19834', 'NA19835', 'NA19900', 'NA19901', 'NA19904', 'NA19908', 'NA19909', 'NA19913', 'NA19914', 'NA19916', 'NA19917', 'NA19920', 'NA19921', 'NA19922', 'NA19923', 'NA19982', 'NA19984', 'NA20126', 'NA20127', 'NA20274', 'NA20276', 'NA20278', 'NA20281', 'NA20282', 'NA20287', 'NA20289', 'NA20291', 'NA20294', 'NA20296', 'NA20298', 'NA20299', 'NA20314', 'NA20317', 'NA20318', 'NA20320', 'NA20321', 'NA20332', 'NA20334', 'NA20339', 'NA20340', 'NA20342', 'NA20346', 'NA20348', 'NA20351', 'NA20355', 'NA20356', 'NA20357', 'NA20359', 'NA20362', 'NA20412']
ind_IDs['AMR'] = ['HG00551', 'HG00553', 'HG00554', 'HG00637', 'HG00638', 'HG00640', 'HG00641', 'HG00731', 'HG00732', 'HG00734', 'HG00736', 'HG00737', 'HG00739', 'HG00740', 'HG00742', 'HG00743', 'HG01047', 'HG01048', 'HG01049', 'HG01051', 'HG01052', 'HG01054', 'HG01055', 'HG01058', 'HG01060', 'HG01061', 'HG01063', 'HG01064', 'HG01066', 'HG01067', 'HG01069', 'HG01070', 'HG01072', 'HG01073', 'HG01075', 'HG01077', 'HG01079', 'HG01080', 'HG01082', 'HG01083', 'HG01085', 'HG01086', 'HG01088', 'HG01089', 'HG01092', 'HG01094', 'HG01095', 'HG01097', 'HG01098', 'HG01101', 'HG01102', 'HG01104', 'HG01105', 'HG01107', 'HG01108', 'HG01110', 'HG01111', 'HG01112', 'HG01113', 'HG01119', 'HG01121', 'HG01122', 'HG01124', 'HG01125', 'HG01130', 'HG01131', 'HG01133', 'HG01134', 'HG01136', 'HG01137', 'HG01139', 'HG01140', 'HG01142', 'HG01148', 'HG01149', 'HG01161', 'HG01162', 'HG01164', 'HG01167', 'HG01168', 'HG01170', 'HG01171', 'HG01173', 'HG01174', 'HG01176', 'HG01177', 'HG01182', 'HG01183', 'HG01187', 'HG01188', 'HG01190', 'HG01191', 'HG01197', 'HG01198', 'HG01200', 'HG01204', 'HG01205', 'HG01241', 'HG01242', 'HG01247', 'HG01248', 'HG01250', 'HG01251', 'HG01253', 'HG01254', 'HG01256', 'HG01257', 'HG01259', 'HG01260', 'HG01269', 'HG01271', 'HG01272', 'HG01275', 'HG01277', 'HG01280', 'HG01281', 'HG01284', 'HG01286', 'HG01302', 'HG01303', 'HG01305', 'HG01308', 'HG01311', 'HG01312', 'HG01323', 'HG01325', 'HG01326', 'HG01341', 'HG01342', 'HG01344', 'HG01345', 'HG01348', 'HG01350', 'HG01351', 'HG01353', 'HG01354', 'HG01356', 'HG01357', 'HG01359', 'HG01360', 'HG01362', 'HG01363', 'HG01365', 'HG01366', 'HG01369', 'HG01372', 'HG01374', 'HG01375', 'HG01377', 'HG01378', 'HG01383', 'HG01384', 'HG01389', 'HG01390', 'HG01392', 'HG01393', 'HG01395', 'HG01396', 'HG01398', 'HG01402', 'HG01403', 'HG01405', 'HG01412', 'HG01413', 'HG01414', 'HG01431', 'HG01432', 'HG01435', 'HG01437', 'HG01438', 'HG01440', 'HG01441', 'HG01443', 'HG01444', 'HG01447', 'HG01455', 'HG01456', 'HG01459', 'HG01461', 'HG01462', 'HG01464', 'HG01465', 'HG01468', 'HG01474', 'HG01479', 'HG01485', 'HG01486', 'HG01488', 'HG01489', 'HG01491', 'HG01492', 'HG01494', 'HG01495', 'HG01497', 'HG01498', 'HG01550', 'HG01551', 'HG01556', 'HG01565', 'HG01566', 'HG01571', 'HG01572', 'HG01577', 'HG01578', 'HG01892', 'HG01893', 'HG01917', 'HG01918', 'HG01920', 'HG01921', 'HG01923', 'HG01924', 'HG01926', 'HG01927', 'HG01932', 'HG01933', 'HG01935', 'HG01936', 'HG01938', 'HG01939', 'HG01941', 'HG01942', 'HG01944', 'HG01945', 'HG01947', 'HG01948', 'HG01950', 'HG01951', 'HG01953', 'HG01954', 'HG01961', 'HG01965', 'HG01967', 'HG01968', 'HG01970', 'HG01971', 'HG01973', 'HG01974', 'HG01976', 'HG01977', 'HG01979', 'HG01980', 'HG01982', 'HG01991', 'HG01992', 'HG01997', 'HG02002', 'HG02003', 'HG02006', 'HG02008', 'HG02089', 'HG02090', 'HG02102', 'HG02104', 'HG02105', 'HG02146', 'HG02147', 'HG02150', 'HG02252', 'HG02253', 'HG02259', 'HG02260', 'HG02262', 'HG02265', 'HG02266', 'HG02271', 'HG02272', 'HG02274', 'HG02275', 'HG02277', 'HG02278', 'HG02285', 'HG02286', 'HG02291', 'HG02292', 'HG02298', 'HG02299', 'HG02301', 'HG02304', 'HG02312', 'HG02345', 'HG02348', 'HG02425', 'NA19648', 'NA19649', 'NA19651', 'NA19652', 'NA19654', 'NA19655', 'NA19657', 'NA19658', 'NA19661', 'NA19663', 'NA19664', 'NA19669', 'NA19670', 'NA19676', 'NA19678', 'NA19679', 'NA19681', 'NA19682', 'NA19684', 'NA19716', 'NA19717', 'NA19719', 'NA19720', 'NA19722', 'NA19723', 'NA19725', 'NA19726', 'NA19728', 'NA19729', 'NA19731', 'NA19732', 'NA19734', 'NA19735', 'NA19740', 'NA19741', 'NA19746', 'NA19747', 'NA19749', 'NA19750', 'NA19752', 'NA19755', 'NA19756', 'NA19758', 'NA19759', 'NA19761', 'NA19762', 'NA19764', 'NA19770', 'NA19771', 'NA19773', 'NA19774', 'NA19776', 'NA19777', 'NA19779', 'NA19780', 'NA19782', 'NA19783', 'NA19785', 'NA19786', 'NA19788', 'NA19789', 'NA19792', 'NA19794', 'NA19795']
ind_IDs['SAS'] = ['HG01583', 'HG01586', 'HG01589', 'HG01593', 'HG02490', 'HG02491', 'HG02493', 'HG02494', 'HG02597', 'HG02600', 'HG02601', 'HG02603', 'HG02604', 'HG02648', 'HG02649', 'HG02651', 'HG02652', 'HG02654', 'HG02655', 'HG02657', 'HG02658', 'HG02660', 'HG02661', 'HG02681', 'HG02682', 'HG02684', 'HG02685', 'HG02687', 'HG02688', 'HG02690', 'HG02691', 'HG02694', 'HG02696', 'HG02697', 'HG02699', 'HG02700', 'HG02724', 'HG02725', 'HG02727', 'HG02728', 'HG02731', 'HG02733', 'HG02734', 'HG02736', 'HG02737', 'HG02774', 'HG02775', 'HG02778', 'HG02780', 'HG02783', 'HG02784', 'HG02786', 'HG02787', 'HG02789', 'HG02790', 'HG02792', 'HG02793', 'HG03006', 'HG03007', 'HG03009', 'HG03012', 'HG03015', 'HG03016', 'HG03018', 'HG03019', 'HG03021', 'HG03022', 'HG03228', 'HG03229', 'HG03234', 'HG03235', 'HG03237', 'HG03238', 'HG03488', 'HG03490', 'HG03491', 'HG03585', 'HG03589', 'HG03593', 'HG03594', 'HG03595', 'HG03598', 'HG03600', 'HG03603', 'HG03604', 'HG03607', 'HG03611', 'HG03615', 'HG03616', 'HG03619', 'HG03624', 'HG03625', 'HG03629', 'HG03631', 'HG03634', 'HG03636', 'HG03640', 'HG03642', 'HG03643', 'HG03644', 'HG03645', 'HG03646', 'HG03649', 'HG03652', 'HG03653', 'HG03660', 'HG03663', 'HG03667', 'HG03668', 'HG03672', 'HG03673', 'HG03679', 'HG03680', 'HG03681', 'HG03684', 'HG03685', 'HG03686', 'HG03687', 'HG03689', 'HG03690', 'HG03691', 'HG03692', 'HG03693', 'HG03694', 'HG03695', 'HG03696', 'HG03697', 'HG03698', 'HG03702', 'HG03703', 'HG03705', 'HG03706', 'HG03708', 'HG03709', 'HG03711', 'HG03713', 'HG03714', 'HG03716', 'HG03717', 'HG03718', 'HG03720', 'HG03722', 'HG03727', 'HG03729', 'HG03730', 'HG03731', 'HG03733', 'HG03736', 'HG03738', 'HG03740', 'HG03741', 'HG03742', 'HG03743', 'HG03744', 'HG03745', 'HG03746', 'HG03750', 'HG03752', 'HG03753', 'HG03754', 'HG03755', 'HG03756', 'HG03757', 'HG03760', 'HG03762', 'HG03765', 'HG03767', 'HG03770', 'HG03771', 'HG03772', 'HG03773', 'HG03774', 'HG03775', 'HG03777', 'HG03778', 'HG03779', 'HG03780', 'HG03781', 'HG03782', 'HG03784', 'HG03785', 'HG03786', 'HG03787', 'HG03788', 'HG03789', 'HG03790', 'HG03792', 'HG03793', 'HG03796', 'HG03800', 'HG03802', 'HG03803', 'HG03805', 'HG03808', 'HG03809', 'HG03812', 'HG03814', 'HG03815', 'HG03817', 'HG03821', 'HG03823', 'HG03824', 'HG03826', 'HG03829', 'HG03830', 'HG03832', 'HG03833', 'HG03836', 'HG03837', 'HG03838', 'HG03844', 'HG03846', 'HG03848', 'HG03849', 'HG03850', 'HG03851', 'HG03854', 'HG03856', 'HG03857', 'HG03858', 'HG03861', 'HG03862', 'HG03863', 'HG03864', 'HG03866', 'HG03867', 'HG03868', 'HG03869', 'HG03870', 'HG03871', 'HG03872', 'HG03873', 'HG03874', 'HG03875', 'HG03882', 'HG03884', 'HG03885', 'HG03886', 'HG03887', 'HG03888', 'HG03890', 'HG03894', 'HG03895', 'HG03896', 'HG03897', 'HG03898', 'HG03899', 'HG03900', 'HG03902', 'HG03905', 'HG03907', 'HG03908', 'HG03910', 'HG03911', 'HG03913', 'HG03914', 'HG03916', 'HG03917', 'HG03919', 'HG03920', 'HG03922', 'HG03925', 'HG03926', 'HG03928', 'HG03931', 'HG03934', 'HG03937', 'HG03940', 'HG03941', 'HG03943', 'HG03944', 'HG03945', 'HG03947', 'HG03949', 'HG03950', 'HG03951', 'HG03953', 'HG03955', 'HG03960', 'HG03963', 'HG03965', 'HG03967', 'HG03968', 'HG03969', 'HG03971', 'HG03973', 'HG03974', 'HG03976', 'HG03977', 'HG03978', 'HG03985', 'HG03986', 'HG03989', 'HG03990', 'HG03991', 'HG03995', 'HG03998', 'HG03999', 'HG04001', 'HG04002', 'HG04003', 'HG04006', 'HG04014', 'HG04015', 'HG04017', 'HG04018', 'HG04019', 'HG04020', 'HG04022', 'HG04023', 'HG04025', 'HG04026', 'HG04029', 'HG04033', 'HG04035', 'HG04038', 'HG04039', 'HG04042', 'HG04047', 'HG04054', 'HG04056', 'HG04059', 'HG04060', 'HG04061', 'HG04062', 'HG04063', 'HG04070', 'HG04075', 'HG04076', 'HG04080', 'HG04090', 'HG04093', 'HG04094', 'HG04096', 'HG04098', 'HG04099', 'HG04100', 'HG04106', 'HG04107', 'HG04118', 'HG04131', 'HG04134', 'HG04140', 'HG04141', 'HG04144', 'HG04146', 'HG04152', 'HG04153', 'HG04155', 'HG04156', 'HG04158', 'HG04159', 'HG04161', 'HG04162', 'HG04164', 'HG04171', 'HG04173', 'HG04176', 'HG04177', 'HG04180', 'HG04182', 'HG04183', 'HG04185', 'HG04186', 'HG04188', 'HG04189', 'HG04194', 'HG04195', 'HG04198', 'HG04200', 'HG04202', 'HG04206', 'HG04209', 'HG04210', 'HG04211', 'HG04212', 'HG04214', 'HG04216', 'HG04219', 'HG04222', 'HG04225', 'HG04227', 'HG04229', 'HG04235', 'HG04238', 'HG04239', 'NA20845', 'NA20846', 'NA20847', 'NA20849', 'NA20850', 'NA20851', 'NA20852', 'NA20853', 'NA20854', 'NA20856', 'NA20858', 'NA20859', 'NA20861', 'NA20862', 'NA20863', 'NA20864', 'NA20866', 'NA20867', 'NA20868', 'NA20869', 'NA20870', 'NA20872', 'NA20874', 'NA20875', 'NA20876', 'NA20877', 'NA20878', 'NA20881', 'NA20882', 'NA20884', 'NA20885', 'NA20886', 'NA20887', 'NA20888', 'NA20889', 'NA20890', 'NA20891', 'NA20892', 'NA20894', 'NA20895', 'NA20896', 'NA20897', 'NA20899', 'NA20900', 'NA20901', 'NA20902', 'NA20903', 'NA20904', 'NA20905', 'NA20906', 'NA20908', 'NA20910', 'NA20911', 'NA21086', 'NA21087', 'NA21088', 'NA21089', 'NA21090', 'NA21091', 'NA21092', 'NA21093', 'NA21094', 'NA21095', 'NA21097', 'NA21098', 'NA21099', 'NA21100', 'NA21101', 'NA21102', 'NA21103', 'NA21104', 'NA21105', 'NA21106', 'NA21107', 'NA21108', 'NA21109', 'NA21110', 'NA21111', 'NA21112', 'NA21113', 'NA21114', 'NA21115', 'NA21116', 'NA21117', 'NA21118', 'NA21119', 'NA21120', 'NA21122', 'NA21123', 'NA21124', 'NA21125', 'NA21126', 'NA21127', 'NA21128', 'NA21129', 'NA21130', 'NA21133', 'NA21135', 'NA21137', 'NA21141', 'NA21142', 'NA21143', 'NA21144']

# find individual IDs for pedigree and contaminating superpopulations
for p in ['EUR', 'EAS', 'AFR', 'AMR', 'SAS']:
    if pedigree_pop == p:
        pedigree_pop_IDs = ind_IDs[p]
    if contam_pop[0] == p:
        contam_pop0_IDs = ind_IDs[p]
    if contam_pop[1] == p:
        contam_pop1_IDs = ind_IDs[p]


# define paths to VCF source files
if chr_type == 'trunc':
    input_chromosome_file = [data_dir+'ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr2.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr3.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr5.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr6.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr7.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr8.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr9.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr10.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr11.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr12.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr13.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr14.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr15.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr16.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr17.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr18.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr19.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr21.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz', \
                             data_dir+'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.first_50k.vcf.gz']
    chrom_length = [1916348, 1504249, 1361858, 1389101, 1366007, 1649527, 1209368, 1165865, 1234702, 1427390, 1523801, 1689733, 20469242, 21038758, 23039206, 1216505, 1392279, 1650616, 1413255, 1686988, 15209071, 18089371]
    chrom_numvars = [49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750, 49750]
elif chr_type == 'whole':
    # replaced v5 data with v5a, changing file names from ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz to ALL.chr*.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz
    # Notes about update from version v5 to v5a (Feb. 20th, 2015)
    # Some additional annotations are added to the VCF files as listed below. No sites are changed.
    # 1. added rs numbers provided by dbSNP to the ID column for SNPs (lastest rs numbers are not available for chrX, Y and the 3M patched up chunk on chr12 yet) 
    # 2. added esv accessions provided by DGVA to the ID columns for SVs
    # 3. added variant type to the INFO column (VT=SNP, VT=SV, VT=MNP)
    # 4. added EX_TARGET and MULTI_ALLELIC flag to the INFO column
    # 5. striped off SVLEN for all types of SV except MEIs. This was suggested by the SV group as the current SVLEN calculation was not consistent and could be misleading. For MEIs, as there isn't INFO:END, SVLEN is kept to give the length of the insertions.
    # ...
    input_chromosome_file = [data_dir+'ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr2.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr3.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr5.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr6.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr7.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr8.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr9.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr10.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr11.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr12.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr13.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr14.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr15.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr16.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr17.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr18.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr19.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr21.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz', \
                             data_dir+'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz']
    chrom_length =  [249250621, 243199373,  198022430,  191154276,  180915260,  171115067,  159138663,  146364022,  141213431,  135534747,  135006516,  133851895,  115169878,  107349540,  102531392,  90354753,   81195210,   78077248,   59128983,   63025520,   48129895,   51304566]
    chrom_numvars = [6437262,   7047141,    5803225,    5702765,    5238706,    4997829,    4692916,    4572517,    3542185,    3971844,    4024959,    3848775,    2843476,    2641439,    2412457,    2682658,    2317399,    2255683,    1822254,    1803869,    1099164,    1097199]
else:
    print 'ERROR! chromosomal input not defined.'

# define path to map of non-uniform recombination rates
# data from: http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/2011-01_phaseII_B37/
# later obtained from here when hapmap website disappeared: https://github.com/johnbowes/CRAFT-GP/blob/master/source_data/genetic_map_HapMapII_GRCh37/genetic_map_GRCh37_chr*.txt
# International HapMap Consortium. 2007. A second generation human haplotype map of over 3.1 million SNPs. Nature 449, 851-861.
if recomb_dir != '':
    recomb_rates_file = [recomb_dir+'genetic_map_GRCh37_chr1.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr2.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr3.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr4.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr5.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr6.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr7.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr8.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr9.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr10.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr11.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr12.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr13.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr14.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr15.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr16.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr17.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr18.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr19.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr20.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr21.txt', \
                         recomb_dir+'genetic_map_GRCh37_chr22.txt']
                            # genetic_map_GRCh37_chrX.txt
                            # genetic_map_GRCh37_chrX_par1.txt
                            # genetic_map_GRCh37_chrX_par2.txt



# create timestamp
timenow = str(datetime.datetime.now()).replace(' ', '_').replace(':', '-')
# make timestamp totally unique
timenow = timenow+'.'+format(random.randint(0,10000), '05d')



# within each parameter set, assign random values to parameters from their specified ranges 
if len(paramNumReps) == 1:
    param_range = paramNumReps[0]
    paramNumReps = []
    for param_set in range(len(output_label)):
     paramNumReps.append(param_range)        
if mean_coverage == 'pileup':
    param_sets = len(output_label)
    mean_coverage = []
    for k in range(param_sets):
        mean_coverage.append([[], []])
        for ind in [0, 1]:
            mean_coverage[k][ind] = ['pileup', 'pileup']
else:
    for param_set in range(len(output_label)):
        for ind in [0, 1]:
            m_range = mean_coverage[param_set][ind]
            mean_coverage[param_set][ind] = []
            for rep in range(paramNumReps[param_set]):
                mean_coverage[param_set][ind].append(random.uniform(m_range[0], m_range[1]))
for param_set in range(len(output_label)):
    for ind in [0, 1]:
        c_range = contam_rate[param_set][ind]
        contam_rate[param_set][ind] = []
        for rep in range(paramNumReps[param_set]):
            contam_rate[param_set][ind].append(random.uniform(c_range[0], c_range[1]))
        q_range = seq_errorrate[param_set][ind]
        seq_errorrate[param_set][ind] = []
        for rep in range(paramNumReps[param_set]):
            seq_errorrate[param_set][ind].append(random.uniform(q_range[0], q_range[1]))



# print arguments and options, with no buffer, to .args output file  
for param_set in range(len(output_label)):
    outfile_args = open(out_dir+scriptname+'.'+output_label[param_set]+'.'+timenow+'.args',      'w', 0)
    for arg in sys.argv:
        outfile_args.write(arg+'\t')
    outfile_args.write('\n')
    #outfile_args.write('dsHz = '+str(Neanderthal_Hz)+'\n')
    outfile_args.write('output_label = '+str(output_label[param_set])+'\n')
    outfile_args.write('scriptname = '+str(scriptname)+'\n')
    outfile_args.write('shuffle = '+str(shuffle)+'\n')
    outfile_args.write('downsample_rate = '+str(downsample_rate)+'\n')
    outfile_args.write('ds_rate_numSNPs = '+str(ds_rate_numSNPs)+'\n')
    outfile_args.write('data_dir = '+str(data_dir)+'\n')
    outfile_args.write('num_reps = '+str(num_reps)+'\n')
    outfile_args.write('chr_type = '+str(chr_type)+'\n')
    outfile_args.write('popAF_polymorph_filter = '+str(popAF_polymorph_filter)+'\n')
    outfile_args.write('verbose = '+str(verbose)+'\n')
    outfile_args.write('min_qual = '+str(min_qual)+'\n')
    outfile_args.write('self_comparison = '+str(self_comparison)+'\n')
    outfile_args.write('min_AF = '+str(min_AF)+'\n')
    outfile_args.write('pedigree_pop = '+str(pedigree_pop)+'\n')
    outfile_args.write('contam_pop = '+str(contam_pop)+'\n')
    outfile_args.write('contam_numind = '+str(contam_numind)+'\n')
    if recomb_dir != '':
        outfile_args.write('recomb_dir = '+str(recomb_dir)+' (using non-uniform recombination rates)\n')
    else:
        outfile_args.write('recomb_rate = '+str(recomb_rate)+'\n')
    #outfile_args.write('writeAllData = '+str(writeAllData)+'\n')
    outfile_args.write('paramNumReps = '+str(paramNumReps)+'\n')
    for rep in range(paramNumReps[param_set]):
        outfile_args.write('\tparam set '+str(param_set)+' replicate '+str(rep)+':\t'                                               \
                           'contam_rate='+str(contam_rate[param_set][0][rep])+','+str(contam_rate[param_set][1][rep])+'\t'+         \
                           'seq_errorrate='+str(seq_errorrate[param_set][0][rep])+','+str(seq_errorrate[param_set][1][rep])+'\t'+   \
                           'mean_coverage='+str(mean_coverage[param_set][0][rep])+','+str(mean_coverage[param_set][1][rep])+'\n')
    outfile_args.close()



# ==========================================================================================================================================
# HASH IN BED-LIKE TARGET SNP POSITIONS FILE (W/ ALT/REF ALLELE INFO IF INCLUDED)
# ==========================================================================================================================================
use_target_positions = 0
if targets_file != '':
    if verbose == 1:
        print 'hashing list of target positions'

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
    use_target_positions = 1



# read header from input VCF to get IDs of pedigree_pop individuals
#fileop = pysam.Tabixfile(input_chromosome_file[0],'r')
fileop = gzip.open(input_chromosome_file[0],'r')
pedigree_pop_ID_index = [None] * len(pedigree_pop_IDs)
pedigree_pop_ID_lookup = {}
contam_pop0_ID_index = [None] * len(contam_pop0_IDs)
contam_pop0_ID_lookup = {}
contam_pop1_ID_index = [None] * len(contam_pop1_IDs)
contam_pop1_ID_lookup = {}

stop = 0
while stop == 0:
    #print fileop.header
    #line = fileop.header.next()
    line = fileop.next()
    elem = line.strip('\n').split('\t')
    if elem[0].find('##') == 0:
        continue
    elif elem[0].find('#') == 0:
        # determine the index of each pedigree_pop ID
        for g in range(len(pedigree_pop_IDs)):
            pedigree_pop_ID_index[g] = elem.index(pedigree_pop_IDs[g])
            pedigree_pop_ID_lookup[pedigree_pop_ID_index[g]] = pedigree_pop_IDs[g]
        for g in range(len(contam_pop0_IDs)):
            contam_pop0_ID_index[g] = elem.index(contam_pop0_IDs[g])
            contam_pop0_ID_lookup[contam_pop0_ID_index[g]] = contam_pop0_IDs[g]
        for g in range(len(contam_pop1_IDs)):
            contam_pop1_ID_index[g] = elem.index(contam_pop1_IDs[g])
            contam_pop1_ID_lookup[contam_pop1_ID_index[g]] = contam_pop1_IDs[g]
        break
        stop = 1
    else:
        print '    error reading VCF'
        stop = 1

# define lists of IDs for choosing random pedigree and contaminating individuals
pedigree_pop_ID_choices = pedigree_pop_ID_index[:]
random.shuffle(pedigree_pop_ID_choices)
# define lists of IDs for choosing random pedigree and contaminating individuals
contam_pop0_ID_choices = contam_pop0_ID_index[:]
random.shuffle(contam_pop0_ID_choices)
contam_pop1_ID_choices = contam_pop1_ID_index[:]
random.shuffle(contam_pop1_ID_choices)

# choose random individual genomes from contam_pop[0], and remove them from pedigree_pop_ID_choices if necessary
contam0_ind_IDs = []
for c in range(contam_numind[0]):
    contam0_ind_IDs.append(contam_pop0_ID_choices.pop())
    print 'sample 0 contaminating individual #'+str(c+1)+' = '+contam_pop0_ID_lookup[contam0_ind_IDs[-1]]
    # remove the possibility that a contaminating individual genome is chosen as a pedigree individual genome
    if contam_pop0_ID_lookup[contam0_ind_IDs[-1]] in pedigree_pop_ID_choices:
        # delete this individual from the list of pedgree population IDs
        pedigree_pop_ID_choices.remove(contam_pop0_ID_lookup[contam0_ind_IDs[-1]])
# choose random individual genomes from contam_pop[1], and remove them from pedigree_pop_ID_choices if necessary
contam1_ind_IDs = []
for c in range(contam_numind[1]):
    contam1_ind_IDs.append(contam_pop1_ID_choices.pop())
    print 'sample 1 contaminating individual #'+str(c+1)+' = '+contam_pop1_ID_lookup[contam1_ind_IDs[-1]]
    # remove the possibility that a contaminating individual genome is chosen as a pedigree individual genome
    if contam_pop1_ID_lookup[contam1_ind_IDs[-1]] in pedigree_pop_ID_choices:
        # delete this individual from the list of pedgree population IDs
        pedigree_pop_ID_choices.remove(contam_pop1_ID_lookup[contam1_ind_IDs[-1]])

# pre-define unique integer "names" for up to 5000 pedigree individuals to be reproduced
newIDs = range(1000000, 1005000)
random.shuffle(newIDs)



# define comparison labels
labels = []
for k in range(0, len(comparisons)):
    labels.append(comparisons[k][0])
# labels = [  'inbred-inbred',    \
#             'father-father',    \
#             'father-child1',    \
#             'child1-child2',    \
#             'mother-cousin',    \
#             'child2-gchild',    \
#             'gchild-halfsib',   \
#             'cousin-gchild',    \
#             'father-mother'     ]            
#print 'comparisons to make:', labels



# create empty lists of results and calculations output files
outfile                 = [[] for _ in xrange(len(output_label))]
outfile_numvars         = [[] for _ in xrange(len(output_label))]
#outfile_pwdiffs         = [[] for _ in xrange(len(output_label))]
#outfile_pwdiffs_sibs    = [[] for _ in xrange(len(output_label))]
#outfile_pwdiffs_gpgc    = [[] for _ in xrange(len(output_label))]
#outfile_pwdiffs_twins   = [[] for _ in xrange(len(output_label))]
#outfile_pwdiffs_fcous   = [[] for _ in xrange(len(output_label))]
outfile_SNPscovered     = [[] for _ in xrange(len(output_label))]
outfile_labels          = [[] for _ in xrange(len(output_label))]
#if writeAllData == 1:
#    outfile_PWDdata     = [[] for _ in xrange(len(output_label))]



# perform simulations
for counter in range(1, num_reps+1):
    if shuffle == 1:
        random.shuffle(pedigree_pop_ID_choices)
    
    indiv_vars = []
    for k in individuals:
        indiv_vars.append([None] * len(input_chromosome_file))
    #father      = [None] * len(input_chromosome_file)
    #mother      = [None] * len(input_chromosome_file)
    #child1      = [None] * len(input_chromosome_file)
    #child2      = [None] * len(input_chromosome_file)
    #wife        = [None] * len(input_chromosome_file)
    #gchild      = [None] * len(input_chromosome_file)
    #cousin      = [None] * len(input_chromosome_file)
    #husband     = [None] * len(input_chromosome_file)  
    #inbred      = [None] * len(input_chromosome_file)
    #stepmom     = [None] * len(input_chromosome_file)
    #halfsib     = [None] * len(input_chromosome_file)
    
    # choose a specific, unique genome sequence for each role in pedigree (same for all chromosomes)
    print 'Choosing genomes for each individual:' 
    indiv_vars_names = []
    for k in range(0, len(indiv_vars)):
        indiv_vars_names.append(pedigree_pop_ID_choices.pop())
        print 'Individual <', individuals[k], '> source genome ID =', pedigree_pop_ID_lookup[indiv_vars_names[k]], '(', indiv_vars_names[k], ')'
    #father_name = pedigree_pop_ID_choices.pop()
    #print 'father name =', pedigree_pop_ID_lookup[father_name]
    #mother_name = pedigree_pop_ID_choices.pop()
    #print 'mother name =', pedigree_pop_ID_lookup[mother_name]    
    #wife_name = pedigree_pop_ID_choices.pop()
    #print 'wife name =', pedigree_pop_ID_lookup[wife_name]
    #husband_name = pedigree_pop_ID_choices.pop()
    #print 'husband name =', pedigree_pop_ID_lookup[husband_name]
    #stepmom_name = pedigree_pop_ID_choices.pop()
    #print 'stepmom name =', pedigree_pop_ID_lookup[stepmom_name]
    
    
    
    # define unrelated individuals
    unrels = individuals[:]
    for m in range(0, len(relationships)):
        unrels.remove(relationships[m][0])
    print 'UNRELATEDS:', unrels

    
    family_unrel_names = []
    family_unrel_chrs = []
    for chr in chrlist:     # chosen chrs in random order
        # create individual object for this chromosome for each individual
        for k in range(0, len(individuals)):
            #print 'indiv_vars[k][chr] =', indiv_vars[k][chr]
            #print 'indiv_vars_names[k] =', indiv_vars_names[k]
            indiv_vars[k][chr] = mfp.Individual(indiv_vars_names[k])
        for r in unrels: 
            family_unrel_names.append(indiv_vars_names[individuals.index(r)])
            family_unrel_chrs.append(indiv_vars[individuals.index(r)][chr])  
            
        #father[chr]     = mfp.Individual(father_name)
        #mother[chr]     = mfp.Individual(mother_name)
        #wife[chr]       = mfp.Individual(wife_name)
        #husband[chr]    = mfp.Individual(husband_name)
        #stepmom[chr]    = mfp.Individual(stepmom_name)
        #family_unrel_names = [father_name, mother_name, wife_name, husband_name, stepmom_name]
        #family_unrel_chrs = [father[chr], mother[chr], wife[chr] , husband[chr], stepmom[chr]]

        
        # print status info
        print ''
        print 'replicate: '+str(counter)+' of '+str(num_reps)+' (chr '+str(chrlist.index(chr)+1)+' of '+str(len(chrlist))+')'
        #print 'output file: '+out_dir+scriptname+'.'+output_label+'.'+timenow+'.out'
        print 'working on '+input_chromosome_file[chr]+' ...'
        

        # these are initiated here so that it can be passed to update_individuals_from1000g() even if it's empty
        pileup_positions = {}   # pileup_positions[pos] = [[basequal1, basequal2, ...], [basequal1, basequal2, ...] ]
        target_positions = {}
        

        # read base qualities from pileup file if using a pileup file for input
        if qualities_pileup_filename != '':
            pileup_positions = mfp.getQualitiesFromPileup(pileup_positions, qualities_pileup_filelist[chr], self_comparison, min_qual, verbose, chr)            


        # load region-specific recombination rates from file (if not using uniform rate)
        # bx-python was installed this way:
            # cd /emc/data/sameoldmike/software/bx-python/james_taylor-bx-python-da37e3aa45dc/
            # python setup2.py install --user
        recomb_rate_regions = []
        recomb_probs = []
        if recomb_dir != '':
            recomb_rate_regions, recomb_probs = mfp.getRecombinationRates(recomb_rates_file[chr], verbose)


        # load genoptype data from VCF
        if qualities_pileup_filename != '':
            use_pileup_positions = 1
        else:
            use_pileup_positions = 0
        print '    loading individuals data ...'
        pos_rev_lookup, numvars_actual, contam_pop_AF, pedigree_pop_AF = mfp.update_individuals_from1000g(chr,                                  \
                                                                                                          input_chromosome_file[chr],           \
                                                                                                          chrom_length[chr],                    \
                                                                                                          chrom_numvars[chr],                   \
                                                                                                          family_unrel_names,                   \
                                                                                                          family_unrel_chrs,                    \
                                                                                                          downsample_rate,                      \
                                                                                                          ds_rate_numSNPs,                      \
                                                                                                          pedigree_pop_ID_index,                \
                                                                                                          popAF_polymorph_filter,               \
                                                                                                          pileup_positions,                     \
                                                                                                          use_pileup_positions,                 \
                                                                                                          target_positions,                     \
                                                                                                          use_target_positions,                 \
                                                                                                          pedigree_pop,                         \
                                                                                                          contam_pop,                           \
                                                                                                          [contam0_ind_IDs, contam1_ind_IDs],   \
                                                                                                          min_AF)
        print '        done.'
        #print 'contam_pop_AF =', contam_pop_AF
        
        # make chromosome_map (map of variant chromosome positions) from pos_rev_lookup
        chromosome_map = sorted(pos_rev_lookup)
                                
        # perform pedigree reproductions 
        print '    performing reproductions ...'
        for r in range(0, len(relationships)):
            indiv_vars[individuals.index(relationships[r][0])][chr] =  mfp.repro(indiv_vars[individuals.index(relationships[r][1])][chr], indiv_vars[individuals.index(relationships[r][2])][chr], recomb_rate, chromosome_map, chrom_length[chr], newIDs.pop(), recomb_rate_regions, recomb_probs)            

        #child1[chr]     = mfp.repro(father[chr], mother[chr],   recomb_rate, chromosome_map, chrom_length[chr], newIDs.pop(), recomb_rate_regions, recomb_probs)
        #child2[chr]     = mfp.repro(father[chr], mother[chr],   recomb_rate, chromosome_map, chrom_length[chr], newIDs.pop(), recomb_rate_regions, recomb_probs)
        #gchild[chr]     = mfp.repro(child1[chr], wife[chr],     recomb_rate, chromosome_map, chrom_length[chr], newIDs.pop(), recomb_rate_regions, recomb_probs)
        #cousin[chr]     = mfp.repro(child2[chr], husband[chr],  recomb_rate, chromosome_map, chrom_length[chr], newIDs.pop(), recomb_rate_regions, recomb_probs)
        #inbred[chr]     = mfp.repro(child1[chr], child2[chr],   recomb_rate, chromosome_map, chrom_length[chr], newIDs.pop(), recomb_rate_regions, recomb_probs) # offspring of siblings
        #halfsib[chr]    = mfp.repro(child1[chr], stepmom[chr],  recomb_rate, chromosome_map, chrom_length[chr], newIDs.pop(), recomb_rate_regions, recomb_probs)

        
        # calculate/write pairwise differences in order of increasing genetic distance
        print '    calculating pairwise differences ...'    
        # [fraction, numpwdiffs, numSNPscovered, totSNPs, pwdiff_expected_considering_coverage] = pwdiffsresults
        # NOTE: pwdiff_expected [.pwexp] is NOT normalized by length of chromosome
        # NOTE: pw_sum [.out]            is normalized by length of chromosome

        for param_set in range(len(output_label)):       
            for rep in range(paramNumReps[param_set]): 
                
                # open results and calculations output files, without buffer, for write/append
                openOutputFiles()

                # get estimates of pairwise distance
                for c in range(0, len(comparisons)):                
                    if qualities_pileup_filename == '':
                        pwdiffsresults = mfp.pairwise_diff(indiv_vars[individuals.index(comparisons[c][1])][chr], indiv_vars[individuals.index(comparisons[c][2])][chr], mean_coverage[param_set][0][rep], mean_coverage[param_set][1][rep], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, chromosome_map)
                    else:
                        pwdiffsresults = mfp.pileup_PWD(indiv_vars[individuals.index(comparisons[c][1])][chr], indiv_vars[individuals.index(comparisons[c][2])][chr], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, pos_rev_lookup, pileup_positions, verbose)
                    pw_sum = pwdiffsresults[1]
                    writePWDtoFiles(labels.index(comparisons[c][0]))
                    #writePWDtoFiles(labels.index(comparisons[c][0].strip('"')))
                
                """
                # IDENTICAL RELATIONSHIP [self-comparison of an offspring of siblings]
                if qualities_pileup_filename == '':
                    pwdiffsresults = mfp.pairwise_diff(inbred[chr], inbred[chr], mean_coverage[param_set][0][rep], mean_coverage[param_set][1][rep], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, chromosome_map)
                else:
                    pwdiffsresults = mfp.pileup_PWD(inbred[chr], inbred[chr], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, pos_rev_lookup, pileup_positions, verbose)
                pw_sum = pwdiffsresults[1]
                writePWDtoFiles(labels.index('inbred-inbred'))

                # IDENTICAL RELATIONSHIP
                if qualities_pileup_filename == '':
                    pwdiffsresults = mfp.pairwise_diff(father[chr], father[chr], mean_coverage[param_set][0][rep], mean_coverage[param_set][1][rep], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, chromosome_map)
                else:
                    pwdiffsresults = mfp.pileup_PWD(father[chr], father[chr], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, pos_rev_lookup, pileup_positions, verbose)
                pw_sum = pwdiffsresults[1]
                writePWDtoFiles(labels.index('father-father'))

                # PARENT-CHILD RELATIONSHIP
                if qualities_pileup_filename == '':
                    pwdiffsresults = mfp.pairwise_diff(father[chr], child1[chr], mean_coverage[param_set][0][rep], mean_coverage[param_set][1][rep], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, chromosome_map)
                else:
                    pwdiffsresults = mfp.pileup_PWD(father[chr], child1[chr], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, pos_rev_lookup, pileup_positions, verbose)
                pw_sum = pwdiffsresults[1]
                writePWDtoFiles(labels.index('father-child1'))    
        
                # SIBLING RELATIONSHIP
                if qualities_pileup_filename == '':
                    pwdiffsresults = mfp.pairwise_diff(child1[chr], child2[chr], mean_coverage[param_set][0][rep], mean_coverage[param_set][1][rep], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, chromosome_map)
                else:
                    pwdiffsresults = mfp.pileup_PWD(child1[chr], child2[chr], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, pos_rev_lookup, pileup_positions, verbose)
                pw_sum = pwdiffsresults[1]
                writePWDtoFiles(labels.index('child1-child2'))

                # GRANDPARENT-GRANDCHILD RELATIONSHIP
                if qualities_pileup_filename == '':
                    pwdiffsresults = mfp.pairwise_diff(mother[chr], cousin[chr], mean_coverage[param_set][0][rep], mean_coverage[param_set][1][rep], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, chromosome_map)
                else:
                    pwdiffsresults = mfp.pileup_PWD(mother[chr], cousin[chr], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, pos_rev_lookup, pileup_positions, verbose)
                pw_sum = pwdiffsresults[1]
                writePWDtoFiles(labels.index('mother-cousin')) 

                # AVUNCULAR RELATIONSHIP
                if qualities_pileup_filename == '':
                    pwdiffsresults = mfp.pairwise_diff(child2[chr], gchild[chr], mean_coverage[param_set][0][rep], mean_coverage[param_set][1][rep], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, chromosome_map)
                else:
                    pwdiffsresults = mfp.pileup_PWD(child2[chr], gchild[chr], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, pos_rev_lookup, pileup_positions, verbose)
                pw_sum = pwdiffsresults[1]
                writePWDtoFiles(labels.index('child2-gchild'))

                # HALF-SIBLING RELATIONSHIP
                if qualities_pileup_filename == '':
                    pwdiffsresults = mfp.pairwise_diff(gchild[chr], halfsib[chr], mean_coverage[param_set][0][rep], mean_coverage[param_set][1][rep], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, chromosome_map)
                else:
                    pwdiffsresults = mfp.pileup_PWD(gchild[chr], halfsib[chr], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, pos_rev_lookup, pileup_positions, verbose)
                pw_sum = pwdiffsresults[1]
                writePWDtoFiles(labels.index('gchild-halfsib'))        

                # COUSIN RELATIONSHIP
                if qualities_pileup_filename == '':
                    pwdiffsresults = mfp.pairwise_diff(cousin[chr], gchild[chr], mean_coverage[param_set][0][rep], mean_coverage[param_set][1][rep], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, chromosome_map)
                else:
                    pwdiffsresults = mfp.pileup_PWD(cousin[chr], gchild[chr], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, pos_rev_lookup, pileup_positions, verbose)
                pw_sum = pwdiffsresults[1]
                writePWDtoFiles(labels.index('cousin-gchild')) 

                # UNRELATED RELATIONSHIP
                if qualities_pileup_filename == '':
                    pwdiffsresults = mfp.pairwise_diff(father[chr], mother[chr], mean_coverage[param_set][0][rep], mean_coverage[param_set][1][rep], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, chromosome_map)
                else:
                    pwdiffsresults = mfp.pileup_PWD(father[chr], mother[chr], contam_rate[param_set][0][rep], contam_rate[param_set][1][rep], seq_errorrate[param_set][0][rep], seq_errorrate[param_set][1][rep], contam_pop_AF, pedigree_pop_AF, pos_rev_lookup, pileup_positions, verbose)
                pw_sum = pwdiffsresults[1]
                writePWDtoFiles(labels.index('father-mother'))        
                """

                # write labels to .labels file 
                for t in range(len(labels)):
                    outfile_labels[param_set][rep].write(labels[t]+'\t')


                # write an endline to all output files
                writeOutputFileEndlines()

                # close output files
                closeOutputFiles()

            # return each list of open files to an empty list
            resetFileLists()


print 'simulations finished.'




