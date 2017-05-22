=== grups v1.0 DOCUMENTATION ===



1. About the software. grups is a collection of Python and R scripts that implement methods developed in Martin et al. (2017). For citation information, see section XC. The scripts can be used to simulate the pairwise genetic distances between human individuals in a provided test pedigree based on either existing genomic genotype data (VCF file input) or genomic sequence data (BAM file input). Contamination and sequencing error parameters can also be introduced to the simulations. Simulations can be performed using only a subset of the input data (genomic positions passing desired filters, or included on an input list).  


2. Prerequisites. grups v1.0 was developed and tested in the following environment:
Python v2.7.5
Python package numpy v1.11.2
Python package matplotlib v1.5.3
Python package pysam v0.8.4
Python package cython v0.24.1
Python package bx-python v0.5.0
R v3.3.1
R library reshape v0.8.5
R library ggplot2 v2.1.0
samtools v0.1.19


3. Usage. The scripts can be used in several modes. Following are some examples of typical usage:
A) Simulate pairwise genetic distances within a pedigree based on input VCFs
    python pedigree_sims.py 
        label=DemonstrateSeqErr_SNPs300k_cov10X_minAF0.05_q0.0000/DemonstrateSeqErr_SNPs300k_cov10X_minAF0.05_q0.0005/DemonstrateSeqErr_SNPs300k_cov10X_minAF0.05_q0.0010/DemonstrateSeqErr_SNPs300k_cov10X_minAF0.05_q0.0100/DemonstrateSeqErr_SNPs300k_cov10X_minAF0.05_q0.1000 ds_rate_numSNPs=0.047062308 
        min_AF=0.05 
        chr_type=whole 
        chr=1-22 
        include=all 
        reps=1
        out=/home/people/michmar/humans/analyses/pedigree_sims/2016-10-14_demo_seq_err/ 
        pedigree_pop=EUR 
        contam_pop=EUR,EUR 
        data_dir=/home/people/michmar/humans/1000genomes_v3_data_EURonly/ 
        recomb_dir=/home/people/michmar/humans/analyses/recombination_maps/ 
        q_rate=0.0000,0.0000/0.0005,0.0005/0.0010,0.0010/0.0100,0.0100/0.1000,0.1000 
        c_rate=0.0000,0.0000/0.0000,0.0000/0.0000,0.0000/0.0000,0.0000 
        mean_cov=10,10/10,10/10,10/10,10/10,10 

B) Generate and split pileup files for the simulation of pairwise genetic distances within a pedigree based on input BAMs
    1. samtools mpileup -q 30 -Q 30 -f ./1000genomes_v3_refseq/human_g1k_v37.fasta -l ./1000genomes_phase3/variant_sites.chr1-22.only_SNPs.txt ./1000g_v3_BAMs/NA20526.mapped.ILLUMINA.bwa.TSI.low_coverage.20130415.rmdup.bam ./1000g_v3_BAMs/NA20792.mapped.ILLUMINA.bwa.TSI.low_coverage.20120522.rmdup.bam | nice -n19 python ./scripts/PWD_from_stdin_2015-06-05.py chr=1-22 min_depth=2,2 min_qual=30 ignore_dels quiet filter_sites | gzip > ./1000g_v3_BAMs/pairwise_pileups_simspecific/NA20526.rmdup_NA20792.rmdup.PWD_from_stdin_2015-06-05.py.only_targets.only_SNPs.min_depth2_2.pileup.gz &
    2. for y in {1..22}; do zcat ./1000g_v3_BAMs/pairwise_pileups_simspecific/NA20526.rmdup_NA20792.rmdup.PWD_from_stdin_2015-06-05.py.only_targets.only_SNPs.min_depth2_2.pileup.gz | awk -v y="$y" '{if($1 == y) print $0}' | gzip > ./1000g_v3_BAMs/pairwise_pileups_simspecific/NA20526.rmdup_NA20792.rmdup.PWD_from_stdin_2015-06-05.py.only_targets.only_SNPs.min_depth2_2.pileup.chr${y}.gz & done

C) Simulate pairwise genetic distances within a pedigree based on input BAMs
    python pedigree_sims.py label=NA20526rmdup_NA20792rmdup_pileup pileup=../1000g_v3_BAMs/pairwise_pileups_simspecific/NA20526.rmdup_NA20792.rmdup.PWD_from_stdin_2015-06-05.py.only_targets.only_SNPs.min_depth2_2.pileup.gz out=../analyses/pedigree_sims/2016-11-11_1000g_PileupSims_q0_c0/ ds_rate=1.0 c_rate=0.0,0.0 q_rate=0.0,0.0 reps=1 chr_type=whole chr=1-22 include=all min_qual=30 data_dir=../1000genomes_v3_data_EURonly/ pedigree_pop=EUR contam_pop=EUR,EUR recomb_dir=../analyses/recombination_maps/
    
D) Concatenate and plot simulation statistics
    bash pedigree_sims_concat.sh Troc1rmdup_Troc3rmdup_pileup ../analyses/pedigree_sims/2016-11-11_Haaketal_mindepth1_onlyTv_ds1.0_c0_q0/ ConcatenatedReps/
    Rscript plot_pedigree_sims.R data_dir=../analyses/pedigree_sims/2016-11-11_Haaketal_mindepth1_onlyTv_ds1.0_c0_q0/ConcatenatedReps/ regex=*Troc1rmdup_Troc3rmdup_pileup.out range=1.5 violin plotval=0.225893,0.00122087 w=6 h=6 noPrint=1,3,6,7

E) Calculate pairwise genetic distance between two BAM files (pileup input)
    1. for y in {1..20}; do nice -n19 samtools mpileup -q 30 -Q 30 -f ./1000genomes_v3_refseq/human_g1k_v37.fasta -l ./1000genomes_phase3/variant_sites.chr1-22.only_SNPs.txt ./1000g_v3_BAMs/NA20526.mapped.ILLUMINA.bwa.TSI.low_coverage.20130415.rmdup.bam ./1000g_v3_BAMs/NA20792.mapped.ILLUMINA.bwa.TSI.low_coverage.20120522.rmdup.bam | nice -n19 python ./scripts/PWD_from_stdin.py chr=1-22 min_depth=2,2 min_qual=30 ignore_dels quiet >> ./analyses/PWD_distributions/PWD_from_stdin.py.only_targets.only_SNPs.mindepth2_2.NA20526.rmdup_NA20792.rmdup.out ; done &
    2. cat ./analyses/PWD_distributions/PWD_from_stdin_2015-06-05.py.only_targets.only_SNPs.mindepth2_2.NA20526.rmdup_NA20792.rmdup.out | awk 'NF > 0' | awk 'BEGIN{max=0; min=10000;} {sum+=$1; sumsq+=$1*$1; if ($1>max){max=$1}; if ($1<min){min=$1}}END {print "n mean std max min"; print NR " " sum/NR " " sqrt(sumsq/NR - (sum/NR)**2) " " max " " min}'




4. Citing grups: Martin MD, Jay F, Castellano S, Slatkin M. 2017. Molecular Ecology.


5. For more information or to report issues, contact Mike Martin: mike.martin@ntnu.no

