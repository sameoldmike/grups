# grups v1.2 DOCUMENTATION


## 1. About the software 

grups is a collection of Python and R scripts that implement methods developed 
in Martin et al. (2017). For citation information, see section 2. The scripts 
can be used to simulate the pairwise genetic distances between human 
individuals in a provided test pedigree based on either existing genomic 
genotype data (VCF file input) or genomic sequence data (BAM file input). 
Contamination and sequencing error parameters can also be introduced to the 
simulations. Simulations can be performed using only a subset of the input data 
(genomic positions passing desired filters, or included on an input list). A 
massively parallel computing cluster is recommended for performing large 
numbers of pedigree replicates using entire human genomes.


## 2. Citing grups

Martin MD, Jay F, Castellano S, Slatkin M. 2017. Determination of genetic relatedness from low‐coverage human genome sequences using pedigree simulations. Molecular Ecology, 26:4145-4157. doi:10.1111/mec.14188 


## 3. Prerequisites

grups v1.1 was developed and tested in the following environment:

- Python v2.7.5
- Python package numpy v1.11.2
- Python package matplotlib v1.5.3
- Python package pysam v0.8.4
- Python package cython v0.24.1
- Python package bx-python v0.5.0
- R v3.3.1
- R library reshape v0.8.5
- R library ggplot2 v2.1.0
- samtools v0.1.19


## 4. Usage

Generally the first step in an analysis of pairwise genetic distances will be 
to use the script `PWD_from_stdin.py` to directly calculate the pairwise 
genetic distance between the sequence data for two individuals in two BAM files 
(e.g. command **5a1**). Then summary statistics are calculated from the output 
files (e.g. command **5a2**). Next pileup files should be generated (e.g. command 
**5b1**) for each pairwise comparison to be simulated. Then each pairwise pileup 
file must be split into chromosome-specific pileup files (e.g. command **5b2**), 
which are used as input for the pedigree simulations performed by the script 
`pedigree_sims.py` (e.g. command **5c**). In general the pedigree simulations 
command should be performed in hundreds of replicates in order to accurately 
estimate the mean and variance of the relatedness distribution for each 
relationship. Finally, the simulation replicates are concatenated into a single 
input file (e.g. command **5d1**), which is the input for an R script that plots 
the results of the pedigree simulations (command **5d2**) and, optionally, the 
previously calculated summary statistics from direct observations of pairwise 
genetic distance between two BAMs.

Alternatively, more general pedigree simulations of various scenarios can be 
performed without using pileup files as input,both with (`--targets` option) and 
without (e.g. command **5e**) using an input file of target genomic positions 
defining where genetic distance should be assessed. A target regions input file 
can be provided to restrict the calculation of genetic distance to particular 
genomic locations defined in the file. If a pairwise pileup file has also been 
provided, the target positions file may contain a subset of the pileup 
positions, but a position must be in both lists to be included. 

Major base data used in the Martin et al. (2017) study can be found as follows:

* The 1000 Genomes v37 reference genome: [ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz/](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz/)

* HapMap recombination data: [ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz](ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz)


## 5. Usage example commands

The scripts can be used in several modes. Following are some examples of typical usage:

#### 5a. Calculate pairwise genetic distance between two BAM files

##### 5a1

```bash
for y in {1..20}; do
    samtools mpileup -q 30 -Q 30 \
    -f ./1000genomes_v3_refseq/human_g1k_v37.fasta \
    -l ./1000genomes_phase3/variant_sites.chr1-22.only_SNPs.txt \
    ./1000g_v3_BAMs/NA20526.mapped.ILLUMINA.bwa.TSI.low_coverage.20130415.rmdup.bam \
    ./1000g_v3_BAMs/NA20792.mapped.ILLUMINA.bwa.TSI.low_coverage.20120522.rmdup.bam | \
    python PWD_from_stdin.py --chr 1-22 --min_depth 2,2 \
        --min_qual 30 --ignore_dels 1 \ 
        --quiet 1 >> ./analyses/PWD_distributions/PWD_from_stdin.py.only_targets.only_SNPs.mindepth2_2.NA20526.rmdup_NA20792.rmdup.out
done &
```

##### 5a2.

```bash
cat ./analyses/PWD_distributions/PWD_from_stdin_2015-06-05.py.only_targets.only_SNPs.mindepth2_2.NA20526.rmdup_NA20792.rmdup.out | \
    awk 'NF > 0' | \
    awk 'BEGIN{
           max=0;
           min=10000;
         }{
           sum+=$1;
           sumsq+=$1*$1;
           if ($1>max){max=$1};
           if ($1<min){min=$1}
         }END{
           print "n mean std max min";
           print NR " " sum/NR " " sqrt(sumsq/NR - (sum/NR)**2) " " max " " min
         }'
```

#### 5b. Generate and split pileup files for the simulation of pairwise genetic distances within a pedigree of randomly chosen individuals, based on input pairwise pileup files

##### 5b1.

```bash
samtools mpileup -q 30 -Q 30 -f ./1000genomes_v3_refseq/human_g1k_v37.fasta \
    -l ./1000genomes_phase3/variant_sites.chr1-22.only_SNPs.txt \
    /path/to/1000g_data/NA20526.mapped.ILLUMINA.bwa.TSI.low_coverage.20130415.rmdup.bam \
    /path/to/1000g_data/NA20792.mapped.ILLUMINA.bwa.TSI.low_coverage.20120522.rmdup.bam | \
    python PWD_from_stdin.py --chr 1-22 --min_depth 2,2 \
        --min_qual 30 --ignore_dels 1 -- quiet 1 --filter_sites 1 | \
    gzip > ./1000g_v3_BAMs/pairwise_pileups_simspecific/NA20526.rmdup_NA20792.rmdup.PWD_from_stdin_2015-06-05.py.only_targets.only_SNPs.min_depth2_2.pileup.gz &
```

or

```bash
samtools mpileup -q 30 -Q 30 -f ./1000genomes_v3_refseq/human_g1k_v37.fasta \
    -l ./Haaketal2015_data/SuppTable2_targetSNPs.chr_pos.txt \
    ./Haaketal2015_data/I0114.390k.damage_rescaled.bam \
    ./Haaketal2015_data/I0114.390k.damage_rescaled.bam | \
    python PWD_from_stdin.py --chr 1-22 --min_depth 2,2 \
        --min_qual 30 --known_variants 1 \
        --targets ./1000genomes_phase3_v5a/variant_sites.chr1-22.only_SNPs.only_tv.txt \
        --ignore_dels 1 --self_comparison 1 --quiet 1 --filter_sites 1 | \
        gzip > ./Haaketal2015_data/pairwise_pileups_simspecific/ESP2_MDR_ESP2_MDR.PWD_from_stdin_2015-06-05.py.only_targets2.only_SNPs.only_transversions.min_depth2_2.pileup.gz
```


##### 5b2.

```bash
for y in {1..22}; do 
    zcat ./1000g_v3_BAMs/pairwise_pileups_simspecific/NA20526.rmdup_NA20792.rmdup.PWD_from_stdin_2015-06-05.py.only_targets.only_SNPs.min_depth2_2.pileup.gz | \
    awk -v y="$y" '{if($1 == y) print $0}' | \
    gzip > ./1000g_v3_BAMs/pairwise_pileups_simspecific/NA20526.rmdup_NA20792.rmdup.PWD_from_stdin_2015-06-05.py.only_targets.only_SNPs.min_depth2_2.pileup.chr${y}.gz & 
done
```

#### 5c. Simulate pairwise genetic distances within a pedigree based on positions and read depths in an input pairwise pileup file

```bash
python pedigree_sims.py 
    --label NA20526rmdup_NA20792rmdup_pileup 
    --ped example_pedigree.txt
    --pileup ../1000g_v3_BAMs/pairwise_pileups_simspecific/NA20526.rmdup_NA20792.rmdup.PWD_from_stdin_2015-06-05.py.only_targets.only_SNPs.min_depth2_2.pileup.gz 
    --out /path/to/simulations_output/ 
    --ds_rate 1.0 
    --c_rate 0.0,0.0
    --q_rate 0.0,0.0 
    --reps 1
    --chr 1-22 
    --include all 
    --min_qual 30 
    --data_dir /path/to/1000g_data/
    --pedigree_pop EUR 
    --contam_pop AFR,AFR 
    --recomb_dir ./recombination_maps/
```

#### 5d. Concatenate and plot simulation statistics

##### 5d1.

```bash
bash pedigree_sims_concat.sh Troc1rmdup_Troc3rmdup_pileup \
    ../analyses/pedigree_sims/2016-11-11_Haaketal_mindepth1_onlyTv_ds1.0_c0_q0/ \
    ConcatenatedReps/
```

##### 5d2.

```bash
Rscript plot_pedigree_sims.R \
    data_dir=../analyses/pedigree_sims/2016-11-11_Haaketal_mindepth1_onlyTv_ds1.0_c0_q0/ConcatenatedReps/ \
    regex=*Troc1rmdup_Troc3rmdup_pileup.out \
    range=1.5 plotval=0.225893,0.00122087 w=6 h=6 noPrint=1,3,6,7
```
**PLOTTING STILL NEEDS TO BE DEFINED BY LABELS FILE. NOT COMPLETELY IMPLEMENTED.**
    
#### 5e. Simulate pairwise genetic distances within a pedigree of randomly chosen individuals, with genetic distance calculations at all genomic positions contained in the pileup file that meet the read depth/quality criteria.

```bash
python pedigree_sims.py 
    --label=Scenario1_q0.0000/Scenario2_q0.0005/Scenario3_q0.0010/Scenario4_q0.0100/Scenario5_q0.1000 
    --ped example_pedigree.txt
    --ds_rate_numSNPs 0.047062308 
    --min_AF 0.05 
    --chr 1-22 
    --include all 
    --reps 1
    --out /path/to/simulations_output/ 
    --pedigree_pop EUR 
    --contam_pop EUR,EUR 
    --data_dir /path/to/1000g_data/
    --recomb_dir ./recombination_maps/ 
    --q_rate 0.0000,0.0000/0.0005,0.0005/0.0010,0.0010/0.0100,0.0100/0.1000,0.1000 
    --c_rate 0.0,0.0/0.0,0.0/0.0,0.0/0.0,0.0/0.0,0.0  
    --mean_cov 10,10/10,10/10,10/10,10/10,10         
```

## 6. Additional resources

For more information check: [https://github.com/sameoldmike/grups/](https://github.com/sameoldmike/grups/).

To report issues, contact Mike Martin: mike.martin@ntnu.no.  When reporting 
issues, please be as specific as possible regarding software environment, 
version numbers, and output. 

