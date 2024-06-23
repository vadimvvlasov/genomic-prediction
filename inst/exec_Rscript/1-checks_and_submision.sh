#!/bin/bash

### Full path to the location of the executable Rscript `gp.R`` which should co-locate with this script: `0-checks_and_submission.sh`` as well as `1-gp_slurm_job.sh`.
DIR_SRC=$(dirname $0)
cd $DIR_SRC
DIR_SRC=$(pwd)
##################################
### Load the conda environment ###
##################################
echo "Loading the conda environment..."
module load Miniconda3/22.11.1-1
conda init bash
source ~/.bashrc
if [ $(conda env list | grep "^genomic_selection " | wc -l) -gt 0 ]
then
    conda activate genomic_selection
else
    echo "Installing the conda environment..."
    conda env create -f ${DIR_SRC}/../../conda.yml
    conda activate genomic_selection
fi
#######################################
### Install gp if not installed yet ###
#######################################
echo "Checking if the R library gp exists and installing it, if not..."
Rscript -e 'if (!require("gp", character.only = TRUE)) {install.packages("devtools", repos="https://cloud.r-project.org"); devtools::install_github("jeffersonfparil/gp")}'
################################################################
### TOP-LEVEL SLURM ARRAY JOB SUBMISSION SCRIPT
### Please edit the input variables below to match your dataset:
################################################################
### Input variables (use the absolute path to files to be precise)
### (1) R matrix object with n rows corresponding to samples, and p columns corresponding to the markers or loci. 
###     - The genotype data can be coded as any numeric range of values, e.g. (0,1,2), (-1,0,1), and (0.00,0.25,0.50,0.75,1.00) or as biallelic characters, e.g. for diploids: "AA", "AB", "BB", and for tetraploids: "AAAA", "AAAB", "AABB", "ABBB", and "BBBB".. It is recommended that this data should be filtered and imputed beforehand.
###     - The rows are expected to have names of the samples corresponding to the names in the phenotype file.
###     - The columns are expected to contain the loci names but does need to follow a specific format: chromosome name and position separated by a tab character (`\t`) and an optional allele identifier, e.g. `chr-1\t12345\tallele_A`
GENOTYPE_DATA_RDS=${DIR_SRC}/input/test_geno.Rds
### (2) Tab-delimited phenotype file where column 1: sample names, column 2: population name, columns 3 and so on refers to the phenotype values of one trait per column.
###     - Headers for the columns should be named appropriately, e.g. ID, POP, TRAIT1, TRAIT2, etc.
###     - Missing values are allowed for samples whose phenotypes will be predicted by the best model identified within the population they belong to.
###     - Missing values may be coded as empty cells, -, NA, na, NaN, missing, and/or MISSING.
PHENOTYPE_DATA_TSV=${DIR_SRC}/input/test_pheno.tsv
### (3) Number of folds for k-fold cross-validation.
KFOLDS=5
### (4) Number of replications of the k-fold cross-validation each representing a random sorting of the samples hence yielding different ways of partitioning the data.
NREPS=3
### (5) Output directory where the output/ folder will be created
DIR_OUT=${DIR_SRC}

### Check if the genotype file exists
if [ ! -f $GENOTYPE_DATA_RDS ]
then
    echo "Error: The genotype file: $GENOTYPE_DATA_RDS does not exist. Are you specifying the full path? Is the name correct?"
    rm $0
    rm 2-gp_slurm_job.sh
    exit 101
else
    echo "Passed: The genotype file: $GENOTYPE_DATA_RDS exists."
fi
### Check if the phenotype file exists
if [ ! -f $PHENOTYPE_DATA_TSV ]
then
    echo "Error: The phenotype file: $PHENOTYPE_DATA_TSV does not exist. Are you specifying the full path? Is the name correct?"
    rm $0
    rm 2-gp_slurm_job.sh
    exit 102
else
    echo "Passed: The phenotype file: $PHENOTYPE_DATA_TSV exists."
fi
### Check if the genotype file is a valid Rds file
echo 'args = commandArgs(trailingOnly=TRUE)
geno = suppressWarnings(tryCatch(readRDS(args[1]), error=function(e){print("Error loading genotype file.")}))
' > test_geno_rds.R
if [ $(Rscript test_geno_rds.R $GENOTYPE_DATA_RDS | grep -i "error" | wc -l) -eq 1 ]
then
    echo "Error: The genotype file: $GENOTYPE_DATA_RDS is not an Rds file."
    rm $0
    rm 2-gp_slurm_job.sh
    exit 103
else
    echo "Passed: The genotype file: $GENOTYPE_DATA_RDS is an Rds file."
fi
rm test_geno_rds.R
### Check if the phenotype file is formatted according to the required specifications
echo 'args = commandArgs(trailingOnly=TRUE)
pheno = suppressWarnings(tryCatch(read.delim(args[1], sep="\t", header=TRUE), error=function(e){print("Error loading phenotype file.")}))
' > test_pheno_rds.R
if [ $(Rscript test_pheno_rds.R $PHENOTYPE_DATA_TSV | grep -i "error" | wc -l) -eq 1 ]
then
    echo "Error: The phenotype file: $GENOTYPE_DATA_RDS is not formatted according to specifications. It should be tab-delimited and a header line must be present."
    rm $0
    rm 2-gp_slurm_job.sh
    exit 104
else
    echo "Passed: The phenotype file: $GENOTYPE_DATA_RDS is tab-delimited with a header line."
fi
rm test_pheno_rds.R
### Check if the output directory belongs to the user
if [ ! -w $DIR_OUT ]
then
    echo "Error: You do not have permission to write in the output directory: $DIR_OUT. Please use an output directory you have write acess to."
    rm $0
    rm 2-gp_slurm_job.sh
    exit 106
else
    echo "Passed: You have permission to write in the output directory: $DIR_OUT."
fi
### Check if the configured slurm array job script exists
if [ ! -f ${DIR_SRC}/2-gp_slurm_job.sh ]
then
    echo "Error: The executable code directory: $DIR_SRC does not contain the script: 2-gp_slurm_job.sh. Are you sure this is the genomic_selection repo directory?"
    rm $0
    rm 2-gp_slurm_job.sh
    exit 107
else
    echo "Passed: The executable code directory: $DIR_SRC contains the script: 2-gp_slurm_job.sh."
fi
### Check if the executable Rscript exists
if [ ! -f ${DIR_SRC}/gp.R ]
then
    echo "Error: The executable code directory: $DIR_SRC does not contain the script: gp.R. Are you sure this is the genomic_selection repo directory?"
    rm $0
    rm 2-gp_slurm_job.sh
    exit 108
else
    echo "Passed: The executable code directory: $DIR_SRC contains the Rscript: gp.R."
fi
### Initialise the output directory which will contain all the output Rds files across populations and traits
if [ ! -d ${DIR_OUT}/output ]
then
    mkdir ${DIR_OUT}/output
fi
### Submit an array of jobs equivalent to the number of traits in the phenotype file
cd $DIR_SRC/
N_TRAITS=$(echo $(head -n1 $PHENOTYPE_DATA_TSV | awk '{print NF}') - 2 | bc)
N_POPS=$(tail -n+2 $PHENOTYPE_DATA_TSV | cut -f2 - | sort | uniq | wc -l)
sbatch --array 1-$(echo "${N_TRAITS} * ${N_POPS}" | bc) \
    2-gp_slurm_job.sh \
        ${GENOTYPE_DATA_RDS} \
        ${PHENOTYPE_DATA_TSV} \
        ${KFOLDS} \
        ${NREPS} \
        ${DIR_SRC} \
        ${DIR_OUT}