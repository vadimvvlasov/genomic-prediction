#!/bin/bash

### Load the conda environment
module load Miniconda3/22.11.1-1
conda init bash
source ~/.bashrc
conda activate genomic_selection

################################################################
### TOP-LEVEL SLURM ARRAY JOB SUBMISSION SCRIPT
### Please edit the input variables below to match your dataset:
################################################################

### Input variables (use the absolute path to files to be precise)
### (1) R matrix object with n rows corresponding to samples, and p columns corresponding to the markers or loci. 
###     Should have no missing data or else will be imputed via mean value imputation.
GENOTYPE_DATA_RDS='/group/pasture/Jeff/lucerne/workdir/FINAL-IMPUTED-noTrailingAllele-filteredSNPlist.Rds'
### (2) Tab-delimited phenotype file where column 1: sample names, column 2: population name, columns 3 and so on refers to the phenotype values of one trait per column.
###     Headers for the columns should be named appropriately, e.g. ID, POP, TRAIT1, TRAIT2, etc.
###     Missing values are allowed for samples whose phenotypes will be predicted by the best model identified within the population they belong to.
###     Missing values may be coded as empty cells, -, NA, na, NaN, missing, and/or MISSING.
PHENOTYPE_DATA_TSV='/group/pasture/Jeff/lucerne/workdir/Lucerne_PhenomicsDB_2024-05-27-BiomassPredicted.tsv'
### (3) Number of folds for k-fold cross-validation.
KFOLDS=5
### (4) Number of replications of the k-fold cross-validation each representing a random sorting of the samples hence yielding different ways of partitioning the data.
NREPS=5
### (5) Full path to the location of the executable Rscript gp.R
DIR='/group/pasture/Jeff/gp/exec'

### Check if the genotype file exists
if [ ! -f $GENOTYPE_DATA_RDS ]
then
    echo "Error: The genotype file: $GENOTYPE_DATA_RDS does not exist. Are you specifying the full path? Is the name correct?"
    exit 101
fi
### Check if the phenotype file exists
if [ ! -f $PHENOTYPE_DATA_TSV ]
then
    echo "Error: The phenotype file: $PHENOTYPE_DATA_TSV does not exist. Are you specifying the full path? Is the name correct?"
    exit 102
fi
### Check if the genotype file is a valid Rds file
echo 'args = commandArgs(trailingOnly=TRUE)
geno = suppressWarnings(tryCatch(readRDS(args[1]), error=function(e){print("Error loading genotype file.")}))
' > test_geno_rds.R
if [ $(Rscript test_geno_rds.R $GENOTYPE_DATA_RDS | grep -i "error" | wc -l) -eq 1 ]
then
    echo "Error: The genotype file: $GENOTYPE_DATA_RDS is not an Rds file."
    exit 103
fi
rm test_geno_rds.R
### Check if the phenotype file is formatted according to the required specifications
echo 'args = commandArgs(trailingOnly=TRUE)
pheno = suppressWarnings(tryCatch(read.delim(args[1], sep="\t", header=TRUE), error=function(e){print("Error loading phenotype file.")}))
' > test_pheno_rds.R
if [ $(Rscript test_pheno_rds.R $PHENOTYPE_DATA_TSV | grep -i "error" | wc -l) -eq 1 ]
then
    echo "Error: The phenotype file: $GENOTYPE_DATA_RDS is not formatted according to specifications. It should be tab-delimited and a header line must be present."
    exit 104
fi
rm test_pheno_rds.R
### Check if the genomic_selection repo folder exists
if [ ! -d $DIR ]
then
    echo "Error: The genotype_selection directory: $DIR does not exist. Are you specifying the full path? Is the name correct?"
    exit 105
fi
### Check if the genomic_selection repo belongs to the user
if [ ! -w $DIR ]
then
    echo "Error: You do not have permission to write in the genotype_selection directory: $DIR. Did you clone the genomic_selection repository into a directory you have write-access to?"
    exit 106
fi
### Check if the genomic_selection repo has contains the slurm array job submission script
if [ ! -f ${DIR}/01_gp_slurm_job.sh ]
then
    echo "Error: The genotype_selection directory: $DIR does not contain the script: 01_gp_slurm_job.sh. Are you sure this is the genomic_selection repo directory?"
    exit 107
fi
### Initialise the output directory which will contain all the output Rds files across populations and traits
if [ ! -d ${DIR}/output ]
then
    mkdir ${DIR}/output
fi
### Submit an array of jobs equivalent to the number of traits in the phenotype file
cd $DIR/
N_TRAITS=$(echo $(head -n1 $PHENOTYPE_DATA_TSV | awk '{print NF}') - 2 | bc)
N_POPS=$(tail -n+2 $PHENOTYPE_DATA_TSV | cut -f2 - | sort | uniq | wc -l)
sbatch --array 1-$(echo "${N_TRAITS} * ${N_POPS}" | bc) \
    01_gp_slurm_job-LUCERNE.sh \
        ${GENOTYPE_DATA_RDS} \
        ${PHENOTYPE_DATA_TSV} \
        ${KFOLDS} \
        ${NREPS} \
        ${DIR}
