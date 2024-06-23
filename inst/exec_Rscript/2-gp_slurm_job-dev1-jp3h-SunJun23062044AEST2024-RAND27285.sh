#!/bin/bash
#SBATCH --job-name='test'
#SBATCH --account='dbiopast2'   ### EDIT ME: Pick the appropriate account name, e.g. dbiopast1 or dbiopast2
#SBATCH --ntasks=1              ### LEAVE ME:Request a single task as we will be submitting this as an array job where each job corresponds to a trait
#SBATCH --cpus-per-task=8      ### EDIT ME: Parallelisation across replications, folds and models (more cpu means faster execution time but probably longer time to wait for the Slurm scheduler to find resources to allocate to the job)
#SBATCH --mem=64G              ### EDIT ME: Proportional to the input data (will need to test the appropriate memory required, hint use `seff ${JOBID}`)
#SBATCH --time=0-2:0:00         ### EDIT ME: Proportional to the input data, number of folds, replications, and models to be used
###################################################################################################
### Edit the Slurm settings above to match your requirements. 
###################################################################################################

###################################################################################################
### The variables below will be exported from `00_gs_slurm_job_wrapper.sh`:
###################################################################################################
### Input variables (use the absolute path to files to be precise)
### (1) R matrix object with n rows corresponding to samples, and p columns corresponding to the markers or loci. 
###     Should have no missing data or else will be imputed via mean value imputation.
GENOTYPE_DATA_RDS=$1
### (2) Tab-delimited phenotype file where column 1: sample names, column 2: population name, columns 3 and so on refers to the phenotype values of one trait per column.
###     Headers for the columns should be named appropriately, e.g. ID, POP, TRAIT1, TRAIT2, etc.
###     Missing values are allowed for samples whose phenotypes will be predicted by the best model identified within the population they belong to.
###     Missing values may be coded as empty cells, -, NA, na, NaN, missing, and/or MISSING.
PHENOTYPE_DATA_TSV=$2
### (3) Number of folds for k-fold cross-validation.
KFOLDS=$3
### (4) Number of replications of the k-fold cross-validation each representing a random sorting of the samples hence yielding different ways of partitioning the data.
NREPS=$4
### (5) Full path to the location of the executable Rscript gp.R
DIR_SRC=$5
### (5) Full path to the location of the executable Rscript gp.R
DIR_OUT=$6

###################################################################################################
### Edit the code below, if and only if you have read the documentation or familiar with `src/*.R`:
###################################################################################################
### Define the trait and population to include
N_POPS=$(tail -n+2 $PHENOTYPE_DATA_TSV | cut -f2 - | sort | uniq | wc -l)
TRAIT_IDX=$(echo "((${SLURM_ARRAY_TASK_ID}-1) / ${N_POPS}) + 1" | bc)
POP_IDX=$(echo "${SLURM_ARRAY_TASK_ID} % ${N_POPS}" | bc)
if [ "${POP_IDX}" -eq 0 ]
then
    POP_IDX=${N_POPS}
fi
COLUMN_ID=$(echo 2 + ${TRAIT_IDX} | bc)
TRAIT=$(head -n1 $PHENOTYPE_DATA_TSV | cut -f${COLUMN_ID})
POP=$(tail -n+2 $PHENOTYPE_DATA_TSV | cut -f2 - | sort | uniq | head -n${POP_IDX} | tail -n1)
### Skip leave-one-population-out cross-validation if there is only one population
if [ "${N_POPS}" -eq 1 ]
then
    BOOL_ACROSS=FALSE
else
    if [ "${POP_IDX}" -eq 1 ]
    then
        BOOL_ACROSS=TRUE
    else
        BOOL_ACROSS=FALSE
    fi
fi
### Output directories
DIR_OUT_MAIN=${DIR_OUT}/output
DIR_OUT_SUB=${DIR_OUT_MAIN}/output-${TRAIT}-${POP}
if [ ! -d DIR_OUT_MAIN ]
then
    mkdir $DIR_OUT_MAIN
fi
mkdir $DIR_OUT_SUB
### Log messages
echo JOB_${SLURM_ARRAY_TASK_ID}-TRAIT_${TRAIT}-POP_${POP} > ${DIR_OUT_SUB}/job_info-${TRAIT}-${POP}.log
echo "==========================================
-------------------------------------------
	    Job Info
-------------------------------------------
SLURM_JOB_ID = $SLURM_JOB_ID
SLURM_JOB_NAME = $SLURM_JOB_NAME
SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST
SLURM_SUBMIT_HOST = $SLURM_SUBMIT_HOST
SLURM_SUBMIT_DIR = $SLURM_SUBMIT_DIR
SLURM_NTASKS = $SLURM_NTASKS
SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID
SLURM_MEM_PER_NODE = $(echo "$SLURM_MEM_PER_NODE / (2^10)" | bc) GB
SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK
-------------------------------------------
	    Variables
-------------------------------------------
GENOTYPE_DATA_RDS	: $GENOTYPE_DATA_RDS
PHENOTYPE_DATA_TSV	: $PHENOTYPE_DATA_TSV
KFOLDS		        : $KFOLDS
NREPS		        : $NREPS
TRAIT		        : $TRAIT
POPULATION          : $POP
-------------------------------------------
	    Output directory
-------------------------------------------
${DIR_OUT_SUB}
==========================================" >> ${DIR_OUT_SUB}/job_info-${TRAIT}-${POP}.log

### Load the conda environment
module load Miniconda3/22.11.1-1
conda init bash
source ~/.bashrc
conda activate genomic_selection

### Run within and across population replicated k-fold cross-validation and prediction of missing phenotypes
time \
Rscript ${DIR_SRC}/gp.R \
    --fname-geno $GENOTYPE_DATA_RDS \
    --fname-pheno $PHENOTYPE_DATA_TSV \
    --population $POP \
    --dir-output $DIR_OUT_SUB \
    --pheno-idx-col-y $COLUMN_ID \
    --bool-within TRUE \
    --bool-across $BOOL_ACROSS \
    --n-folds $KFOLDS \
    --n-reps $NREPS \
    --bool-parallel TRUE \
    --max-mem-Gb $(echo "$SLURM_MEM_PER_NODE / (2^10)" | bc) \
    --n-threads $SLURM_CPUS_PER_TASK \
    --verbose TRUE >> ${DIR_OUT_SUB}/job_info-${TRAIT}-${POP}.log
### Clean-up
mv ${DIR_OUT_SUB}/GENOMIC_PREDICTIONS_OUTPUT-*.Rds ${DIR_OUT_MAIN}
mv ${DIR_OUT_SUB}/job_info-${TRAIT}-${POP}.log ${DIR_OUT_MAIN}
rm -R ${DIR_OUT_SUB}