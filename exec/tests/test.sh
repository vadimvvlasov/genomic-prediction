#!/bin/bash
DIR='exec/tests/'
cd $DIR
time \
Rscript ../gp.R \
--fname-geno grape.Rds \
--fname-pheno grape_pheno.tsv \
--population g_1 \
--dir-output outdir \
--pheno-idx-col-y 4 \
--bool-within TRUE \
--bool-across TRUE \
--n-folds 5 \
--n-reps 2 \
--bool-parallel TRUE \
--max-mem-Gb 60 \
--n-threads 32 \
--verbose TRUE

cd exec
squeue -u $USER | sort
grep -i "err" slurm-*
grep "gpError" output/output-*/job_info-*.log
grep -A10 "Finished after" output/job_info-*.log
# head -n51 output/output-*/job_info-*.log | bat --wrap never
# bat --wrap never output/output-*/job_info-*.log
# grep "SLURM_ARRAY_TASK_ID = 36" output/output-*/job_info-*.log
# grep "SLURM_ARRAY_TASK_ID = 120" output/job_info-*.log

########################################
### COMPLETE SET OF INPUT PARAMETERS ###
########################################
# Rscript ../gp.R \
# --fname-geno= \
# --fname-pheno= \
# --population= \
# --fname-covar= \
# --dir-output= \
# --geno-fname-snp-list= \
# --geno-ploidy= \
# --geno-bool-force-biallelic= \
# --geno-bool-retain-minus-one-alleles-per-locus= \
# --geno-min-depth= \
# --geno-max-depth= \
# --geno-maf= \
# --geno-sdev-min= \
# --geno-max-n-alleles= \
# --geno-max-sparsity-per-locus= \
# --geno-frac-topmost-sparse-loci-to-remove= \
# --geno-n-topmost-sparse-loci-to-remove= \
# --geno-max-sparsity-per-sample= \
# --geno-frac-topmost-sparse-samples-to-remove= \
# --geno-n-topmost-sparse-samples-to-remove= \
# --pheno-sep= \
# --pheno-header= \
# --pheno-idx-col-id= \
# --pheno-idx-col-pop= \
# --pheno-idx-col-y= \
# --pheno-na-strings= \
# --pheno-bool-remove-NA= \
# --bool-within= \
# --bool-across= \
# --n-folds= \
# --n-reps= \
# --vec-models-to-test= \
# --bool-parallel= \
# --max-mem-Gb= \
# --n-threads= \
# --verbose= \
