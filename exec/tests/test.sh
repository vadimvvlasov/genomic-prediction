#!/bin/bash

#############
### GRAPE ###
#############
DIR='gp/exec/tests/'
cd $DIR
time \
Rscript ../gp.R \
--fname-geno grape.rds \
--fname-pheno grape_pheno.txt \
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

################
### RYEGRASS ###
################
DIR='/group/pasture/Jeff/gp/exec/tests'
cd $DIR
fname_geno='/group/pasture/Jeff/ryegrass/workdir/STR_NUE_WUE_HS-1717536141.3435302.3200855812-IMPUTED.tsv'
fname_pheno='/group/pasture/Jeff/ryegrass/workdir/Lucerne_PhenomicsDB_2024-05-27-Biomass.tsv'
n_traits=$(head -n1 $fname_pheno | awk '{print NF}')
touch logfile
time \
for idx_pheno in $(seq 3 $n_traits)
do
    time \
    Rscript ../gp.R \
        --fname-geno $fname_geno \
        --fname-pheno $fname_pheno \
        --population pop_1 \
        --dir-output outdir \
        --pheno-idx-col-y $idx_pheno \
        --bool-within TRUE \
        --bool-across FALSE \
        --n-folds 5 \
        --n-reps 1 \
        --bool-parallel TRUE \
        --max-mem-Gb 60 \
        --n-threads 32 \
        --verbose TRUE >> logfile
done
tail logfile
grep -A1 "ERROR:" logfile
ls -lhtr outdir/*.Rds


###############
### LUCERNE ###
###############
DIR='/group/pasture/Jeff/gp/exec/tests'
cd $DIR
fname_geno='/group/pasture/Jeff/lucerne/workdir/FINAL-IMPUTED-noTrailingAllele-filteredSNPlist.Rds'
fname_pheno='/group/pasture/Jeff/lucerne/workdir/Lucerne_PhenomicsDB_2024-05-27-BiomassPredicted.tsv'
n_traits=$(head -n1 $fname_pheno | awk '{print NF}')
touch outdir/lucerne/logfile
time \
for idx_pheno in $(seq 3 $n_traits)
do
    time \
    Rscript ../gp.R \
        --fname-geno $fname_geno \
        --fname-pheno $fname_pheno \
        --population "DB-MS-31-22-001" \
        --dir-output outdir/lucerne \
        --pheno-idx-col-y $idx_pheno \
        --bool-within TRUE \
        --bool-across TRUE \
        --n-folds 5 \
        --n-reps 1 \
        --bool-parallel TRUE \
        --max-mem-Gb 360 \
        --n-threads 32 \
        --verbose TRUE >> outdir/lucerne/logfile
done
tail outdir/lucerne/logfile
grep -A1 "ERROR:" outdir/lucerne/logfile
ls -lhtr outdir/lucerne/*.Rds


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
