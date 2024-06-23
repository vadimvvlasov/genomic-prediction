#!/bin/bash
### Unique name of the run given a specific instance of config.txt
RUN_NAME=$(hostname)-${USER}-$(date | sed 's/ //g' | sed 's/://g')-RAND${RANDOM}
### Extract config variables
CONFIG_GENO=$(sed "s/\"/'/g" config.txt | sed -n '1p')
CONFIG_PHENO=$(sed "s/\"/'/g" config.txt | sed -n '2p')
CONFIG_KFOLDS=$(sed "s/\"/'/g" config.txt | sed -n '3p')
CONFIG_NREPS=$(sed "s/\"/'/g" config.txt | sed -n '4p')
CONFIG_DIR_OUT=$(sed "s/\"/'/g" config.txt | sed -n '5p')
CONFIG_JOB_NAME=$(sed "s/\"/'/g" config.txt | sed -n '6p')
CONFIG_ACCOUNT_NAME=$(sed "s/\"/'/g" config.txt | sed -n '7p')
CONFIG_NTASKS=$(sed "s/\"/'/g" config.txt | sed -n '8p')
CONFIG_NCPUS=$(sed "s/\"/'/g" config.txt | sed -n '9p')
CONFIG_MEM=$(sed "s/\"/'/g" config.txt | sed -n '10p')
CONFIG_TIME_LIMIT=$(sed "s/\"/'/g" config.txt | sed -n '11p')
### Create the checks and submission scripts using the config variables
sed "s|GENOTYPE_DATA_RDS=\${DIR_SRC}/input/test_geno.Rds|$CONFIG_GENO|g" 1-checks_and_submision.sh | \
    sed "s|PHENOTYPE_DATA_TSV=\${DIR_SRC}/input/test_pheno.tsv|$CONFIG_PHENO|g" | \
    sed "s|KFOLDS=5|$CONFIG_KFOLDS|g" | \
    sed "s|NREPS=3|$CONFIG_NREPS|g" | \
    sed "s|DIR_OUT=\${DIR_SRC}|$CONFIG_DIR_OUT|g" | \
    sed "s|2-gp_slurm_job.sh|2-gp_slurm_job-${RUN_NAME}.sh|g" \
> 1-checks_and_submision-${RUN_NAME}.sh
### Create the slurm job scripts using the config variables
sed "s|SBATCH --job-name='GS'|$CONFIG_JOB_NAME|g" 2-gp_slurm_job.sh | \
    sed "s|SBATCH --account='dbiopast1'|$CONFIG_ACCOUNT_NAME|g" | \
    sed "s|SBATCH --ntasks=1|$CONFIG_NTASKS|g" | \
    sed "s|SBATCH --cpus-per-task=16|$CONFIG_NCPUS|g" | \
    sed "s|SBATCH --mem=100G|$CONFIG_MEM|g" | \
    sed "s|SBATCH --time=1-0:0:00|$CONFIG_TIME_LIMIT|g" \
> 2-gp_slurm_job-${RUN_NAME}.sh
### Check input and submit the slurm job
chmod +x 1-checks_and_submision-${RUN_NAME}.sh
chmod +x 2-gp_slurm_job-${RUN_NAME}.sh
./1-checks_and_submision-${RUN_NAME}.sh
### Clean-up after tests
# rm *RAND*.sh