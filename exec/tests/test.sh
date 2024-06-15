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

# ### Extract mean prediction-to-observation correlations within each population
# cd output/
# R
# vec_fnames_Rds = list.files(path=".", pattern="*.Rds")
# vec_traits = c()
# vec_populations = c()
# vec_training_size = c()
# vec_validation_size = c()
# vec_n_folds = c()
# vec_n_reps = c()
# vec_correlations = c()
# for (fname in vec_fnames_Rds) {
#     # fname = vec_fnames_Rds[1]
#     list_output = readRDS(fname)
#     vec_traits = c(vec_traits, list_output$TRAIT_NAME)
#     vec_populations = c(vec_populations, unique(list_output$METRICS_WITHIN_POP$pop_training))
#     vec_training_size = c(vec_training_size, mean(list_output$METRICS_WITHIN_POP$n_training))
#     vec_validation_size = c(vec_validation_size, mean(list_output$METRICS_WITHIN_POP$n_validation, na.rm=TRUE))
#     vec_n_folds = c(vec_n_folds, max(list_output$METRICS_WITHIN_POP$fold, na.rm=TRUE))
#     vec_n_reps = c(vec_n_reps, max(list_output$METRICS_WITHIN_POP$rep, na.rm=TRUE))
#     vec_idx = which(list_output$METRICS_WITHIN_POP$model == "Bayes_A")
#     vec_correlations = c(vec_correlations, mean(list_output$METRICS_WITHIN_POP$corr[vec_idx], na.rm=TRUE))
# }
# df = data.frame(
#     trait=vec_traits,
#     population=vec_populations,
#     training_size=vec_training_size,
#     validation_size=vec_validation_size,
#     n_fold=vec_n_folds,
#     n_rep=vec_n_reps,
#     correlation=vec_correlations)
# df = df[order(df$trait), ]
# vec_idx_sort_by_trait_name_length = unlist(lapply(df$trait, FUN=function(x){length(unlist(strsplit(x, "")))}))
# df = df[order(vec_idx_sort_by_trait_name_length, decreasing=FALSE), ]
# write.table(df, file="SUMMARY_WITHIN_POP_CORRELATIONS-Bayes_A.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# ### Extract lucerne additive effects for Bayes B in traits: 
# cd output/
# R
# vec_fnames_Rds = list.files(path=".", pattern="*.Rds")
# vec_requested_dates = c(
#     "biomass_SUM_2024a_20240213",
#     "biomass_AUT_2024_20240312",
#     "biomass_LSP_2023_20231023")

# for (date in vec_requested_dates) {
#     # date = vec_requested_dates[1]
#     vec_idx = which(grepl(date, vec_fnames_Rds))
#     for (idx in vec_idx) {
#         # idx = vec_idx[1]
#         list_output = readRDS(vec_fnames_Rds[idx])
#         df = data.frame(loc=gsub("\t", "|", names(list_output$ADDITIVE_GENETIC_EFFECTS$Bayes_B$b)), b=list_output$ADDITIVE_GENETIC_EFFECTS$Bayes_B$b)
#         trait_pop = paste(rev(rev(unlist(strsplit(vec_fnames_Rds[idx], "-"))[-1])[-1]), collapse="-")
#         fname_out = paste0("MARKER_EFFECTS-Bayes_B-", trait_pop, ".tsv")
#         write.table(df, file=fname_out, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#     }
# }

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
