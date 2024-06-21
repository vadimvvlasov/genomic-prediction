devtools::load_all()
# suppressPackageStartupMessages(library("gp"))
suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser(description="Genomic prediction cross-validation within and across populations. Please find the documentation within and interactive R session.")
parser$add_argument("-g", "--fname-geno",                                dest="fname_geno",                                      type="character")
parser$add_argument("-p", "--fname-pheno",                               dest="fname_pheno",                                     type="character")
parser$add_argument("-n", "--population",                                dest="population",                                      type="character")
parser$add_argument("-c", "--fname-covar",                               dest="fname_covar",                                     type="character", default=NULL)
parser$add_argument("-o", "--dir-output",                                dest="dir_output",                                      type="character", default=NULL)
parser$add_argument("-s", "--geno-fname-snp-list",                       dest="geno_fname_snp_list",                             type="character", default=NULL)
parser$add_argument("-x", "--geno-ploidy",                               dest="geno_ploidy",                                     type="integer", default=NULL)
parser$add_argument("--geno-bool-force-biallelic",                       dest="geno_bool_force_biallelic",                       type="logical", default=TRUE)
parser$add_argument("--geno-bool-retain-minus-one-alleles-per-locus",    dest="geno_bool_retain_minus_one_alleles_per_locus",    type="logical", default=TRUE)
parser$add_argument("--geno-min-depth",                                  dest="geno_min_depth",                                  type="integer", default=0)
parser$add_argument("--geno-max-depth",                                  dest="geno_max_depth",                                  type="integer", default=.Machine$integer.max)
parser$add_argument("--geno-maf",                                        dest="geno_maf",                                        type="double", default=0.01)
parser$add_argument("--geno-sdev-min",                                   dest="geno_sdev_min",                                   type="double", default=0.0001)
parser$add_argument("--geno-max-n-alleles",                              dest="geno_max_n_alleles",                              type="integer", default=NULL)
parser$add_argument("--geno-max-sparsity-per-locus",                     dest="geno_max_sparsity_per_locus",                     type="double", default=NULL)
parser$add_argument("--geno-frac-topmost-sparse-loci-to-remove",         dest="geno_frac_topmost_sparse_loci_to_remove",         type="double", default=NULL)
parser$add_argument("--geno-n-topmost-sparse-loci-to-remove",            dest="geno_n_topmost_sparse_loci_to_remove",            type="integer", default=NULL)
parser$add_argument("--geno-max-sparsity-per-sample",                    dest="geno_max_sparsity_per_sample",                    type="double", default=NULL)
parser$add_argument("--geno-frac-topmost-sparse-samples-to-remove",      dest="geno_frac_topmost_sparse_samples_to_remove",      type="double", default=NULL)
parser$add_argument("--geno-n-topmost-sparse-samples-to-remove",         dest="geno_n_topmost_sparse_samples_to_remove",         type="integer", default=NULL)
parser$add_argument("--pheno-sep",                                       dest="pheno_sep",                                       type="character", default="\t")
parser$add_argument("--pheno-header",                                    dest="pheno_header",                                    type="character", default=TRUE)
parser$add_argument("--pheno-idx-col-id",                                dest="pheno_idx_col_id",                                type="integer", default=1)
parser$add_argument("--pheno-idx-col-pop",                               dest="pheno_idx_col_pop",                               type="integer", default=2)
parser$add_argument("--pheno-idx-col-y",                                 dest="pheno_idx_col_y",                                 type="integer", default=3)
parser$add_argument("--pheno-na-strings",                                dest="pheno_na_strings",                                type="character", default=c("", "-", "NA", "na", "NaN", "missing", "MISSING"))
parser$add_argument("--pheno-bool-remove-outliers",                      dest="pheno_bool_remove_outliers",                      type="logical", default=TRUE)
parser$add_argument("--pheno-bool-remove-NA",                            dest="pheno_bool_remove_NA",                            type="logical", default=FALSE)
parser$add_argument("--bool-within",                                     dest="bool_within",                                     type="logical", default=TRUE)
parser$add_argument("--bool-across",                                     dest="bool_across",                                     type="logical", default=FALSE)
parser$add_argument("--n-folds",                                         dest="n_folds",                                         type="integer", default=10)
parser$add_argument("--n-reps",                                          dest="n_reps",                                          type="integer", default=10)
parser$add_argument("--vec-models-to-test",                              dest="vec_models_to_test",                              type="character", default=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C","gBLUP"))
parser$add_argument("--bool-parallel",                                   dest="bool_parallel",                                   type="logical", default=TRUE)
parser$add_argument("--max-mem-Gb",                                      dest="max_mem_Gb",                                      type="double", default=15)
parser$add_argument("--n-threads",                                       dest="n_threads",                                       type="integer", default=2)
parser$add_argument("--verbose",                                         dest="verbose",                                         type="logical", default=TRUE)
args = parser$parse_args()
time_ini = Sys.time()
print("               ,@@@@@@@,")
print("       ,,,.   ,@@@@@@/@@,  .oo8888o.")
print("    ,&%%&%&&%,@@@@@/@@@@@@,8888|88/8o")
print("   ,%&|%&&%&&%,@@@|@@@/@@@88|88888/88'")
print("   %&&%&%&/%&&%@@|@@/ /@@@88888|88888'")
print("   %&&%/ %&%%&&@@| V /@@' `88|8 `/88'")
print("   `&%| ` /%&'    |.|        | '|8'")
print("       |o|        | |         | |")
print("       |.|        | |         | |")
print("    ||/ ._|//_/__/  ,|_//__||/.  |_//__/_")
print("Performing genomic prediction cross-validation using")
print(paste0("     - genotype file: ", args$fname_geno))
print(paste0("     - phenotype file: ", args$fname_pheno))
print(paste0("     - covariate file: ", args$fname_covar))
print(paste0("     - population: ", args$population))
print(paste0("     - within population CV: ", args$bool_within))
print(paste0("     - across populations CV: ", args$bool_across))
print(paste0("     - with a total of ", args$n_threads, " threads available and "))
print(paste0("       a total memory of ", args$max_mem_Gb, " Gb."))
print(paste0("Start time: ", time_ini))
print("Input parameters:")
print(args)
fname_out_Rds = gp::gp(args=args)
time_fin = Sys.time()
time_duration_minutes = as.numeric(difftime(time_fin, time_ini, units="min"))
print("-----------------------------------------------------------")
print("-----------------------------------------------------------")
print("-----------------------------------------------------------")
print(paste0("End time: ", time_fin))
print(paste0(" Finished after ", time_duration_minutes, " minutes"))
if (methods::is(fname_out_Rds, "gpError")) {
    print("ERROR:")
    print(fname_out_Rds)
} else {
    print(paste0(" Output Rds file: ", fname_out_Rds))
}
print("        |   ^__^")
print("         |_ (oo)|_______")
print("            (__)|       )|/|")
print("                ||----w |")
print("                ||     ||")

### TEST ON LUCERNE
# args = list(
#     fname_geno='/group/pasture/Jeff/lucerne/workdir/FINAL-IMPUTED-noTrailingAllele-filteredSNPlist.Rds',
#     fname_pheno='/group/pasture/Jeff/lucerne/workdir/Lucerne_PhenomicsDB_2024-05-27-BiomassPredicted.tsv',
#     population="DB-MS-31-22-001",
#     fname_covar=NULL,
#     dir_output="outdir/lucerne",
#     geno_fname_snp_list=NULL,
#     geno_ploidy=NULL,
#     geno_bool_force_biallelic=TRUE,
#     geno_bool_retain_minus_one_alleles_per_locus=FALSE,
#     geno_min_depth=0,
#     geno_max_depth=.Machine$integer.max,
#     geno_maf=0.01,
#     geno_sdev_min=0.0001,
#     geno_max_n_alleles=NULL,
#     geno_max_sparsity_per_locus=NULL,
#     geno_frac_topmost_sparse_loci_to_remove=NULL,
#     geno_n_topmost_sparse_loci_to_remove=NULL,
#     geno_max_sparsity_per_sample=NULL,
#     geno_frac_topmost_sparse_samples_to_remove=NULL,
#     geno_n_topmost_sparse_samples_to_remove=NULL,
#     pheno_sep="\t",
#     pheno_header=TRUE,
#     pheno_idx_col_id=1,
#     pheno_idx_col_pop=2,
#     pheno_idx_col_y=3,
#     pheno_na_strings=c("", "-", "NA", "na", "NaN", "missing", "MISSING"),
#     pheno_bool_remove_NA=FALSE,
#     bool_within=TRUE,
#     bool_across=TRUE,
#     n_folds=5,
#     n_reps=1,
#     vec_models_to_test=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C","gBLUP"),
#     bool_parallel=TRUE,
#     max_mem_Gb=360,
#     n_threads=32,
#     verbose=TRUE
# )

### TEST ON GRAPE
# args = list(
#     fname_geno='grape.Rds',
#     fname_pheno='grape_pheno.txt',
#     population="g_1",
#     fname_covar=NULL,
#     dir_output="outdir",
#     geno_fname_snp_list=NULL,
#     geno_ploidy=NULL,
#     geno_bool_force_biallelic=TRUE,
#     geno_bool_retain_minus_one_alleles_per_locus=FALSE,
#     geno_min_depth=0,
#     geno_max_depth=.Machine$integer.max,
#     geno_maf=0.01,
#     geno_sdev_min=0.0001,
#     geno_max_n_alleles=NULL,
#     geno_max_sparsity_per_locus=NULL,
#     geno_frac_topmost_sparse_loci_to_remove=NULL,
#     geno_n_topmost_sparse_loci_to_remove=NULL,
#     geno_max_sparsity_per_sample=NULL,
#     geno_frac_topmost_sparse_samples_to_remove=NULL,
#     geno_n_topmost_sparse_samples_to_remove=NULL,
#     pheno_sep="\t",
#     pheno_header=TRUE,
#     pheno_idx_col_id=1,
#     pheno_idx_col_pop=2,
#     pheno_idx_col_y=4,
#     pheno_na_strings=c("", "-", "NA", "na", "NaN", "missing", "MISSING"),
#     pheno_bool_remove_NA=FALSE,
#     bool_within=TRUE,
#     bool_across=TRUE,
#     n_folds=5,
#     n_reps=2,
#     vec_models_to_test=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C","gBLUP"),
#     bool_parallel=TRUE,
#     max_mem_Gb=60,
#     n_threads=32,
#     verbose=TRUE
# )
