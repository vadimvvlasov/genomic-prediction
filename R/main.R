### Parse use inputs
suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser(description="A workflow for genomic prediction/selection using multiple models with cross-validation performance metrics. Outputs a Rds file containing a list of lists each corresponding to a trait requested by `--phenotype-column-data`. Each list or trait consists of the following dataframes: METRICS_WITHIN_POP, YPRED_WITHIN_POP, METRICS_ACROSS_POP, YPRED_ACROSS_POP, and GENOMIC_PREDICTIONS.")
parser$add_argument("-g", "--genotype-file", dest="fname_rds_or_vcf", type="character", help="Genotype file in Rds or vcf format. The Rds file needs to contain a n samples x p loci numeric matrix with row and column names corresponding to the samples and loci names (any locus names format or chromosome or scaffold names, position and allele separated by underscores, e.g. 'chr1_12345_A'). The vcf file needs to have either the genotype call (**GT**) or allele depth (**AD**) field, where the AD field has precedence if both are present (i.e., allele frequencies rather than binary genotypes are assumed by default). The sample names should be the same as the names in the phenotype file.")
parser$add_argument("-p", "--phenotype-file", dest="fname_pheno", type="character", help="Phenotype data in text-delimited format, e.g. tab-delimited (`*.txt`) and comma-separated (`*.csv`), which may or may not have a header line. At least three columns are required: a column for the names of the samples, a column for the population or grouping ID of each sample, and a column for the numeric phenotype data (which can have missing data coded as whatever is convenient). Samples with known genotype and unknown phenotype data will have their phenotypes predicted. The best genomic prediction model (identified after k-fold cross-validation) for the population the samples belong to will be used. However, if the samples do not belong to populations with known genotype and phenotype data, then the best overall model will be used.")
parser$add_argument("-c", "--covariate-file", dest="fname_covar_rds", type="character", default=NULL, help="Covariate matrix where the row names overlap with the row names (sample names) of the genotype and phenotype files. Standard normalised per column, and set to be the first columns of the explanatory matrix for ridge, lasso, and elastic-net regression, while set as fixed effects (flat priors) in Bayesian models. [Default=NULL]")
parser$add_argument("--retain-minus-one-alleles-per-locus", dest="retain_minus_one_alleles_per_locus", type="logical", default=TRUE, help="Remove the minor allele from the genotype matrix? [Default=TRUE]")
parser$add_argument("--minimum-allele-frequency", dest="maf", type="double", default=0.01, help="Minimum allele frequency threshold for filtering out fixed or nearly fixed loci [Default=0.01]")
parser$add_argument("--phenotype-delimiter", dest="sep", type="character", default=",", help="Delimiter of the phenotype file [Default=',']")
parser$add_argument("--phenotype-header", dest="header", type="logical", default=TRUE, help="Does the phenotype file have a header line? [Default=TRUE]")
parser$add_argument("--phenotype-column-id", dest="idx_col_id", type="integer", default=1, help="Which column contains the sample names? [Default=1]")
parser$add_argument("--phenotype-column-pop", dest="idx_col_pop", type="integer", default=2, help="Which column contains the population or group names? [Default=2]")
parser$add_argument("--phenotype-column-data", dest="idx_col_y", type="character", default="3", help="Which column/s contain the phenotype data? If more than 1 column is requested, then separate them with commas, e.g. '3,4,5'.. [Default='3']")
parser$add_argument("--phenotype-missing-strings", dest="na_strings", type="character", default=c("", "-", "NA", "na", "NaN", "missing", "MISSING"), help="How were the missing data in the phenotype file coded? [Default=',-,NA,na,NaN,missing,MISSING'")
parser$add_argument("--populations-to-include", dest="vec_pop", type="character", default=c("all"), help="Which populations do you wish to include? [Default='all']")
parser$add_argument("--skip-leave-one-population-out-cross-validation", dest="skip_lopo_cv", type="logical", default=FALSE, help="Do you wish to skip leave-one-populatiop-out cross-validataion? [Default=FALSE]")
parser$add_argument("--models", dest="models_to_test", type="character", default=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C", "gBLUP"), help="Genomic prediction models to use [Default='ridge,lasso,elastic_net,Bayes_A,Bayes_B,Bayes_C'")
parser$add_argument("--k-folds", dest="k_folds", type="integer", default=10, help="Number of folds to perform in within population k-fold cross-validation [Default=10]")
parser$add_argument("--n-reps", dest="n_reps", type="integer", default=10, help="Number of reshuffled replications per fold of within population k-fold cross-validation [Default=3]")
parser$add_argument("--n-threads", dest="n_threads", type="integer", default=8, help="Number of computing threads to use in within and across population cross-validations [Default=8]")
parser$add_argument("--mem-mb", dest="mem_mb", type="integer", default=64000, help="Total memory in megabytes available for within and across population cross-validations [Default=64000]")
parser$add_argument("-o", "--output-file-prefix", dest="output_file_prefix", type="character", default="", help="Prefix of the filename of the output Rds file containing a list of metrics and predicted phenotype dataframes. The base name of the input genotype file will be used by default. [Default='']")
parser$add_argument("--verbose", dest="verbose", type="logical", default=TRUE, help="Do you wish to print detailed progress reports of the workflow? [Default=TRUE]")
args = parser$parse_args()
# args = list(
#     fname_rds_or_vcf="/group/pasture/Jeff/lucerne/workdir/FINAL-IMPUTED.Rds",
#     fname_pheno="/group/pasture/Jeff/lucerne/workdir/phenotype_for_groundTruthDryWeightGP-20240508.tsv",
#     fname_covar_rds=NULL,
#     retain_minus_one_alleles_per_locus=TRUE,
#     maf=0.05,
#     sep="\t",
#     header=TRUE,
#     idx_col_id=1,
#     idx_col_pop=2,
#     idx_col_y="8",
#     na_strings=c("", "-", "NA", "na", "NaN", "missing", "MISSING"),
#     vec_pop=c("DB-MS-31-22-001"),
#     skip_lopo_cv=TRUE,
#     models_to_test=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C", "gBLUP"),
#     k_folds=2,
#     n_reps=10,
#     n_threads=2,
#     mem_mb=64000,
#     output_file_prefix="output_file_prefix",
#     verbose=TRUE
# )
print("#########################################################################################################################")
print("♥ (｡•̀ᴗ-)✧(✿◠‿◠)(◕‿◕✿)(ง'̀-'́)ง(つ▀¯▀)つ(｡♥‿♥｡)(っ◔◡◔)っ ♥(｡•̀ᴗ-)✧(✿◠‿◠)(◕‿◕✿)(ง'̀-'́)ง(つ▀¯▀)つ(｡♥‿♥｡)(っ◔◡◔)っ ♥")
print("          GENOMIC PREDICTION CROSS-VALIDATION WITHIN (K-FOLD) AND ACROSS (LEAVE-ONE-POPULATION-OUT) POPULATIONS          ")
print("♥ (◕‿◕✿)(｡♥‿♥｡)(✿◠‿◠)(つ▀¯▀)つ(っ◔◡◔)っ ♥(｡•̀ᴗ-)✧(ง'̀-'́)ง(｡•̀ᴗ-)✧(｡•̀ᴗ-)✧(ง'̀-'́)ง(｡•̀ᴗ-)✧(っ◔◡◔)っ ♥(つ▀¯▀)つ(✿◠‿◠) ♥")
print("#########################################################################################################################")
### Load functions
args_dummy = commandArgs(trailingOnly=FALSE)
dir_src = dirname(gsub("^--file=", "", args_dummy[grepl("--file=", args_dummy)]))
# dir_src = "/group/pasture/Jeff/genomic_selection/src"
source(file.path(dir_src, "load.R"))
source(file.path(dir_src, "cross_validation.R"))
print("##################################################")
print("Check input data and parameters")

### Check if the genotype and phenotype files exist
if (file.exists(args$fname_rds_or_vcf) == FALSE) {
    print(paste0("[Error] The genotype file: ", args$fname_rds_or_vcf, " does not exist. Exiting."))
    quit()
}
if (file.exists(args$fname_pheno) == FALSE) {
    print(paste0("[Error] The phenotype file: ", args$fname_pheno, " does not exist. Exiting."))
    quit()
}
### Parse the list of column indexes referring to the traits which will be iteratively processed
vec_idx_col_y = c()
for (x in unlist(strsplit(args$idx_col_y, ","))) {
    vec_range = sort(as.numeric(unlist(strsplit(x, "-"))))
    if (length(vec_range)==1) {
        vec_idx_col_y = c(vec_idx_col_y, vec_range)
    } else {
        ini = head(vec_range, 1)
        fin = tail(vec_range, 1)
        vec_idx_col_y = c(vec_idx_col_y, ini:fin)
    }
}
vec_idx_col_y = sort(unique(vec_idx_col_y))
### Check if we have all the required columns in the phenotype file
tmp_pheno = read.table(args$fname_pheno, sep=args$sep, header=args$header, na.strings=args$na_strings)
idx_col_max = max(c(args$idx_col_id, args$idx_col_pop, vec_idx_col_y))
if (ncol(tmp_pheno) < idx_col_max) {
    print(paste0("[Error] Incompatible column/s specified for the phenotype file: ", args$fname_pheno, ". Please check the following phenotype-related parameters:"))
    print(paste0("--phenotype-file = ", args$fname_pheno))
    print(paste0("--phenotype-delimiter = ", args$sep))
    print(paste0("--phenotype-header = ", args$header))
    print(paste0("--phenotype-column-id = ", args$idx_col_id))
    print(paste0("--phenotype-column-pop = ", args$idx_col_pop))
    print(paste0("--phenotype-column-data = ", args$idx_col_y))
    print(paste0("--phenotype-missing-strings = ", paste(args$na_strings, collapse=", ")))
    print("Exiting.")
    quit()
}
### Check if the phenotype columns are numerics
for (j in vec_idx_col_y) {
    list_y_pop = fn_load_phenotype(fname_csv_txt=args$fname_pheno,
                                sep=args$sep,
                                header=args$header,
                                idx_col_id=args$idx_col_id,
                                idx_col_pop=args$idx_col_pop,
                                idx_col_y=j,
                                na.strings=args$na_strings)
    tryCatch(sum(list_y_pop$y), error=function(e) {
        print(paste0("[Error] Error parsing the phenotype column: ", j, ". Do you need a header line? Are the missing values correctly coded?"))
        print(paste0("--phenotype-file = ", args$fname_pheno))
        print(paste0("--phenotype-delimiter = ", args$sep))
        print(paste0("--phenotype-header = ", args$header))
        print(paste0("--phenotype-column-id = ", args$idx_col_id))
        print(paste0("--phenotype-column-pop = ", args$idx_col_pop))
        print(paste0("--phenotype-column-data = ", args$idx_col_y))
        print(paste0("--phenotype-missing-strings = ", paste(args$na_strings, collapse=", ")))
        print("Exiting.")
        quit()
    })
}
### Check if the output prefix contains a folder name which exists
if (dirname(args$output_file_prefix) != "") {
    if (dir.exists(dirname(args$output_file_prefix)) == FALSE) {
        print(paste0("[Error] The directory of the output file requested does not exist. Please modify --output-file-prefix: ", args$output_file_prefix, ". Exiting."))
        quit()
    }
}
### Parse missing phenotype data codes and model names
fn_sign = read.delim(file.path(dir_src, "models.R"))$X
idx = which(grepl("function\\(G, y, idx_training, idx_validation, other_params=", fn_sign) & !grepl("^#", fn_sign))
vec_available_models = gsub("fn_", "", unlist(lapply(strsplit(fn_sign[idx], " = "), FUN=function(x){x[1]})))
args$na_strings = unlist(strsplit(gsub(" ", "", args$na_strings), ","))
args$models_to_test = unlist(strsplit(gsub(" ", "", args$models_to_test), ","))
idx_unknown_model = which(!(args$models_to_test %in% vec_available_models))
if(length(idx_unknown_model) > 0) {
    print(paste0("[Error] Unknown model: ", args$models_to_test[idx_unknown_model]))
    print("Please select from: ")
    cat(paste0(paste(paste0("     - ", vec_available_models)), collapse="\n"), "\n")
    quit()
}
### Load genotype matrix
print("##################################################")
print("Load genotype & filter by minimum allele frequency")
X = fn_load_genotype(fname_rds_or_vcf=args$fname_rds_or_vcf,
                     retain_minus_one_alleles_per_locus=args$retain_minus_one_alleles_per_locus)
### Deal with missing data (stop and impute using an external tool or continue using mean value imputation)
if (sum(is.na(X)) > 0) {
    print("Please consider stopping this analysis to impute missing data using an appropriate imputation tool for your dataset.")
    print("Missing data detected and will be naively imputed using mean genotypes per allele/locus.")
    print("You have 30 seconds to stare at this message and cancel anytime before or after this pause.")
    Sys.sleep(30)
    for (j in 1:ncol(X)) {
        idx = which(is.na(X[,j]))
        if (length(idx) > 0) {
            mu = mean(X[,j], na.rm=TRUE)
            if (is.na(mu)) {
                mu = mean(X, na.rm=TRUE)
            }
            X[idx, j] = mu
        }
    }
}
### Filter loci by minimum allele frequency
X = fn_filter_loci(G=X, maf=args$maf, sdev_min=0.001)
if (ncol(X) < 100) {
    print("There are less than 100 loci left after filtering.")
    print(paste0("Please consider relaxing the minimum allele frequency: --minimum-allele-frequency=", args.maf, "."))
    quit()
}
### Initialise temporary output directory
print("##################################################")
print("Create a temporary output directory:")
dir_geno = dirname(args$fname_rds_or_vcf) ### Note: "." for the current directory and will only be an empty string if it is empty itself which should return and quit above
if (dir_geno == ".") {
    dir_geno = getwd()
}
dir_tmp = file.path(dir_geno, paste0("tmp_gsout-", gsub(" ", "-", gsub(":", "", gsub("[.]", "", paste0(as.character(Sys.time()), sample.int(1e6, size=1)))))))
dir.create(dir_tmp)
print(paste0("    - ", dir_tmp))
### Perform genomic prediction cross-validation within and across populations as well as predcitions per se of genotypes missing phenotype data
print("##################################################")
print("Genomic prediction cross-validations per trait:")
out = list()
for (idx_col_y in vec_idx_col_y) {
    gp = fn_within_across_perse_genomic_prediction(G=X, idx_col_y=idx_col_y, args=args, dir_tmp=dir_tmp)
    trait_name = gp$TRAIT_NAME
    eval(parse(text=paste0("out$`", trait_name, "` = gp")))
}
### Each item of the out list is:
# list(TRAIT_NAME          = character,
#      SUMMARY             = data.frame,
#      METRICS_WITHIN_POP  = data.frame,
#      YPRED_WITHIN_POP    = data.frame,
#      METRICS_ACROSS_POP  = data.frame,
#      YPRED_ACROSS_POP    = data.frame,
#      GENOMIC_PREDICTIONS = data.frame)
### Output Rds of list of data frames
suffix = paste0(gsub("/", "", format(Sys.time(), format="%D%H%M%S")), round(runif(n=1, min=1e6, max=9e6)))
if (args$output_file_prefix == "") {
    output_file_rds = paste0(gsub(".vcf.gz$", "", ignore.case=TRUE, gsub(".vcf$", "", ignore.case=TRUE, gsub(".rds$", "", ignore.case=TRUE, args$fname_rds_or_vcf))), "-GENOMIC_PREDICTION_CV_GEBV-", suffix, ".Rds")
} else {
    output_file_rds = paste0(args$output_file_prefix, "-", suffix, ".Rds")
}
saveRDS(out, file=output_file_rds)
### Clean-up temporary directory
unlink(dir_tmp, recursive=TRUE)
### Final message with a summary table of the best model per population per trait
print("##################################################")
print("Finished successfully!")
print("##################################################")
print("Summary table of the genomic prediction accuracies within populations:")
summary_table = do.call(rbind, lapply(out, FUN=function(x) {
    data.frame(trait=x$TRAIT_NAME, x$SUMMARY)
}))
rownames(summary_table) = NULL
print(summary_table)
print("##################################################")
print("Please find the output Rds file containing a list of lists of data frames (top level list refer to each trait used) consisting of:")
print(paste0("     - within population ", args$k_folds, "-fold cross-validation metrics and predicted-expected phenotypes (2 data frames)"))
print("     - leave-one-population-out cross-validation metrics and predicted-expected phenotypes, if there are more than 1 population with known genotype and phenotype data supplied (another 2 data frames)")
print("     - predicted phenotypes from samples with known genotype and unknown phenotype data (1 data frame)")
print(paste0("OUTPUT FILENAME: ", output_file_rds))
print("##################################################")
