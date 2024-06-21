dir = args(trailingOnly=TRUE)
# dir = "/group/pasture/Jeff/gp/inst/exec_Rscript/output/grape"
vec_fnames_Rds = list.files(path=dir, pattern="*.Rds")
out = list()
for (fname in vec_fnames_Rds) {
    # fname = vec_fnames_Rds[1]
    x = readRDS(fname)
    trait_name = x$TRAIT_NAME
    pop_name = x$POPULATION
    eval(parse(text=paste0("out$`", trait_name, "-", pop_name, "` = x")))
}

best_model_within = do.call(rbind, lapply(out, FUN=function(x) {
    df = x$METRICS_WITHIN_POP
    df_corr_mu = aggregate(corr ~ model + pop_validation, data=df, FUN=mean, na.rm=TRUE); colnames(df_corr_mu) = c("model", "pop_validation", "mu")
    df_corr_sd = aggregate(corr ~ model + pop_validation, data=df, FUN=sd, na.rm=TRUE); colnames(df_corr_sd) = c("model", "pop_validation", "sd")
    df_corr = merge(df_corr_mu, df_corr_sd, by=c("model", "pop_validation"))
    pop_validation = sort(unique(df_corr$pop_validation))
    df_out = data.frame(trait=x$TRAIT_NAME, pop_validation=pop_validation, model="", corr=NA, corr_sd=NA)
    for (i in 1:length(pop_validation)) {
        subdf = df_corr[df_corr$pop_validation==pop_validation[i], ]
        idx = which(subdf$mu == max(subdf$mu))[1]
        df_out$pop_validation[i] = pop_validation[i]
        df_out$model[i] = as.character(subdf$model[idx])
        df_out$corr[i] = subdf$mu[idx]
        df_out$corr_sd[i] = subdf$sd[idx]
    }
    ### Overall best model across all the populations in terms of power to correctly identify the top 10%
    df_corr_mu = aggregate(corr ~ model, data=df, FUN=mean, na.rm=TRUE); colnames(df_corr_mu) = c("model", "mu")
    df_corr_sd = aggregate(corr ~ model, data=df, FUN=sd, na.rm=TRUE); colnames(df_corr_sd) = c("model", "sd")
    df_corr = merge(df_corr_mu, df_corr_sd, by="model")
    df_corr = df_corr[which(df_corr$mu == max(df_corr$mu, na.rm=TRUE))[1], ]
    df_out = rbind(df_out, data.frame(trait=x$TRAIT_NAME, pop_validation="overall", model=df_corr$model, corr=df_corr$mu, corr_sd=df_corr$sd))
    return(df_out)
}))
rownames(best_model_within) = NULL
# print(best_model_within)
### Extract the best models in terms of power to detect the true top 10% samples
best_model_within_by_power = do.call(rbind, lapply(out, FUN=function(x) {
    df = x$METRICS_WITHIN_POP
    df_power_10_mu = aggregate(power_t10 ~ model + pop_validation, data=df, FUN=mean, na.rm=TRUE); colnames(df_power_10_mu) = c("model", "pop_validation", "mu")
    df_power_10_sd = aggregate(power_t10 ~ model + pop_validation, data=df, FUN=sd, na.rm=TRUE); colnames(df_power_10_sd) = c("model", "pop_validation", "sd")
    df_power_10 = merge(df_power_10_mu, df_power_10_sd, by=c("model", "pop_validation"))
    pop_validation = sort(unique(df_power_10$pop_validation))
    df_out = data.frame(trait=x$TRAIT_NAME, pop_validation=pop_validation, model="", power_t10=NA, power_t10_sd=NA)
    for (i in 1:length(pop_validation)) {
        subdf = df_power_10[df_power_10$pop_validation==pop_validation[i], ]
        idx = which(subdf$mu == max(subdf$mu))[1]
        df_out$pop_validation[i] = pop_validation[i]
        df_out$model[i] = as.character(subdf$model[idx])
        df_out$power_t10[i] = subdf$mu[idx]
        df_out$power_t10_sd[i] = subdf$sd[idx]
        ### Determine if the highest power is zero then set the model to NA
        if (max(subdf$mu) == 0.0) {
            df_out$model[i] = NA
        }
    }
    ### Overall best model across all the populations in terms of power to correctly identify the top 10%
    df_power_10_mu = aggregate(power_t10 ~ model, data=df, FUN=mean, na.rm=TRUE); colnames(df_power_10_mu) = c("model", "mu")
    df_power_10_sd = aggregate(power_t10 ~ model, data=df, FUN=sd, na.rm=TRUE); colnames(df_power_10_sd) = c("model", "sd")
    df_power_10 = merge(df_power_10_mu, df_power_10_sd, by="model")
    df_power_10 = df_power_10[which(df_power_10$mu == max(df_power_10$mu, na.rm=TRUE))[1], ]
    df_out = rbind(df_out, data.frame(trait=x$TRAIT_NAME, pop_validation="overall", model=df_power_10$model, power_t10=df_power_10$mu, power_t10_sd=df_power_10$sd))
    return(df_out)
}))
rownames(best_model_within_by_power) = NULL
# print(best_model_within_by_power)
### Extract the table of the best models according to leave-one-population-out cross-validations
best_model_across = do.call(rbind, lapply(out, FUN=function(x) {
    if (!is.na(x$METRICS_ACROSS_POP_LOPO[1])[1]) {
        df = x$METRICS_ACROSS_POP_LOPO
        pop_validation = sort(unique(df$pop_validation))
        if (is.null(pop_validation)) {
            df_out = data.frame(trait=x$TRAIT_NAME, pop_validation="", model="", corr=NA)
        } else {
            df_out = data.frame(trait=x$TRAIT_NAME, pop_validation=pop_validation, model="", corr=NA)
            for (i in 1:length(pop_validation)) {
                subdf = df[df$pop_validation==pop_validation[i], ]
                idx = which(subdf$corr == max(subdf$corr))[1]
                df_out$pop_validation[i] = pop_validation[i]
                df_out$model[i] = as.character(subdf$model[idx])
                df_out$corr[i] = subdf$corr[idx]
            }
        }
    } else {
        df_out = data.frame(trait=x$TRAIT_NAME, pop_validation="", model="", corr=NA)
    }
    return(df_out)
}))
best_model_across = best_model_across[!is.na(best_model_across$corr), , drop=FALSE]
rownames(best_model_across) = NULL
# print(best_model_across)
### Extract heritability estimates from the best model (in terms of corr)
heritability_estimates = do.call(rbind, lapply(out, FUN=function(x){
    df_h2 = aggregate(h2 ~ model + pop_validation, FUN=mean, data=x$METRICS_WITHIN_POP)
    df_corr = aggregate(corr ~ model + pop_validation, FUN=mean, data=x$METRICS_WITHIN_POP)
    df = merge(df_h2, df_corr, by=c("model", "pop_validation"))
    vec_pop = unique(df$pop_validation)
    vec_model = c()
    vec_h2 = c()
    for (p in vec_pop) {
        subdf = df[df$pop_validation==p, ]
        idx = which(subdf$corr==max(subdf$corr))[1]
        vec_model = c(vec_model, as.character(subdf$model[idx]))
        vec_h2 = c(vec_h2, subdf$h2[idx])
    }
    return(data.frame(trait=x$TRAIT_NAME, pop_validation=vec_pop, model=vec_model, h2=vec_h2))
}))
rownames(heritability_estimates) = NULL
# print(heritability_estimates)
### Performance of each model per trait across populations
performances_per_model_across_traits = reshape(do.call(rbind, lapply(out, FUN=function(x) {
    df = x$METRICS_WITHIN_POP
    agg = aggregate(corr ~ model, data=df, FUN=mean)
    agg$trait = x$TRAIT_NAME
    return(agg)

})), idvar="model", timevar="trait", direction="wide")
rownames(performances_per_model_across_traits) = NULL
# print(performances_per_model_across_traits)
### Extract predictions of entries without phenotype data across traits
predictions = do.call(rbind, lapply(out, FUN=function(x) {
    if (is.null(x$GENOMIC_PREDICTIONS)) {
        df = data.frame(trait=x$TRAIT_NAME, pop_validation=NA, entry=NA, y_pred=NA, top_model=NA, corr_from_kfoldcv=NA)
    } else {
        df = data.frame(trait=x$TRAIT_NAME, x$GENOMIC_PREDICTIONS)
    }
    df
}))
rownames(predictions) = NULL
# print(predictions)
### Output csv files
filename_best_model_within = file.path("output", paste0(prefix, "BEST_MODELS_WITHIN.csv"))
filename_best_model_within_by_power = file.path("output", paste0(prefix, "BEST_MODELS_WITHIN_BY_POWER.csv"))
filename_best_model_across = file.path("output", paste0(prefix, "BEST_MODELS_ACROSS.csv"))
filename_heritability_estimates = file.path("output", paste0(prefix, "HERITABILITY_ESTIMATES.csv"))
filename_performances_per_model_across_traits = file.path("output", paste0(prefix, "MODEL_PERFORMANCES_PER_TRAIT.csv"))
filename_predictions = file.path("output", paste0(prefix, "GENOMIC_PREDICTIONS.csv"))
write.table(best_model_within, file=filename_best_model_within, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")
write.table(best_model_within_by_power, file=filename_best_model_within_by_power, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")
write.table(best_model_across, file=filename_best_model_across, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")
write.table(heritability_estimates, file=filename_heritability_estimates, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")
write.table(performances_per_model_across_traits, file=filename_performances_per_model_across_traits, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")
write.table(predictions, file=filename_predictions, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=",")
print("Please find the output files:")
print(paste0("     - ", filename_best_model_within))
print(paste0("     - ", filename_best_model_within_by_power))
print(paste0("     - ", filename_best_model_across))
print(paste0("     - ", filename_heritability_estimates))
print(paste0("     - ", filename_performances_per_model_across_traits))
print(paste0("     - ", filename_predictions))
print("You may also find the R shiny app useful in plotting the results:")
print("     - go to http://shiny.science.depi.vic.gov.au/jp3h/test/")
print("     - enter the directory containing the GS_OUTPUT-*.Rds files within BASC (does not work on local files)")
print("     - select the plots you wish to make")
print("     - save them as desired")

