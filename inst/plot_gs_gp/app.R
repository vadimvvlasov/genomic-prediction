library(shiny)
library(shinyWidgets)
library(plotly)
library(bslib)
library(shinyFiles)

###########
### I/O ###
###########
fn_io_server = function(dir=NULL) {
	# dirname_root = "/"
	if (is.null(dir)) {
		### Dummy dataset to show at initialisation
      list_list_output = list(
		dummyTrait_dummyPop=list(
			TRAIT_NAME="dummyTrait", 
			POPULATION="dummyPop",
			METRICS_WITHIN_POP=data.frame(trait="dummyTrait", pop_training="dummyPop", pop_validation="dummyPop", model="dummyModel", corr=rnorm(20)),
			METRICS_ACROSS_POP_LOPO=data.frame(trait="dummyTrait", pop_training=paste0("dummyPop", 1:2), pop_validation=paste0("dummyPop", 2:1),  model=paste0("dummyModel", 1:5), corr=rnorm(10)),
			METRICS_ACROSS_POP_PAIRWISE=data.frame(trait="dummyTrait", pop_training=paste0("dummyPop", 1:2), pop_validation=paste0("dummyPop", 2:1), model=paste0("dummyModel", 1:5), corr=rnorm(10)),
			YPRED_WITHIN_POP=data.frame(id=rep(paste0("id", 1:10), times=5), model=rep(paste0("dummyModel", 1:5), each=10), y_true=rnorm(50), y_pred=rnorm(50)),
			YPRED_ACROSS_POP_LOPO=data.frame(id=paste0("id", 1:50), model=paste0("dummyModel", 1:5), y_true=rnorm(50), y_pred=rnorm(50), pop_training=paste0("dummyPop", 1:2)),
			YPRED_ACROSS_POP_PAIRWISE=data.frame(id=paste0("id", 1:50), model=paste0("dummyModel", 1:5), y_true=rnorm(50), y_pred=rnorm(50), pop_training=paste0("dummyPop", 1:2), pop_validation=paste0("dummyPop", 2:1)),
			GENOMIC_PREDICTION=data.frame(pop=paste0("dummyPop", 1:2), id=paste0("dummyEntry", 51:60), y_pred=rnorm(20), top_model="dummyModel", corr_from_kfoldcv=runif(20))
		)
      )
	  return(list_list_output)
	}
	vec_Rds = file.path(dir, list.files(path=dir, pattern=".*[Rr][Dd][Ss]$"))
	# vec_Rds = list.files(path="/group/pasture/Jeff/gp/inst/exec_Rscript/output/ryegrass", pattern=".*[Rr][Dd][Ss]$")
	list_list_output = list()
	for (i in 1:length(vec_Rds)) {
		# i = 1
		if (i == 1) {
			nth = "1st"
		} else if (i == 2) {
			nth = "2nd"
		} else if (i == 3) {
			nth = "3rd"
		} else {
			nth = paste(i, "th")
		}
		list_output = readRDS(vec_Rds[i])
		shiny::validate(
			shiny::need(class(list_output) == "list", paste0("Error: the ", nth, " Rds input is not the correct genomic_selection pipeline output (not a list)!")),
			shiny::need(length(list_output) == 12, paste0("Error: the ", nth, " Rds input is not the correct genomic_selection pipeline output (does not have 7 elements)!"))
		)
		shiny::validate(
			shiny::need(!is.null(list_output$TRAIT_NAME), paste0("Missing filed: TRAIT_NAME")),
			shiny::need(!is.null(list_output$POPULATION), paste0("Missing filed: POPULATION")),
			shiny::need(!is.null(list_output$METRICS_WITHIN_POP), paste0("Missing filed: METRICS_WITHIN_POP")),
			shiny::need(!is.null(list_output$YPRED_WITHIN_POP), paste0("Missing filed: YPRED_WITHIN_POP")),
			shiny::need(!is.null(list_output$METRICS_ACROSS_POP_BULK), paste0("Missing filed: METRICS_ACROSS_POP_BULK")),
			shiny::need(!is.null(list_output$YPRED_ACROSS_POP_BULK), paste0("Missing filed: YPRED_ACROSS_POP_BULK")),
			shiny::need(!is.null(list_output$METRICS_ACROSS_POP_PAIRWISE), paste0("Missing filed: METRICS_ACROSS_POP_PAIRWISE")),
			shiny::need(!is.null(list_output$YPRED_ACROSS_POP_PAIRWISE), paste0("Missing filed: YPRED_ACROSS_POP_PAIRWISE")),
			shiny::need(!is.null(list_output$METRICS_ACROSS_POP_LOPO), paste0("Missing filed: METRICS_ACROSS_POP_LOPO")),
			shiny::need(!is.null(list_output$YPRED_ACROSS_POP_LOPO), paste0("Missing filed: YPRED_ACROSS_POP_LOPO")),
			shiny::need(!is.null(list_output$GENOMIC_PREDICTIONS), paste0("Missing filed: GENOMIC_PREDICTIONS")),
			shiny::need(!is.null(list_output$ADDITIVE_GENETIC_EFFECTS), paste0("Missing filed: ADDITIVE_GENETIC_EFFECTS"))
		)
		trait = list_output$TRAIT_NAME
		pop = list_output$POPULATION
		eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "` = list_output")))
	}
	return(list_list_output)
}

##################
### WITHIN POP ###
##################
fn_within_ui = function() {
  card(
    card_header(h1(strong("Input files and main plot"), style="font-size:21px; text-align:left")),
    min_height="840px",
    layout_sidebar(
      sidebar=sidebar(
        width=500,
        shinyDirButton('dir', label='Select output directory', title='Please select the directory containing the genomic prediction output', multiple=FALSE),
        shinyWidgets::pickerInput(inputId="within_trait", label="Filter by trait:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="within_pop", label="Filter by population:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="within_model", label="Filter by model:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="within_metric", label="Use the genomic prediction accuracy metric:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="within_group_by", label="Group by:", choices=c("trait", "pop_training", "model"), selected=c("trait", "pop", "model"), multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="within_sort_by", label="Sort by:", choices=c("increasing_mean", "decreasing_mean", "increasing_median", "decreasing_median", "alphabetical"), selected="increasing_mean", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::materialSwitch(inputId="within_box_with_labels", label="Mean-labelled boxplot", value=FALSE, status="primary", right=FALSE)
      ),
      mainPanel(
        width=750,
        plotlyOutput(outputId="within_plot_violin_or_box"),
        verbatimTextOutput(outputId="within_data_sizes"),
		downloadButton("within_download", "Download")
      )
    )
  )
}

fn_within_prep_server = function(input, list_list_output) {
	### Define the grouping, i.e. trait and/or pop and/or model
		METRICS_WITHIN_POP = NULL
		for (list_output in list_list_output) {
		df_metrics_within_pop = list_output$METRICS_WITHIN_POP
		df_metrics_within_pop$trait = list_output$TRAIT_NAME
		if (is.null(METRICS_WITHIN_POP)) {
			METRICS_WITHIN_POP = df_metrics_within_pop
		} else {
			METRICS_WITHIN_POP = rbind(METRICS_WITHIN_POP, df_metrics_within_pop)
		}
		}
		idx = (METRICS_WITHIN_POP$trait %in% input$within_trait) & (METRICS_WITHIN_POP$pop_training %in% input$within_pop) & (METRICS_WITHIN_POP$model %in% input$within_model)
		grouping = c()
		x_label_one_level = c()
		x_label_one_level_names = c()
		x_label_multiple_levels = c()
		for (g in input$within_group_by) {
		g_names = as.character(eval(parse(text=paste0("METRICS_WITHIN_POP$", g, "[idx]"))))
		if ((length(unique(g_names)) == 1) & (length(list_list_output)>1)) {
			x_label_one_level = c(x_label_one_level, g)
			x_label_one_level_names = c(x_label_one_level_names, g_names[1])
			next
		} else {
			x_label_multiple_levels = c(x_label_multiple_levels, g)
		}
		if (length(grouping)==0) {
			grouping = g_names
		} else {
			grouping = paste0(grouping, "_x_", g_names)
		}
		}
		if (length(x_label_one_level)==0) {
		title = paste(x_label_multiple_levels, collapse=" x ")
		} else {
		title = paste0(
			paste(x_label_multiple_levels, collapse=" x "),
			"\n",
			paste(
			paste0(x_label_one_level, ": ", x_label_one_level_names), 
			collapse=" | "
			)
		)
		}
		### Sort by mean, median or alphabetical (default)
		shiny::validate(
		shiny::need(grouping > 0, "Error: Please select a grouping factor ('Group by:' selection) with more than two levels selected.")
		)
		df = data.frame(
		grouping=grouping,
		metric=eval(parse(text=paste0("METRICS_WITHIN_POP$", input$within_metric, "[idx]")))
		)
		if (input$within_sort_by == "increasing_mean") {
		df_agg = aggregate(metric ~ grouping, data=df, FUN=mean, na.rm=TRUE)
		df$grouping = factor(df$grouping, levels=as.character(df_agg$grouping[order(df_agg$metric, decreasing=FALSE)]))
		}
		if (input$within_sort_by == "decreasing_mean") {
		df_agg = aggregate(metric ~ grouping, data=df, FUN=mean, na.rm=TRUE)
		df$grouping = factor(df$grouping, levels=as.character(df_agg$grouping[order(df_agg$metric, decreasing=TRUE)]))
		}
		if (input$within_sort_by == "increasing_median") {
		df_agg = aggregate(metric ~ grouping, data=df, FUN=median, na.rm=TRUE)
		df$grouping = factor(df$grouping, levels=as.character(df_agg$grouping[order(df_agg$metric, decreasing=FALSE)]))
		}
		if (input$within_sort_by == "decreasing_median") {
		df_agg = aggregate(metric ~ grouping, data=df, FUN=median, na.rm=TRUE)
		df$grouping = factor(df$grouping, levels=as.character(df_agg$grouping[order(df_agg$metric, decreasing=TRUE)]))
		}
		return(df)
}

fn_within_violin_server = function(input, df) {
		### Plot
		if (!input$within_box_with_labels) {
			p = plot_ly(data=df,
				y=~metric,
				x=~grouping,
				type='violin',
				box=list(visible=TRUE),
				meanline = list(visible=TRUE),
				split=~grouping
			)
		} else {
			list_bp = boxplot(metric ~ grouping, data=df, plot=FALSE)
			n_groups = ncol(list_bp$stats)
			df_agg = aggregate(metric ~ grouping, data=df, FUN=mean, na.rm=TRUE)
			colnames(df_agg)[2] = "metric_mean"
			df = merge(merge(df, df_agg, sort=FALSE, by="grouping"),
				data.frame(grouping=list_bp$names, pos_x=seq(0, n_groups+1, length=n_groups+1)[1:n_groups], pos_y=1.05*max(df$metric, na.rm=TRUE)),
				sort=FALSE, by="grouping")
			p = plot_ly(data=df,
				y=~metric,
				x=~grouping,
				type='box',
				boxmean=TRUE,
				split=~grouping
			)
			p = p %>% add_annotations(
				x=~grouping,
				y=~pos_y, 
				text=~round(metric_mean, 4),
				font=list(size=11),
				showarrow=FALSE
			)
		}
		p = p %>% layout(
			title=title,
			yaxis=list(title=input$within_metric),
			xaxis=list(title="")
		)
		p = p %>% config(toImageButtonOptions = list(format = "svg"))
		return(p)
}

####################################################################################################
### Front-end
####################################################################################################
ui <- page_fillable(
  titlePanel("PLOT(genomic_selection)"),
  fn_within_ui()
)
####################################################################################################
### Back-end
####################################################################################################
server = function(input, output, session) {
	### Extract data and update the selections
	data = reactive({
		dirname_root = "/group/pasture/Jeff/gp/inst/exec_Rscript/output"
		shinyFiles::shinyDirChoose(input, id='dir', roots=c(root=dirname_root), filetypes=c('rds', 'Rds', 'RDS'), session=session)
		dir = as.character(parseDirPath(c(root=dirname_root), input$dir))
		list_list_output = fn_io_server(dir=dir)
		### Update selections
		vec_traits = c()
		vec_populations = c()
		vec_models = c()
		vec_metrics = c()
		vec_entry_names = c()
		vec_pop_names = c()
		for (x in list_list_output) {
		# x = list_list_output[[1]]
		vec_traits = c(vec_traits, as.character(x$TRAIT_NAME))
		vec_populations = c(vec_populations, as.character(x$POPULATION))
		vec_models = c(vec_models, unique(as.character(x$METRICS_WITHIN_POP$model)))
		df_tmp = x$METRICS_WITHIN_POP
		for (j in 1:ncol(df_tmp)) {
			if (is.numeric(df_tmp[,j])) {
			vec_metrics = c(vec_metrics, colnames(df_tmp)[j])
			}
		}
		}
		vec_traits = sort(unique(vec_traits))
		vec_populations = sort(unique(vec_populations))
		vec_models = sort(unique(vec_models))
		vec_metrics = sort(unique(vec_metrics))
		updatePickerInput(session, "within_trait", choices=vec_traits, selected=vec_traits[1])
		updatePickerInput(session, "within_pop", choices=vec_populations, selected=vec_populations[1])
		updatePickerInput(session, "within_model", choices=vec_models, selected=vec_models)
		updatePickerInput(session, "within_metric", choices=vec_metrics, selected=vec_metrics[1])
		updatePickerInput(session, "within_scat_trait_x_pop", choices=names(list_list_output), selected=names(list_list_output)[1])
		updatePickerInput(session, "within_scat_model", choices=vec_models, selected=vec_models[1])
		return(list_list_output)
	})
	### Within population: violin or box plot
	output$within_plot_violin_or_box = renderPlotly({
		list_list_output = data()
		fn_within_violin_server(input=input, df=fn_within_prep_server(input, list_list_output))
	})
	### Within population: download selected raw and aggregated data based on user selections
	output$within_download = downloadHandler(
		filename = function() {
			paste0("gp_withinPop_tables-", format(Sys.time(),'%Y%m%d%a%H%M%S'), ".zip")
		},
		content = function(file) {
			list_list_output = data()
			df = fn_within_prep_server(input, list_list_output)
			df_agg_mean = aggregate(metric ~ grouping, data=df, FUN=mean, na.rm=TRUE); colnames(df_agg_mean)[2] = "mean"
			df_agg_median = aggregate(metric ~ grouping, data=df, FUN=median, na.rm=TRUE); colnames(df_agg_median)[2] = "median"
			df_agg_min = aggregate(metric ~ grouping, data=df, FUN=min, na.rm=TRUE); colnames(df_agg_min)[2] = "min"
			df_agg_max = aggregate(metric ~ grouping, data=df, FUN=max, na.rm=TRUE); colnames(df_agg_max)[2] = "max"
			df_agg_var = aggregate(metric ~ grouping, data=df, FUN=var, na.rm=TRUE); colnames(df_agg_var)[2] = "var"
			df_agg_nan = aggregate(metric ~ grouping, data=df, FUN=function(x){sum(is.na(x))}); colnames(df_agg_nan)[2] = "n_missing"
			df_stats = merge(df_agg_mean, 
				merge(df_agg_median, 
				merge(df_agg_min, 
				merge(df_agg_max, 
				merge(df_agg_var, df_agg_nan, by="grouping"), by="grouping"), by="grouping"), by="grouping"), by="grouping")
			dir_tmp = tempdir()
			vec_fnames_tsv = unlist(lapply(list(df=df, df_stats=df_stats), FUN=function(x){
				if (ncol(x) == 2) {
					fname_tsv = "raw_selected_data.tsv"
					write.table(df, file=file.path(dir_tmp, fname_tsv), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
				} else {
					fname_tsv = "aggregated_stats.tsv"
					write.table(df_stats, file=file.path(dir_tmp, fname_tsv), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
				}
				return(fname_tsv)
			}))
			zip::zip(zipfile=file, files=vec_fnames_tsv, root=dir_tmp)
		}
	)
}
####################################################################################################
### Serve the app
####################################################################################################
options(shiny.maxRequestSize = 50 * 1024^2)
shinyApp(ui = ui, server = server)
