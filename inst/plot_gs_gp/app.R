library(bslib)
library(shiny)
library(shinyWidgets)
library(plotly)
library(shinyFiles)
library(shinycssloaders)

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
		# print(i)
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
			shiny::need(length(list_output) >= 12, paste0("Error: the ", nth, " Rds input is not the correct genomic_selection pipeline output (does not have 12 elements)!"))
		)
		shiny::validate(
			shiny::need(!is.null(list_output$TRAIT_NAME), paste0("Missing filed: TRAIT_NAME")),
			shiny::need(!is.null(list_output$POPULATION), paste0("Missing filed: POPULATION")),
			shiny::need(!is.null(list_output$METRICS_WITHIN_POP), paste0("Missing filed: METRICS_WITHIN_POP")),
			shiny::need(!is.null(list_output$YPRED_WITHIN_POP), paste0("Missing filed: YPRED_WITHIN_POP")),
			# shiny::need(!is.null(list_output$METRICS_ACROSS_POP_BULK), paste0("Missing filed: METRICS_ACROSS_POP_BULK")),
			# shiny::need(!is.null(list_output$YPRED_ACROSS_POP_BULK), paste0("Missing filed: YPRED_ACROSS_POP_BULK")),
			# shiny::need(!is.null(list_output$METRICS_ACROSS_POP_PAIRWISE), paste0("Missing filed: METRICS_ACROSS_POP_PAIRWISE")),
			# shiny::need(!is.null(list_output$YPRED_ACROSS_POP_PAIRWISE), paste0("Missing filed: YPRED_ACROSS_POP_PAIRWISE")),
			# shiny::need(!is.null(list_output$METRICS_ACROSS_POP_LOPO), paste0("Missing filed: METRICS_ACROSS_POP_LOPO")),
			# shiny::need(!is.null(list_output$YPRED_ACROSS_POP_LOPO), paste0("Missing filed: YPRED_ACROSS_POP_LOPO")),
			# shiny::need(!is.null(list_output$GENOMIC_PREDICTIONS), paste0("Missing filed: GENOMIC_PREDICTIONS")),
			shiny::need(!is.null(list_output$ADDITIVE_GENETIC_EFFECTS), paste0("Missing field: ADDITIVE_GENETIC_EFFECTS"))
		)
		trait = list_output$TRAIT_NAME
		pop = list_output$POPULATION
		if (!is.null(eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "`"))))) {
			if (!is.na(head(list_output$METRICS_ACROSS_POP_BULK, n=1)[1])) {
				eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "`$METRICS_ACROSS_POP_BULK = list_output$METRICS_ACROSS_POP_BULK")))
				eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "`$YPRED_ACROSS_POP_BULK = list_output$YPRED_ACROSS_POP_BULK")))
				eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "`$METRICS_ACROSS_POP_PAIRWISE = list_output$METRICS_ACROSS_POP_PAIRWISE")))
				eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "`$YPRED_ACROSS_POP_PAIRWISE = list_output$YPRED_ACROSS_POP_PAIRWISE")))
				eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "`$METRICS_ACROSS_POP_LOPO = list_output$METRICS_ACROSS_POP_LOPO")))
				eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "`$YPRED_ACROSS_POP_LOPO = list_output$YPRED_ACROSS_POP_LOPO")))
			}
			if (!is.na(head(list_output$METRICS_WITHIN_POP, n=1)[1])) {
				eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "`$METRICS_WITHIN_POP = list_output$METRICS_WITHIN_POP")))
				eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "`$YPRED_WITHIN_POP = list_output$YPRED_WITHIN_POP")))
				eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "`$GENOMIC_PREDICTIONS = list_output$GENOMIC_PREDICTIONS")))
				eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "`$ADDITIVE_GENETIC_EFFECTS = list_output$ADDITIVE_GENETIC_EFFECTS")))
			}
		} else {
			eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "` = list_output")))
		}
	}
	return(list_list_output)
}

#########################
### WITHIN POPULATION ###
#########################
fn_within_violin_ui = function() {
	card(
		card_header(h1(strong("Input files and within population prediction accuracies"), style="font-size:21px; text-align:left")),
		min_height="840px",
		layout_sidebar(
			sidebar=sidebar(
				width=500,
				shinyDirButton('dir', label='Select input directory', title='Please select the directory containing the genomic prediction output', multiple=FALSE),
				shinyWidgets::pickerInput(inputId="within_vec_traits", label="Filter by trait:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="within_vec_populations", label="Filter by population:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="within_vec_models", label="Filter by model:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="within_metric", label="Use the genomic prediction accuracy metric:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="within_vec_group_by", label="Group by:", choices=c("trait", "pop_training", "model"), selected=c("trait", "model"), multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="within_sort_by", label="Sort by:", choices=c("increasing_mean", "decreasing_mean", "increasing_median", "decreasing_median", "alphabetical"), selected="increasing_mean", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::materialSwitch(inputId="within_bool_box_with_labels", label="Mean-labelled boxplot", value=FALSE, status="primary", right=FALSE)
			),
			mainPanel(
				width=750,
				shinycssloaders::withSpinner(plotlyOutput(outputId="within_plot_violin_or_box")),
				br(),
				br(),
				downloadButton("within_download_metrics", "Download summary and raw data tables previewed below"),
				br(),
				br(),
				div("Preview of summary statistics of the current selection:"),
				shinycssloaders::withSpinner(tableOutput(outputId="within_data_tables"))
			)
		)
	)
}

fn_within_barplot_ui = function() {
	card(
		card_header(h1(strong("Within population bar plot"), style="font-size:21px; text-align:left")),
		min_height="700px",
		layout_sidebar(
			sidebar=sidebar(
				width=500,
				shinyWidgets::pickerInput(inputId="within_barplot_vec_traits", label="Filter by trait:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="within_barplot_vec_populations", label="Filter by population:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="within_barplot_metric", label="Use the genomic prediction accuracy metric:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="within_barplot_model", label="Model:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::materialSwitch(inputId="within_barplot_bool_group_pop", label="Group by population", value=FALSE, status="primary", right=TRUE)
			),
			mainPanel(
				width=750,
				shinycssloaders::withSpinner(plotlyOutput(outputId="within_plot_bar"))
			)
		)
	)
}

fn_within_scatter_hist_ui = function() {
	card(
		card_header(h1(strong("Within population observed vs predicted phenotypes"), style="font-size:21px; text-align:left")),
		min_height="1200px",
		layout_sidebar(
			sidebar=sidebar(
				width=500,
				shinyWidgets::pickerInput(inputId="within_scat_trait_x_pop", label="Trait-Pop combination:", choices="", options=list(`live-search`=TRUE)),
				shinyWidgets::pickerInput(inputId="within_scat_model", label="Model:", choices="", options=list(`live-search`=TRUE)),
				shinyWidgets::materialSwitch(inputId="within_scat_points", label="Plot points", value=FALSE, status="primary", right=TRUE),
				sliderInput(inputId="within_bins", label="Number of bins:", min=1, max=100, value=10)
			),
			mainPanel(
				width=750,
				shinycssloaders::withSpinner(plotlyOutput(outputId="within_plot_scatter")),
				br(),
				br(),
				br(),
				br(),
				shinycssloaders::withSpinner(plotlyOutput(outputId="within_plot_histogram")),
				downloadButton("within_download_predictions", "Download observed and predicted phenotypes given the current selection above")
			)
		)
	)
}

fn_within_violin_server = function(input, list_list_output) {
	### Define the grouping, i.e. trait and/or pop and/or model
	METRICS_WITHIN_POP = NULL
	for (list_output in list_list_output) {
		df_metrics_within_pop = list_output$METRICS_WITHIN_POP
		if (is.na(df_metrics_within_pop[1])[1]) {next}
		df_metrics_within_pop$trait = list_output$TRAIT_NAME
		if (is.null(METRICS_WITHIN_POP)) {
			METRICS_WITHIN_POP = df_metrics_within_pop
		} else {
			METRICS_WITHIN_POP = rbind(METRICS_WITHIN_POP, df_metrics_within_pop)
		}
	}
	if (is.null(METRICS_WITHIN_POP)) {
		df = data.frame(ID="No within populations data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
		return(p)
	}
	vec_idx = which(
		(METRICS_WITHIN_POP$trait %in% input$within_vec_traits) &
		(METRICS_WITHIN_POP$pop_training %in% input$within_vec_populations) &
		(METRICS_WITHIN_POP$model %in% input$within_vec_models))
	grouping = c()
	x_label_one_level = c()
	x_label_one_level_names = c()
	x_label_multiple_levels = c()
	for (g in input$within_vec_group_by) {
		g_names = as.character(eval(parse(text=paste0("METRICS_WITHIN_POP$", g, "[vec_idx]"))))
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
	### Sort by mean, median or alphabetical (default)
	shiny::validate(shiny::need(grouping > 0, "Error: Please select a grouping factor ('Group by:' selection) with more than two levels selected."))
	df = data.frame(
		grouping=grouping,
		metric=eval(parse(text=paste0("METRICS_WITHIN_POP$", input$within_metric, "[vec_idx]")))
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
	### Plot
	if (!input$within_bool_box_with_labels) {
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
		yaxis=list(title=input$within_metric),
		xaxis=list(title="")
	)
	p = p %>% config(toImageButtonOptions = list(format = "svg"))
	return(p)
}

fn_within_table_metrics_server = function(input, list_list_output) {
	METRICS_WITHIN_POP = NULL
	for (list_output in list_list_output) {
		df_metrics_within_pop = list_output$METRICS_WITHIN_POP
		if (is.na(df_metrics_within_pop[1])[1]) {next}
		df_metrics_within_pop$trait = list_output$TRAIT_NAME
		if (is.null(METRICS_WITHIN_POP)) {
			METRICS_WITHIN_POP = df_metrics_within_pop
		} else {
			METRICS_WITHIN_POP = rbind(METRICS_WITHIN_POP, df_metrics_within_pop)
		}
	}
	if (is.null(METRICS_WITHIN_POP)) {
		return(list(df_stats=NULL, df_raw=NULL))
	}
	vec_idx = which((METRICS_WITHIN_POP$trait %in% input$within_vec_traits) &
		(METRICS_WITHIN_POP$pop_validation %in% input$within_vec_populations) &
		(METRICS_WITHIN_POP$model %in% input$within_vec_models))
	df_raw = METRICS_WITHIN_POP[vec_idx, , drop=FALSE]

	df_agg_mean = eval(parse(text=paste0("aggregate(", input$within_metric, " ~ ", paste(input$within_vec_group_by, collapse="+"), ", data=df_raw, FUN=mean, na.rm=TRUE)"))); colnames(df_agg_mean)[ncol(df_agg_mean)] = paste0("mean_", input$within_metric)
	df_agg_median = eval(parse(text=paste0("aggregate(", input$within_metric, " ~ ", paste(input$within_vec_group_by, collapse="+"), ", data=df_raw, FUN=median, na.rm=TRUE)"))); colnames(df_agg_median)[ncol(df_agg_median)] = paste0("median_", input$within_metric)
	df_agg_min = eval(parse(text=paste0("aggregate(", input$within_metric, " ~ ", paste(input$within_vec_group_by, collapse="+"), ", data=df_raw, FUN=min, na.rm=TRUE)"))); colnames(df_agg_min)[ncol(df_agg_min)] = paste0("min_", input$within_metric)
	df_agg_max = eval(parse(text=paste0("aggregate(", input$within_metric, " ~ ", paste(input$within_vec_group_by, collapse="+"), ", data=df_raw, FUN=max, na.rm=TRUE)"))); colnames(df_agg_max)[ncol(df_agg_max)] = paste0("max_", input$within_metric)
	df_agg_var = eval(parse(text=paste0("aggregate(", input$within_metric, " ~ ", paste(input$within_vec_group_by, collapse="+"), ", data=df_raw, FUN=var, na.rm=TRUE)"))); colnames(df_agg_var)[ncol(df_agg_var)] = paste0("var_", input$within_metric)
	df_agg_nan = eval(parse(text=paste0("aggregate(", input$within_metric, " ~ ", paste(input$within_vec_group_by, collapse="+"), ", data=df_raw, FUN=function(x){sum(is.na(x))})"))); colnames(df_agg_nan)[ncol(df_agg_nan)] = paste0("n_missing_", input$within_metric)
	df_stats = merge(df_agg_mean, 
		merge(df_agg_median, 
		merge(df_agg_min, 
		merge(df_agg_max, 
		merge(df_agg_var, df_agg_nan, 
			by=input$within_vec_group_by),
			by=input$within_vec_group_by),
			by=input$within_vec_group_by),
			by=input$within_vec_group_by),
			by=input$within_vec_group_by)
	return(list(df_stats=df_stats, df_raw=df_raw))
}

fn_within_df_stats_for_barplot_server = function(input, list_list_output) {
	# list_list_output = fn_io_server(dir="/group/pasture/Jeff/lucerne/workdir/gs/output_ground_truth_biomass_traits/output")
	# names(list_list_output)
	# input = list(
	# 	within_barplot_vec_traits=unique(unlist(lapply(list_list_output, FUN=function(x){x$TRAIT_NAME}))), 
	# 	within_barplot_vec_populations=unique(unlist(lapply(list_list_output, FUN=function(x){x$POPULATION})))

	METRICS_WITHIN_POP = NULL
	for (list_output in list_list_output) {
		df_metrics_within_pop = list_output$METRICS_WITHIN_POP
		if (is.na(df_metrics_within_pop[1])[1]) {next}
		df_metrics_within_pop$trait = list_output$TRAIT_NAME
		if (is.null(METRICS_WITHIN_POP)) {
			METRICS_WITHIN_POP = df_metrics_within_pop
		} else {
			METRICS_WITHIN_POP = rbind(METRICS_WITHIN_POP, df_metrics_within_pop)
		}
	}
	if (is.null(METRICS_WITHIN_POP)) {
		return(list(df_stats=NULL, df_raw=NULL))
	}
	vec_idx = which((METRICS_WITHIN_POP$trait %in% input$within_barplot_vec_traits) &
		(METRICS_WITHIN_POP$pop_validation %in% input$within_barplot_vec_populations))
	df_raw = METRICS_WITHIN_POP[vec_idx, , drop=FALSE]

	vec_within_barplot_vec_group_by = c("trait", "pop_training", "model")
	df_agg_mean = eval(parse(text=paste0("aggregate(", input$within_barplot_metric, " ~ ", paste(vec_within_barplot_vec_group_by, collapse="+"), ", data=df_raw, FUN=mean, na.rm=TRUE)"))); colnames(df_agg_mean)[ncol(df_agg_mean)] = paste0("mean_", input$within_barplot_metric)
	df_agg_median = eval(parse(text=paste0("aggregate(", input$within_barplot_metric, " ~ ", paste(vec_within_barplot_vec_group_by, collapse="+"), ", data=df_raw, FUN=median, na.rm=TRUE)"))); colnames(df_agg_median)[ncol(df_agg_median)] = paste0("median_", input$within_barplot_metric)
	df_agg_min = eval(parse(text=paste0("aggregate(", input$within_barplot_metric, " ~ ", paste(vec_within_barplot_vec_group_by, collapse="+"), ", data=df_raw, FUN=min, na.rm=TRUE)"))); colnames(df_agg_min)[ncol(df_agg_min)] = paste0("min_", input$within_barplot_metric)
	df_agg_max = eval(parse(text=paste0("aggregate(", input$within_barplot_metric, " ~ ", paste(vec_within_barplot_vec_group_by, collapse="+"), ", data=df_raw, FUN=max, na.rm=TRUE)"))); colnames(df_agg_max)[ncol(df_agg_max)] = paste0("max_", input$within_barplot_metric)
	df_agg_var = eval(parse(text=paste0("aggregate(", input$within_barplot_metric, " ~ ", paste(vec_within_barplot_vec_group_by, collapse="+"), ", data=df_raw, FUN=var, na.rm=TRUE)"))); colnames(df_agg_var)[ncol(df_agg_var)] = paste0("var_", input$within_barplot_metric)
	df_agg_nan = eval(parse(text=paste0("aggregate(", input$within_barplot_metric, " ~ ", paste(vec_within_barplot_vec_group_by, collapse="+"), ", data=df_raw, FUN=function(x){sum(is.na(x))})"))); colnames(df_agg_nan)[ncol(df_agg_nan)] = paste0("n_missing_", input$within_barplot_metric)
	df_stats = merge(df_agg_mean, 
		merge(df_agg_median, 
		merge(df_agg_min, 
		merge(df_agg_max, 
		merge(df_agg_var, df_agg_nan, 
			by=vec_within_barplot_vec_group_by),
			by=vec_within_barplot_vec_group_by),
			by=vec_within_barplot_vec_group_by),
			by=vec_within_barplot_vec_group_by),
			by=vec_within_barplot_vec_group_by)
	df_stats = df_stats[order(df_stats$trait, decreasing=FALSE), , drop=FALSE]
	return(df_stats)
}

fn_within_barplot_server = function(input, list_list_output) {
	df_stats = fn_within_df_stats_for_barplot_server(input=input, list_list_output=list_list_output)
	df_stats = df_stats[df_stats$model %in% input$within_barplot_model, , drop=FALSE]
	if (!input$within_barplot_bool_group_pop) {
		p = eval(parse(text=paste0("plotly::plot_ly(data=df_stats, x=~trait, y=~mean_", input$within_barplot_metric, ", type='bar', name=~pop_training, error_y=~list(array=sqrt(var_", input$within_barplot_metric, "), color='#000000'))")))
	} else {
		p = eval(parse(text=paste0("plotly::plot_ly(data=df_stats, x=~pop_training, y=~mean_", input$within_barplot_metric, ", type='bar', name=~trait, error_y=~list(array=sqrt(var_", input$within_barplot_metric, "), color='#000000'))")))
	}
	p = p %>% config(toImageButtonOptions = list(format = "svg"))
	return(p)


}

fn_within_scatterplot_server = function(input, list_list_output) {
	list_output = list_list_output[[input$within_scat_trait_x_pop]]
	if (is.null(list_output)) {
		df = data.frame(ID="No within populations data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
		return(p)
	}
    idx = which(list_output$YPRED_WITHIN_POP$model == input$within_scat_model)
    df_tmp = data.frame(ID=list_output$YPRED_WITHIN_POP$id[idx], Observed=list_output$YPRED_WITHIN_POP$y_true[idx], Predicted=list_output$YPRED_WITHIN_POP$y_pred[idx])
    if (input$within_scat_points) {
      df = df_tmp
      df$Predicted_sd = 0
    } else {
      df_observed = aggregate(Observed ~ ID, data=df_tmp, FUN=mean, na.rm=TRUE); colnames(df_observed) = c("ID", "Observed")
      df_predicted = aggregate(Predicted ~ ID, data=df_tmp, FUN=mean, na.rm=TRUE); colnames(df_predicted) = c("ID", "Predicted")
      df_predicted_sd = aggregate(Predicted ~ ID, data=df_tmp, FUN=sd, na.rm=TRUE); colnames(df_predicted_sd) = c("ID", "Predicted_sd")
      df = merge(merge(df_observed, df_predicted, by="ID", all.x=TRUE), df_predicted_sd, by="ID", all.x=TRUE)
    }
    fit = lm(Predicted ~ Observed, data=df); new_x=seq(from=min(df$Observed,na.rm=TRUE), to=max(df$Observed,na.rm=TRUE), length=1000)
    p = plotly::plot_ly(data=df, x=~Observed, y=~Predicted, type="scatter", mode='markers', text=~ID, error_y=~list(array=Predicted_sd, thickness=0.5, color="gray"))
    p = p %>% add_trace(x=new_x, y=predict(fit, newdata=data.frame(Observed=new_x)),
      mode='lines', 
      hoverinfo='text',
      text=~paste(
        "Slope:", round(coef(fit)[2], 4), 
        "<br>Intercept:", round(coef(fit)[1], 4),
        "<br>R-squared:", round(100*summary(fit)$r.squared), "%",
        "<br>Correaltion:", round(100*cor(df$Observed, df$Predicted)), "%",
        "<br>n = ", nrow(df),
        "<br>folds = ", length(unique(list_output$METRICS_WITHIN_POP$k))),
      error_y=NULL
    )
    p = p %>% layout(
      title=paste0(
        "Trait: ", list_output$TRAIT_NAME, 
        " | Population: ", unique(list_output$METRICS_WITHIN_POP$pop_validation), 
        " | Model: ", input$within_scat_model
      ),
      showlegend=FALSE
    )
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
	return(p)
}

fn_within_histogram_server = function(input, list_list_output) {
    list_output = list_list_output[[input$within_scat_trait_x_pop]]
	if (is.null(list_output)) {
		df = data.frame(ID="No within populations data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
		return(p)
	}
    idx = which(list_output$YPRED_WITHIN_POP$model == input$within_scat_model)
    df = data.frame(ID=list_output$YPRED_WITHIN_POP$id[idx], Observed=list_output$YPRED_WITHIN_POP$y_true[idx], Predicted=list_output$YPRED_WITHIN_POP$y_pred[idx])
    p = plot_ly(alpha=0.5, type="histogram", nbinsx=input$within_bins)
    p = p %>% add_histogram(data=df, x=~Observed, name="Observed")
    p = p %>% add_histogram(data=df, x=~Predicted, name="Predicted")
    p = p %>% layout(
      title="Phenotype distribution",
      barmode="overlay", 
      xaxis=list(title=list_output$TRAIT_NAME)
    )
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
	return(p)
}

fn_within_table_predictions_server = function(input, list_list_output) {
	list_output = list_list_output[[input$within_scat_trait_x_pop]]
	if (is.null(list_output)) {
		return(NULL)
	}
    vec_idx = which(list_output$YPRED_WITHIN_POP$model == input$within_scat_model)
    df = data.frame(
		pop_training=list_output$YPRED_WITHIN_POP$pop_training[vec_idx],
		pop_validation=list_output$YPRED_WITHIN_POP$pop_validation[vec_idx],
		model=list_output$YPRED_WITHIN_POP$model[vec_idx],
		id=list_output$YPRED_WITHIN_POP$id[vec_idx],
		observed=list_output$YPRED_WITHIN_POP$y_true[vec_idx],
		predicted=list_output$YPRED_WITHIN_POP$y_pred[vec_idx])
	return(df)
}

###################################
### ACROSS POPULATIONS PAIRWISE ###
###################################
fn_across_pair_barplot_ui = function() {
	card(
		card_header(h1(strong("Across populations: pairwise-population prediction accuracies"), style="font-size:21px; text-align:left")),
		min_height="750px",
		layout_sidebar(
			sidebar=sidebar(
				width=500,
				shinyWidgets::pickerInput(inputId="across_pair_trait", label="Filter by trait:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_pair_metric", label="Use the genomic prediction accuracy metric:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::materialSwitch(inputId="across_pair_bool_group_model", label="Group by model", value=FALSE, status="primary", right=TRUE),
				shinyWidgets::pickerInput(inputId="across_pair_pop_training", label="Training population:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_pair_vec_pop_validation", label="Validation population/s:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_pair_vec_models", label="Filter by model:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE))
			),
			mainPanel(
				width=750,
				shinycssloaders::withSpinner(plotlyOutput(outputId="across_pair_plot_bar")),
				downloadButton("across_pair_download_metrics", "Download subset and full data tables previewed below"),
				br(),
				br(),
				div("Preview of data filtered by current selection:"),
				shinycssloaders::withSpinner(tableOutput(outputId="across_pair_data_tables"))
			)
		)
	)
}

fn_across_pair_scatter_hist_ui = function() {
	card(
		card_header(h1(strong("Across populations: pairwise-population observed vs predicted phenotypes"), style="font-size:21px; text-align:left")),
		min_height="1200px",
		layout_sidebar(
			sidebar=sidebar(
				width=500,
				shinyWidgets::pickerInput(inputId="across_pair_scat_trait", label="Trait (Note: Ignore the population label as the first population contains the PAIRWISE output):", choices="", options=list(`live-search`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_pair_scat_pop_training", label="Training population:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_pair_scat_pop_validation", label="Validation population:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_pair_scat_model", label="Model:", choices="", options=list(`live-search`=TRUE)),
				sliderInput(inputId="across_pair_bins", label="Number of bins:", min=1, max=100, value=10)
			),
			mainPanel(
				width=750,
				shinycssloaders::withSpinner(plotlyOutput(outputId="across_pair_plot_scatter")),
				br(),
				br(),
				br(),
				br(),
				shinycssloaders::withSpinner(plotlyOutput(outputId="across_pair_plot_histogram")),
				downloadButton("across_pair_download_predictions", "Download observed and predicted phenotypes given the current selection above")
			)
		)
	)
}

fn_across_pair_barplot_server = function(input, list_list_output) {
	### Extract the data frame of prediction metrics from pairwise-population cross-validation
	df_metrics_across_pop = NULL
	for (list_output in list_list_output) {
		if (list_output$TRAIT_NAME != input$across_pair_trait) {
			next
		}
		if (is.na(list_output$METRICS_ACROSS_POP_PAIRWISE[1])[1]) {
			next
		}
		if (is.null(df_metrics_across_pop)) {
			df_metrics_across_pop = list_output$METRICS_ACROSS_POP_PAIRWISE
		} else {
			df_metrics_across_pop = rbind(df_metrics_across_pop, list_output$METRICS_ACROSS_POP_PAIRWISE)
		}
	}
	if (is.null(df_metrics_across_pop)) {
		df = data.frame(ID="No across populations data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
		return(p)
	}
	### Retain only the data set trained with the input population and validated in the defined validation populations
	vec_pairwise_pop_validation = input$across_pair_vec_pop_validation[!(input$across_pair_vec_pop_validation %in% input$across_pair_pop_training)]
	vec_idx = which(
		(df_metrics_across_pop$pop_training == input$across_pair_pop_training) &
		(df_metrics_across_pop$pop_validation %in% vec_pairwise_pop_validation)
	)
	df_metrics_across_pop = df_metrics_across_pop[vec_idx, ]
	vec_pairwise_pop_validation = sort(unique(df_metrics_across_pop$pop_validation))
	if (nrow(df_metrics_across_pop) == 0) {
		df = data.frame(ID="No leave-one-population-out data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
	} else {
		### Prepare 1 training vs multiple validation populations vs multiple models matrix of GP metric data
		n = length(vec_pairwise_pop_validation)
		m = length(input$across_pair_vec_models)
		M = matrix(NA, nrow=n, ncol=m)
		rownames(M) = vec_pairwise_pop_validation
		colnames(M) = input$across_pair_vec_models
		for (i in 1:n) {
			for (j in 1:m) {
				vec_idx = which((df_metrics_across_pop$pop_validation == vec_pairwise_pop_validation[i]) & (df_metrics_across_pop$model == input$across_pair_vec_models[j]))
				M[i, j] = eval(parse(text=paste0("df_metrics_across_pop$", input$across_pair_metric, "[vec_idx]")))
			}
		}
		vec_validation_pops = sort(unique(as.character(df_metrics_across_pop$pop_validation)))
		vec_models = sort(unique(as.character(df_metrics_across_pop$model)))
		if (!input$across_pair_bool_group_model) {
			df = data.frame(vec_pairwise_pop_validation, M)
			colnames(df) = c("pop_validation", colnames(M))
			### Plot
			eval(parse(text=paste0("p = plotly::plot_ly(data=df, x=~pop_validation, y=~", input$across_pair_vec_models[1], ", type='bar', name='", input$across_pair_vec_models[1], "')")))
			vec_metric = eval(parse(text=paste0("df$", input$across_pair_vec_models[1])))
			if (length(input$across_pair_vec_models) > 1) {
				for (i in 2:length(input$across_pair_vec_models)) {
					eval(parse(text=paste0("p = p %>% plotly::add_trace(y=~", input$across_pair_vec_models[i], ", name='", input$across_pair_vec_models[i], "')")))
					vec_metric = c(vec_metric, eval(parse(text=paste0("df$", input$across_pair_vec_models[i]))))
				}
			}
		} else {
			df = data.frame(input$across_pair_vec_models, t(M))
			colnames(df) = c("model", gsub("-" , ".", rownames(M)))
			### Plot
			eval(parse(text=paste0("p = plotly::plot_ly(data=df, x=~model, y=~", gsub("-", ".", vec_pairwise_pop_validation)[1], ", type='bar', name='", gsub("-", ".", vec_pairwise_pop_validation)[1], "')")))
			vec_metric = eval(parse(text=paste0("df$", gsub("-", ".", vec_pairwise_pop_validation)[1])))
			if (length(gsub("-", ".", vec_pairwise_pop_validation)) > 1) {
				for (i in 2:length(gsub("-", ".", vec_pairwise_pop_validation))) {
					eval(parse(text=paste0("p = p %>% plotly::add_trace(y=~", gsub("-", ".", vec_pairwise_pop_validation)[i], ", name='", gsub("-", ".", vec_pairwise_pop_validation)[i], "')")))
					vec_metric = c(vec_metric, eval(parse(text=paste0("df$", gsub("-", ".", vec_pairwise_pop_validation)[i]))))
				}
			}
		}
		p = p %>% plotly::layout(
			title=paste0("Trait: ", input$across_pair_trait, 
				" | Trained on ", input$across_pair_pop_training,
				" | Mean ", input$across_pair_metric, " = ", round(mean(vec_metric, na.rm=TRUE), 4)), 
			barmode="group",
			yaxis=list(title=input$across_pair_metric)
		)
	}
	p = p %>% config(toImageButtonOptions = list(format = "svg"))
	return(p)
}

fn_across_pair_table_metrics_server = function(input, list_list_output) {
	METRICS_ACROSS_POP_PAIRWISE = NULL
	for (list_output in list_list_output) {
		df_metrics_across_pop_pairwise = list_output$METRICS_ACROSS_POP_PAIRWISE
		if (is.na(df_metrics_across_pop_pairwise[1])[1]) {next}
		df_metrics_across_pop_pairwise = cbind(data.frame(trait=rep(list_output$TRAIT_NAME, nrow(df_metrics_across_pop_pairwise))), df_metrics_across_pop_pairwise)
		# df_metrics_across_pop_pairwise$trait = list_output$TRAIT_NAME
		if (is.null(METRICS_ACROSS_POP_PAIRWISE)) {
			METRICS_ACROSS_POP_PAIRWISE = df_metrics_across_pop_pairwise
		} else {
			METRICS_ACROSS_POP_PAIRWISE = rbind(METRICS_ACROSS_POP_PAIRWISE, df_metrics_across_pop_pairwise)
		}
	}
	vec_idx_row = which((METRICS_ACROSS_POP_PAIRWISE$trait %in% input$across_pair_trait) &
		(METRICS_ACROSS_POP_PAIRWISE$pop_training %in% input$across_pair_pop_training) &
		(METRICS_ACROSS_POP_PAIRWISE$pop_validation %in% input$across_pair_vec_pop_validation) &
		(METRICS_ACROSS_POP_PAIRWISE$model %in% input$across_pair_vec_models))
	vec_idx_col = which(colnames(METRICS_ACROSS_POP_PAIRWISE) %in% c(input$across_pair_metric, "trait", "pop_training", "pop_validation", "model"))
	return(list(
		df_sub=METRICS_ACROSS_POP_PAIRWISE[vec_idx_row, vec_idx_col, drop=FALSE],
		df_full=METRICS_ACROSS_POP_PAIRWISE
	))
}

fn_across_pair_scatter_server = function(input, list_list_output) {
	list_output = list_list_output[[input$across_pair_scat_trait]]
	# if (!is.na(list_output$YPRED_ACROSS_POP_PAIRWISE[1])[1]) {
	if (is.null(list_output$YPRED_ACROSS_POP_PAIRWISE)) {
		df = data.frame(ID="No across populations data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
		return(p)
	}
	vec_idx = which(
		(list_output$YPRED_ACROSS_POP_PAIRWISE$model == input$across_pair_scat_model) & 
		(list_output$YPRED_ACROSS_POP_PAIRWISE$pop_training == input$across_pair_scat_pop_training) &
		(list_output$YPRED_ACROSS_POP_PAIRWISE$pop_validation == input$across_pair_scat_pop_validation))
	POP = paste0("trained: ", list_output$YPRED_ACROSS_POP_PAIRWISE$pop_training[vec_idx], "\nvalidated: ", list_output$YPRED_ACROSS_POP_PAIRWISE$pop_validation[vec_idx])
	df = data.frame(ID=list_output$YPRED_ACROSS_POP_PAIRWISE$id[vec_idx], POP=POP, Observed=list_output$YPRED_ACROSS_POP_PAIRWISE$y_true[vec_idx], Predicted=list_output$YPRED_ACROSS_POP_PAIRWISE$y_pred[vec_idx])
	fit = lm(Predicted ~ Observed, data=df); new_x=seq(from=min(df$Observed,na.rm=TRUE), to=max(df$Observed,na.rm=TRUE), length=1000)
	p = plotly::plot_ly(data=df, x=~Observed, y=~Predicted, type="scatter", mode='markers', text=~ID, color=~POP, showlegend=TRUE)
	p = p %>% add_trace(x=new_x, y=predict(fit, newdata=data.frame(Observed=new_x)),
		color=NULL,
		name="linear fit",
		mode='lines', 
		hoverinfo='text',
		text=~paste(
			"Slope:", round(coef(fit)[2], 4), 
			"<br>Intercept:", round(coef(fit)[1], 4),
			"<br>R-squared:", round(100*summary(fit)$r.squared), "%",
			"<br>Correaltion:", round(100*cor(df$Observed, df$Predicted)), "%",
			"<br>n = ", nrow(df),
			"<br>folds = ", length(unique(list_output$METRICS_WITHIN_POP$k))),
		error_y=NULL
	)
	p = p %>% layout(
		title=paste0(
			"Trait: ", list_output$TRAIT_NAME, 
			" | Model: ", input$across_pair_scat_model
		)
	)
	p = p %>% config(toImageButtonOptions = list(format = "svg"))
	return(p)
}

fn_across_pair_histogram_server = function(input, list_list_output) {
    list_output = list_list_output[[input$across_pair_scat_trait]]
    # if (!is.na(list_output$YPRED_ACROSS_POP_PAIRWISE[1])[1]) {
    if (is.null(list_output$YPRED_ACROSS_POP_PAIRWISE)) {
		df = data.frame(ID="No across populations data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
		return(p)
	}
	vec_idx = which(
	(list_output$YPRED_ACROSS_POP_PAIRWISE$model == input$across_pair_scat_model) & 
	(list_output$YPRED_ACROSS_POP_PAIRWISE$pop_training == input$across_pair_scat_pop_training) &
	(list_output$YPRED_ACROSS_POP_PAIRWISE$pop_validation == input$across_pair_scat_pop_validation))
	df = data.frame(ID=list_output$YPRED_ACROSS_POP_PAIRWISE$id[vec_idx], Observed=list_output$YPRED_ACROSS_POP_PAIRWISE$y_true[vec_idx], Predicted=list_output$YPRED_ACROSS_POP_PAIRWISE$y_pred[vec_idx])
	p = plot_ly(alpha=0.5, type="histogram", nbinsx=input$across_pair_bins)
	p = p %>% add_histogram(data=df, x=~Observed, name="Observed")
	p = p %>% add_histogram(data=df, x=~Predicted, name="Predicted")
	p = p %>% layout(
	title="Phenotype distribution",
	barmode="overlay", 
	xaxis=list(title=list_output$TRAIT_NAME)
	)
	p = p %>% config(toImageButtonOptions = list(format = "svg"))
	return(p)
}

fn_across_pair_table_predictions_server = function(input, list_list_output) {
	list_output = list_list_output[[input$across_pair_scat_trait]]
	vec_idx = which(
		(list_output$YPRED_ACROSS_POP_PAIRWISE$model == input$across_pair_scat_model) & 
		(list_output$YPRED_ACROSS_POP_PAIRWISE$pop_training == input$across_pair_scat_pop_training) &
		(list_output$YPRED_ACROSS_POP_PAIRWISE$pop_validation == input$across_pair_scat_pop_validation))
	df = data.frame(
		pop_training=input$across_pair_scat_pop_training,
		pop_validation=input$across_pair_scat_pop_validation,
		model=input$across_pair_scat_model,
		id=list_output$YPRED_ACROSS_POP_PAIRWISE$id[vec_idx],
		observed=list_output$YPRED_ACROSS_POP_PAIRWISE$y_true[vec_idx],
		predicted=list_output$YPRED_ACROSS_POP_PAIRWISE$y_pred[vec_idx])
	return(df)
}

##########################################################
### ACROSS POPULATIONS LEAVE-ONE-POPULATION-OUT (LOPO) ###
##########################################################
fn_across_lopo_barplot_ui = function() {
	card(
		card_header(h1(strong("Across populations: leave-one-population-out prediction accuracies"), style="font-size:21px; text-align:left")),
		min_height="750px",
		layout_sidebar(
			sidebar=sidebar(
				width=500,
				shinyWidgets::pickerInput(inputId="across_lopo_trait", label="Filter by trait:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_lopo_metric", label="Use the genomic prediction accuracy metric:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_lopo_pop_validation", label="Validation populations:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_lopo_vec_models", label="Filter by model:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE))
			),
			mainPanel(
				width=750,
				shinycssloaders::withSpinner(plotlyOutput(outputId="across_lopo_plot_bar")),
				downloadButton("across_lopo_download_metrics", "Download subset and full data tables previewed below"),
				br(),
				br(),
				div("Preview of data filtered by current selection:"),
				shinycssloaders::withSpinner(tableOutput(outputId="across_lopo_data_tables"))
			)
		)
	)
}

fn_across_lopo_scatter_hist_ui = function() {
	card(
		card_header(h1(strong("Across populations: leave-one-population-out observed vs predicted phenotypes"), style="font-size:21px; text-align:left")),
		min_height="1200px",
		layout_sidebar(
			sidebar=sidebar(
				width=500,
				shinyWidgets::pickerInput(inputId="across_lopo_scat_trait", label="Trait (Note: Ignore the population label as the first population contains the LOPO output):", choices="", options=list(`live-search`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_lopo_scat_pop_validation", label="Validation population:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_lopo_scat_model", label="Model:", choices="", options=list(`live-search`=TRUE)),
				sliderInput(inputId="across_lopo_bins", label="Number of bins:", min=1, max=100, value=10)
			),
			mainPanel(
				width=750,
				shinycssloaders::withSpinner(plotlyOutput(outputId="across_lopo_plot_scatter")),
				br(),
				br(),
				br(),
				br(),
				shinycssloaders::withSpinner(plotlyOutput(outputId="across_lopo_plot_histogram")),
				downloadButton("across_lopo_download_predictions", "Download observed and predicted phenotypes given the current selection above")
			)
		)
	)
}

fn_across_lopo_barplot_server = function(input, list_list_output) {
	### Extract the data frame of prediction metrics from leave-one-population-out cross-validation
	df_metrics_across_pop = NULL
	for (list_output in list_list_output) {
		if (list_output$TRAIT_NAME != input$across_lopo_trait) {
			next
		}
		if (is.na(list_output$METRICS_ACROSS_POP_LOPO[1])[1]) {
			next
		}
		if (is.null(df_metrics_across_pop)) {
			df_metrics_across_pop = list_output$METRICS_ACROSS_POP_LOPO
		} else {
			df_metrics_across_pop = rbind(df_metrics_across_pop, list_output$METRICS_ACROSS_POP_LOPO)
		}
	}
	### Revise to account for both training and validation populations
	if (is.null(df_metrics_across_pop)) {
		df = data.frame(ID="No across populations data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
		return(p)
	}
	### Retain only the data set trained with the input population and validated in the defined validation populations
	vec_idx = which(df_metrics_across_pop$pop_validation == input$across_lopo_pop_validation)
	df_metrics_across_pop = df_metrics_across_pop[vec_idx, ]
	vec_pop_validation = sort(unique(df_metrics_across_pop$pop_validation))
	if (nrow(df_metrics_across_pop) == 0) {
		df = data.frame(ID="No leave-one-population-out data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
	} else {
		### Prepare 1 training vs multiple validation populations vs multiple models matrix of GP metric data
		n = length(vec_pop_validation)
		m = length(input$across_lopo_vec_models)
		M = matrix(NA, nrow=n, ncol=m)
		rownames(M) = vec_pop_validation
		colnames(M) = input$across_lopo_vec_models
		for (i in 1:n) {
			for (j in 1:m) {
				vec_idx = which((df_metrics_across_pop$pop_validation == vec_pop_validation[i]) & (df_metrics_across_pop$model == input$across_lopo_vec_models[j]))
				M[i, j] = eval(parse(text=paste0("df_metrics_across_pop$", input$across_lopo_metric, "[vec_idx]")))
			}
		}
		vec_validation_pops = sort(unique(as.character(df_metrics_across_pop$pop_validation)))
		vec_models = sort(unique(as.character(df_metrics_across_pop$model)))
		df = data.frame(vec_pop_validation, M)
		colnames(df) = c("pop_validation", colnames(M))
		### Plot
		eval(parse(text=paste0("p = plotly::plot_ly(data=df, x=~pop_validation, y=~", input$across_lopo_vec_models[1], ", type='bar', name='", input$across_lopo_vec_models[1], "')")))
		vec_metric = eval(parse(text=paste0("df$", input$across_lopo_vec_models[1])))
		if (length(input$across_lopo_vec_models) > 1) {
			for (i in 2:length(input$across_lopo_vec_models)) {
				eval(parse(text=paste0("p = p %>% plotly::add_trace(y=~", input$across_lopo_vec_models[i], ", name='", input$across_lopo_vec_models[i], "')")))
				vec_metric = c(vec_metric, eval(parse(text=paste0("df$", input$across_lopo_vec_models[i]))))
			}
		}
		p = p %>% plotly::layout(
			title=paste0("Trait: ", input$across_lopo_trait, 
				" | Trained on ", unique(df_metrics_across_pop$pop_training),
				" | Mean ", input$across_lopo_metric, " = ", round(mean(vec_metric, na.rm=TRUE), 4)), 
			barmode="group",
			yaxis=list(title=input$across_lopo_metric)
		)
	}
	p = p %>% config(toImageButtonOptions = list(format = "svg"))
	return(p)
}

fn_across_lopo_table_metrics_server = function(input, list_list_output) {
	METRICS_ACROSS_POP_LOPO = NULL
	for (list_output in list_list_output) {
		df_metrics_across_pop_lopo = list_output$METRICS_ACROSS_POP_LOPO
		if (is.na(df_metrics_across_pop_lopo[1])[1]) {next}
		df_metrics_across_pop_lopo = cbind(data.frame(trait=rep(list_output$TRAIT_NAME, nrow(df_metrics_across_pop_lopo))), df_metrics_across_pop_lopo)
		# df_metrics_across_pop_lopo$trait = list_output$TRAIT_NAME
		if (is.null(METRICS_ACROSS_POP_LOPO)) {
			METRICS_ACROSS_POP_LOPO = df_metrics_across_pop_lopo
		} else {
			METRICS_ACROSS_POP_LOPO = rbind(METRICS_ACROSS_POP_LOPO, df_metrics_across_pop_lopo)
		}
	}
	vec_idx_row = which((METRICS_ACROSS_POP_LOPO$trait %in% input$across_lopo_trait) &
		(METRICS_ACROSS_POP_LOPO$pop_validation %in% input$across_lopo_pop_validation) &
		(METRICS_ACROSS_POP_LOPO$model %in% input$across_lopo_vec_models))
	vec_idx_col = which(colnames(METRICS_ACROSS_POP_LOPO) %in% c(input$across_lopo_metric, "trait", "pop_training", "pop_validation", "model"))
	return(list(
		df_sub=METRICS_ACROSS_POP_LOPO[vec_idx_row, vec_idx_col, drop=FALSE],
		df_full=METRICS_ACROSS_POP_LOPO
	))
}

fn_across_lopo_scatter_server = function(input, list_list_output) {
	list_output = list_list_output[[input$across_lopo_scat_trait]]
	# if (!is.na(list_output$YPRED_ACROSS_POP_LOPO[1])[1]) {
	if (is.null(list_output$YPRED_ACROSS_POP_LOPO) | is.na(head(list_output$YPRED_ACROSS_POP_LOPO, n=1)[1])) {
		df = data.frame(ID="No across populations data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
		return(p)
	}
	vec_idx = which(
		(list_output$YPRED_ACROSS_POP_LOPO$model == input$across_lopo_scat_model) & 
		(list_output$YPRED_ACROSS_POP_LOPO$pop_validation == input$across_lopo_scat_pop_validation))
	POP = paste0("trained: ", list_output$YPRED_ACROSS_POP_LOPO$pop_training[vec_idx], "\nvalidated: ", list_output$YPRED_ACROSS_POP_LOPO$pop_validation[vec_idx])
	df = data.frame(ID=list_output$YPRED_ACROSS_POP_LOPO$id[vec_idx], POP=POP, Observed=list_output$YPRED_ACROSS_POP_LOPO$y_true[vec_idx], Predicted=list_output$YPRED_ACROSS_POP_LOPO$y_pred[vec_idx])
	fit = lm(Predicted ~ Observed, data=df); new_x=seq(from=min(df$Observed,na.rm=TRUE), to=max(df$Observed,na.rm=TRUE), length=1000)
	p = plotly::plot_ly(data=df, x=~Observed, y=~Predicted, type="scatter", mode='markers', text=~ID, color=~POP, showlegend=TRUE)
	p = p %>% add_trace(x=new_x, y=predict(fit, newdata=data.frame(Observed=new_x)),
		color=NULL,
		name="linear fit",
		mode='lines', 
		hoverinfo='text',
		text=~paste(
			"Slope:", round(coef(fit)[2], 4), 
			"<br>Intercept:", round(coef(fit)[1], 4),
			"<br>R-squared:", round(100*summary(fit)$r.squared), "%",
			"<br>Correaltion:", round(100*cor(df$Observed, df$Predicted)), "%",
			"<br>n = ", nrow(df),
			"<br>folds = ", length(unique(list_output$METRICS_WITHIN_POP$k))),
		error_y=NULL
	)
	p = p %>% layout(
		title=paste0(
			"Trait: ", list_output$TRAIT_NAME, 
			" | Model: ", input$across_lopo_scat_model
		)
	)
	p = p %>% config(toImageButtonOptions = list(format = "svg"))
	return(p)
}

fn_across_lopo_histogram_server = function(input, list_list_output) {
    list_output = list_list_output[[input$across_lopo_scat_trait]]
    # if (!is.na(list_output$YPRED_ACROSS_POP_LOPO[1])[1]) {
    if (!is.null(list_output$YPRED_ACROSS_POP_LOPO)) {
      vec_idx = which(
        (list_output$YPRED_ACROSS_POP_LOPO$model == input$across_lopo_scat_model) & 
        (list_output$YPRED_ACROSS_POP_LOPO$pop_validation == input$across_lopo_scat_pop_validation))
      df = data.frame(ID=list_output$YPRED_ACROSS_POP_LOPO$id[vec_idx], Observed=list_output$YPRED_ACROSS_POP_LOPO$y_true[vec_idx], Predicted=list_output$YPRED_ACROSS_POP_LOPO$y_pred[vec_idx])
	  p = plot_ly(alpha=0.5, type="histogram", nbinsx=input$across_lopo_bins)
		p = p %>% add_histogram(data=df, x=~Observed, name="Observed")
		p = p %>% add_histogram(data=df, x=~Predicted, name="Predicted")
		p = p %>% layout(
		title="Phenotype distribution",
		barmode="overlay", 
		xaxis=list(title=list_output$TRAIT_NAME)
		)
		p = p %>% config(toImageButtonOptions = list(format = "svg"))
    } else {
		df = data.frame(ID="No across populations data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
    }
}

fn_across_lopo_table_predictions_server = function(input, list_list_output) {
	list_output = list_list_output[[input$across_lopo_scat_trait]]
	vec_idx = which(
		(list_output$YPRED_ACROSS_POP_LOPO$model == input$across_lopo_scat_model) & 
		(list_output$YPRED_ACROSS_POP_LOPO$pop_validation == input$across_lopo_scat_pop_validation))
	df = data.frame(
		pop_training=list_output$YPRED_ACROSS_POP_LOPO$pop_training[vec_idx],
		pop_validation=input$across_lopo_scat_pop_validation,
		model=input$across_lopo_scat_model,
		id=list_output$YPRED_ACROSS_POP_LOPO$id[vec_idx],
		observed=list_output$YPRED_ACROSS_POP_LOPO$y_true[vec_idx],
		predicted=list_output$YPRED_ACROSS_POP_LOPO$y_pred[vec_idx])
	return(df)
}

###############################
### ACROSS POPULATIONS BULK ###
###############################
fn_across_bulk_violin_ui = function() {
	card(
		card_header(h1(strong("Across populations (bulk): prediction accuracies"), style="font-size:21px; text-align:left")),
		min_height="840px",
		layout_sidebar(
			sidebar=sidebar(
				width=500,
				shinyWidgets::pickerInput(inputId="across_bulk_vec_traits", label="Filter by trait:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_bulk_vec_models", label="Filter by model:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_bulk_metric", label="Use the genomic prediction accuracy metric:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_bulk_vec_group_by", label="Group by:", choices=c("trait", "model"), selected=c("trait", "pop", "model"), multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_bulk_sort_by", label="Sort by:", choices=c("increasing_mean", "decreasing_mean", "increasing_median", "decreasing_median", "alphabetical"), selected="increasing_mean", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
				shinyWidgets::materialSwitch(inputId="across_bulk_bool_box_with_labels", label="Mean-labelled boxplot", value=FALSE, status="primary", right=FALSE)
			),
			mainPanel(
				width=750,
				shinycssloaders::withSpinner(plotlyOutput(outputId="across_bulk_plot_violin_or_box")),
				br(),
				br(),
				downloadButton("across_bulk_download_metrics", "Download summary and raw data tables previewed below"),
				br(),
				br(),
				div("Preview of summary statistics of the current selection:"),
				shinycssloaders::withSpinner(tableOutput(outputId="across_bulk_data_tables"))
			)
		)
	)
}

fn_across_bulk_scatter_hist_ui = function() {
	card(
		card_header(h1(strong("Across populations (bulk): observed vs predicted phenotypes"), style="font-size:21px; text-align:left")),
		min_height="1200px",
		layout_sidebar(
			sidebar=sidebar(
				width=500,
				shinyWidgets::pickerInput(inputId="across_bulk_scat_trait_x_pop", label="Trait (Note: Ignore the population label as the first population contains the LOPO output):", choices="", options=list(`live-search`=TRUE)),
				shinyWidgets::pickerInput(inputId="across_bulk_scat_model", label="Model:", choices="", options=list(`live-search`=TRUE)),
				shinyWidgets::materialSwitch(inputId="across_bulk_scat_points", label="Plot points", value=FALSE, status="primary", right=TRUE),
				sliderInput(inputId="across_bulk_bins", label="Number of bins:", min=1, max=100, value=10)
			),
			mainPanel(
				width=750,
				shinycssloaders::withSpinner(plotlyOutput(outputId="across_bulk_plot_scatter")),
				br(),
				br(),
				br(),
				br(),
				shinycssloaders::withSpinner(plotlyOutput(outputId="across_bulk_plot_histogram")),
				downloadButton("across_bulk_download_predictions", "Download observed and predicted phenotypes given the current selection above")
			)
		)
	)
}

fn_across_bulk_violin_server = function(input, list_list_output) {
	### Define the grouping, i.e. trait and/or pop and/or model
	METRICS_ACROSS_POP_BULK = NULL
	for (list_output in list_list_output) {
		df_metrics_across_pop_bulk = list_output$METRICS_ACROSS_POP_BULK
		if (is.na(df_metrics_across_pop_bulk[1])[1]) {next}
		df_metrics_across_pop_bulk$trait = list_output$TRAIT_NAME
		if (is.null(METRICS_ACROSS_POP_BULK)) {
			METRICS_ACROSS_POP_BULK = df_metrics_across_pop_bulk
		} else {
			METRICS_ACROSS_POP_BULK = rbind(METRICS_ACROSS_POP_BULK, df_metrics_across_pop_bulk)
		}
	}
	if (is.null(METRICS_ACROSS_POP_BULK)) {
		df = data.frame(ID="No across populations data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
		return(p)
	}
	vec_idx = which(
		(METRICS_ACROSS_POP_BULK$trait %in% input$across_bulk_vec_traits) &
		(METRICS_ACROSS_POP_BULK$model %in% input$across_bulk_vec_models))
	grouping = c()
	x_label_one_level = c()
	x_label_one_level_names = c()
	x_label_multiple_levels = c()
	for (g in input$across_bulk_vec_group_by) {
		g_names = as.character(eval(parse(text=paste0("METRICS_ACROSS_POP_BULK$", g, "[vec_idx]"))))
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
	### Sort by mean, median or alphabetical (default)
	shiny::validate(shiny::need(grouping > 0, "Error: Please select a grouping factor ('Group by:' selection) with more than two levels selected."))
	df = data.frame(
		grouping=grouping,
		metric=eval(parse(text=paste0("METRICS_ACROSS_POP_BULK$", input$across_bulk_metric, "[vec_idx]")))
	)
	if (input$across_bulk_sort_by == "increasing_mean") {
		df_agg = aggregate(metric ~ grouping, data=df, FUN=mean, na.rm=TRUE)
		df$grouping = factor(df$grouping, levels=as.character(df_agg$grouping[order(df_agg$metric, decreasing=FALSE)]))
	}
	if (input$across_bulk_sort_by == "decreasing_mean") {
		df_agg = aggregate(metric ~ grouping, data=df, FUN=mean, na.rm=TRUE)
		df$grouping = factor(df$grouping, levels=as.character(df_agg$grouping[order(df_agg$metric, decreasing=TRUE)]))
	}
	if (input$across_bulk_sort_by == "increasing_median") {
		df_agg = aggregate(metric ~ grouping, data=df, FUN=median, na.rm=TRUE)
		df$grouping = factor(df$grouping, levels=as.character(df_agg$grouping[order(df_agg$metric, decreasing=FALSE)]))
	}
	if (input$across_bulk_sort_by == "decreasing_median") {
		df_agg = aggregate(metric ~ grouping, data=df, FUN=median, na.rm=TRUE)
		df$grouping = factor(df$grouping, levels=as.character(df_agg$grouping[order(df_agg$metric, decreasing=TRUE)]))
	}
	### Plot
	if (!input$across_bulk_bool_box_with_labels) {
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
		yaxis=list(title=input$across_bulk_metric),
		xaxis=list(title="")
	)
	p = p %>% config(toImageButtonOptions = list(format = "svg"))
	return(p)
}

fn_across_bulk_table_metrics_server = function(input, list_list_output) {
	METRICS_ACROSS_POP_BULK = NULL
	for (list_output in list_list_output) {
		df_metrics_across_pop_bulk = list_output$METRICS_ACROSS_POP_BULK
		if (is.na(df_metrics_across_pop_bulk[1])[1]) {next}
		df_metrics_across_pop_bulk$trait = list_output$TRAIT_NAME
		if (is.null(METRICS_ACROSS_POP_BULK)) {
			METRICS_ACROSS_POP_BULK = df_metrics_across_pop_bulk
		} else {
			METRICS_ACROSS_POP_BULK = rbind(METRICS_ACROSS_POP_BULK, df_metrics_across_pop_bulk)
		}
	}
	if (is.null(METRICS_ACROSS_POP_BULK)) {
			return(list(df_stats=NULL, df_raw=NULL))
	}
	vec_idx = which(
		(METRICS_ACROSS_POP_BULK$trait %in% input$across_bulk_vec_traits) &
		(METRICS_ACROSS_POP_BULK$model %in% input$across_bulk_vec_models))
	df_raw = METRICS_ACROSS_POP_BULK[vec_idx, , drop=FALSE]
	df_agg_mean = eval(parse(text=paste0("aggregate(", input$across_bulk_metric, " ~ ", paste(input$across_bulk_vec_group_by, collapse="+"), ", data=df_raw, FUN=mean, na.rm=TRUE)"))); colnames(df_agg_mean)[ncol(df_agg_mean)] = paste0("mean_", input$across_bulk_metric)
	df_agg_median = eval(parse(text=paste0("aggregate(", input$across_bulk_metric, " ~ ", paste(input$across_bulk_vec_group_by, collapse="+"), ", data=df_raw, FUN=median, na.rm=TRUE)"))); colnames(df_agg_median)[ncol(df_agg_median)] = paste0("median_", input$across_bulk_metric)
	df_agg_min = eval(parse(text=paste0("aggregate(", input$across_bulk_metric, " ~ ", paste(input$across_bulk_vec_group_by, collapse="+"), ", data=df_raw, FUN=min, na.rm=TRUE)"))); colnames(df_agg_min)[ncol(df_agg_min)] = paste0("min_", input$across_bulk_metric)
	df_agg_max = eval(parse(text=paste0("aggregate(", input$across_bulk_metric, " ~ ", paste(input$across_bulk_vec_group_by, collapse="+"), ", data=df_raw, FUN=max, na.rm=TRUE)"))); colnames(df_agg_max)[ncol(df_agg_max)] = paste0("max_", input$across_bulk_metric)
	df_agg_var = eval(parse(text=paste0("aggregate(", input$across_bulk_metric, " ~ ", paste(input$across_bulk_vec_group_by, collapse="+"), ", data=df_raw, FUN=var, na.rm=TRUE)"))); colnames(df_agg_var)[ncol(df_agg_var)] = paste0("var_", input$across_bulk_metric)
	df_agg_nan = eval(parse(text=paste0("aggregate(", input$across_bulk_metric, " ~ ", paste(input$across_bulk_vec_group_by, collapse="+"), ", data=df_raw, FUN=function(x){sum(is.na(x))})"))); colnames(df_agg_nan)[ncol(df_agg_nan)] = paste0("n_missing_", input$across_bulk_metric)
	df_stats = merge(df_agg_mean, 
		merge(df_agg_median, 
		merge(df_agg_min, 
		merge(df_agg_max, 
		merge(df_agg_var, df_agg_nan, 
			by=input$across_bulk_vec_group_by),
			by=input$across_bulk_vec_group_by),
			by=input$across_bulk_vec_group_by),
			by=input$across_bulk_vec_group_by),
			by=input$across_bulk_vec_group_by)
	return(list(df_stats=df_stats, df_raw=df_raw))
}

fn_across_bulk_scatterplot_server = function(input, list_list_output) {
	list_output = list_list_output[[input$across_bulk_scat_trait_x_pop]]
	if (is.null(list_output)) {
		df = data.frame(ID="No across populations data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
		return(p)
	}
    idx = which(list_output$YPRED_ACROSS_POP_BULK$model == input$across_bulk_scat_model)
    df_tmp = data.frame(ID=list_output$YPRED_ACROSS_POP_BULK$id[idx], Observed=list_output$YPRED_ACROSS_POP_BULK$y_true[idx], Predicted=list_output$YPRED_ACROSS_POP_BULK$y_pred[idx])
    if (input$across_bulk_scat_points) {
		df = df_tmp
		df$Predicted_sd = 0
    } else {
		df_observed = aggregate(Observed ~ ID, data=df_tmp, FUN=mean, na.rm=TRUE); colnames(df_observed) = c("ID", "Observed")
		df_predicted = aggregate(Predicted ~ ID, data=df_tmp, FUN=mean, na.rm=TRUE); colnames(df_predicted) = c("ID", "Predicted")
		df_predicted_sd = aggregate(Predicted ~ ID, data=df_tmp, FUN=sd, na.rm=TRUE); colnames(df_predicted_sd) = c("ID", "Predicted_sd")
		df = merge(merge(df_observed, df_predicted, by="ID", all.x=TRUE), df_predicted_sd, by="ID", all.x=TRUE)
    }
    fit = lm(Predicted ~ Observed, data=df); new_x=seq(from=min(df$Observed,na.rm=TRUE), to=max(df$Observed,na.rm=TRUE), length=1000)
    p = plotly::plot_ly(data=df, x=~Observed, y=~Predicted, type="scatter", mode='markers', text=~ID, error_y=~list(array=Predicted_sd, thickness=0.5, color="gray"))
    p = p %>% add_trace(x=new_x, y=predict(fit, newdata=data.frame(Observed=new_x)),
      mode='lines', 
      hoverinfo='text',
      text=~paste(
        "Slope:", round(coef(fit)[2], 4), 
        "<br>Intercept:", round(coef(fit)[1], 4),
        "<br>R-squared:", round(100*summary(fit)$r.squared), "%",
        "<br>Correaltion:", round(100*cor(df$Observed, df$Predicted)), "%",
        "<br>n = ", nrow(df),
        "<br>folds = ", length(unique(list_output$METRICS_ACROSS_POP_BULK$k))),
      error_y=NULL
    )
    p = p %>% layout(
      title=paste0(
        "Trait: ", list_output$TRAIT_NAME, 
        " | Model: ", input$across_bulk_scat_model
      ),
      showlegend=FALSE
    )
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
	return(p)
}

fn_across_bulk_histogram_server = function(input, list_list_output) {
    list_output = list_list_output[[input$across_bulk_scat_trait_x_pop]]
	if (is.null(list_output)) {
		df = data.frame(ID="No across populations data available", x=0, y=0)
		p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
		return(p)
	}
    idx = which(list_output$YPRED_ACROSS_POP_BULK$model == input$across_bulk_scat_model)
    df = data.frame(ID=list_output$YPRED_ACROSS_POP_BULK$id[idx], Observed=list_output$YPRED_ACROSS_POP_BULK$y_true[idx], Predicted=list_output$YPRED_ACROSS_POP_BULK$y_pred[idx])
    p = plot_ly(alpha=0.5, type="histogram", nbinsx=input$across_bulk_bins)
    p = p %>% add_histogram(data=df, x=~Observed, name="Observed")
    p = p %>% add_histogram(data=df, x=~Predicted, name="Predicted")
    p = p %>% layout(
      title="Phenotype distribution",
      barmode="overlay", 
      xaxis=list(title=list_output$TRAIT_NAME)
    )
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
}

fn_across_bulk_table_predictions_server = function(input, list_list_output) {
	list_output = list_list_output[[input$across_bulk_scat_trait_x_pop]]
	if (is.na(list_output[1])[1]) {
		return(NULL)
	}
    vec_idx = which(list_output$YPRED_ACROSS_POP_BULK$model == input$across_bulk_scat_model)
    df = data.frame(
		pop_training=list_output$YPRED_ACROSS_POP_BULK$pop_training[vec_idx],
		pop_validation=list_output$YPRED_ACROSS_POP_BULK$pop_validation[vec_idx],
		model=list_output$YPRED_ACROSS_POP_BULK$model[vec_idx],
		id=list_output$YPRED_ACROSS_POP_BULK$id[vec_idx],
		observed=list_output$YPRED_ACROSS_POP_BULK$y_true[vec_idx],
		predicted=list_output$YPRED_ACROSS_POP_BULK$y_pred[vec_idx])
	return(df)
}

####################################################################################################
### Front-end
####################################################################################################
ui <- page_fillable(
	titlePanel("plot(gp)"),
	navset_card_underline(
		nav_panel(h1(strong("WITHIN POPULATION"), style="font-size:15px; text-align:left"),
			fn_within_violin_ui(),
			fn_within_barplot_ui(),
			fn_within_scatter_hist_ui()
		),
		nav_panel(h1(strong("ACROSS POPULATIONS (PAIRWISE)"), style="font-size:15px; text-align:left"),
			fn_across_pair_barplot_ui(),
			fn_across_pair_scatter_hist_ui()
		),
		nav_panel(h1(strong("ACROSS POPULATIONS (LEAVE-ONE-OUT)"), style="font-size:15px; text-align:left"),
			fn_across_lopo_barplot_ui(),
			fn_across_lopo_scatter_hist_ui()
		),
		nav_panel(h1(strong("ACROSS POPULATIONS (BULK)"), style="font-size:15px; text-align:left"),
			fn_across_bulk_violin_ui(),
			fn_across_bulk_scatter_hist_ui()
		)
	)
	# fn_within_violin_ui(),
	# fn_within_scatter_hist_ui(),
	# fn_across_pair_barplot_ui(),
	# fn_across_pair_scatter_hist_ui(),
	# fn_across_lopo_barplot_ui(),
	# fn_across_lopo_scatter_hist_ui(),
	# fn_across_bulk_violin_ui(),
	# fn_across_bulk_scatter_hist_ui()
)
####################################################################################################
### Back-end
####################################################################################################
server = function(input, output, session) {
	##########################################
	### EXTRACT DATA AND UPDATE SELECTIONS ###
	##########################################
	data = reactive({
		# dirname_root = "/group/pasture/Jeff/gp/inst/exec_Rscript/output"
		# dirname_root = "/group/pasture/Jeff/lucerne/workdir/gs/output_ground_truth_biomass_traits/output"
		# dirname_root = "/group/pasture/Jeff/lucerne/workdir/gs/output_remote_sensing_biomass_traits/output"
		# dirname_root = "/group/pc3/Jeff/gs/per_trial_gp/output/"
		dirname_root = "/"
		shinyFiles::shinyDirChoose(input, id='dir', roots=c(root=dirname_root), filetypes=c('rds', 'Rds', 'RDS'), session=session)
		dir = as.character(parseDirPath(c(root=dirname_root), input$dir))
		list_list_output = fn_io_server(dir=dir)
		# list_list_output = fn_io_server(dir=NULL)
		#######################################
		### Update within population selections
		vec_traits = c()
		vec_populations = c()
		vec_models = c()
		vec_metrics = c()
		vec_entry_names = c()
		vec_pop_names = c()
		for (x in list_list_output) {
			# x = list_list_output[[1]]
			if (is.na(head(x$METRICS_WITHIN_POP, n=1)[1])) {
				next
			}
			vec_traits = c(vec_traits, as.character(x$TRAIT_NAME))
			vec_populations = c(vec_populations, unique(c(as.character(x$POPULATION), unique(x$METRICS_WITHIN_POP$pop_training))))
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
		updatePickerInput(session, "within_vec_traits", choices=vec_traits, selected=vec_traits[1])
		updatePickerInput(session, "within_vec_populations", choices=vec_populations, selected=vec_populations[1])
		updatePickerInput(session, "within_vec_models", choices=vec_models, selected=vec_models[1:2])
		updatePickerInput(session, "within_metric", choices=vec_metrics, selected=vec_metrics[1])
		updatePickerInput(session, "within_barplot_vec_traits", choices=vec_traits, selected=head(vec_traits))
		updatePickerInput(session, "within_barplot_vec_populations", choices=vec_populations, selected=head(vec_populations))
		updatePickerInput(session, "within_barplot_model", choices=vec_models, selected=vec_models[1])
		updatePickerInput(session, "within_barplot_metric", choices=vec_metrics, selected=vec_metrics[1])
		updatePickerInput(session, "within_scat_trait_x_pop", choices=names(list_list_output), selected=names(list_list_output)[1])
		updatePickerInput(session, "within_scat_model", choices=vec_models, selected=vec_models[1])
		#############################################################
		### Update across population selections (pairwise-population)
		vec_list_list_output_names = c()
		vec_traits = c()
		vec_pop_training = c()
		vec_pop_validation = c()
		vec_models = c()
		vec_metrics = c()
		for (i in 1:length(list_list_output)) {
			# i = 2
			x = list_list_output[[i]]
			if (!is.na(x$METRICS_ACROSS_POP_PAIRWISE[1])[1]) {
				vec_list_list_output_names = c(vec_list_list_output_names, names(list_list_output)[i])
				vec_traits = c(vec_traits, as.character(x$TRAIT_NAME))
				vec_pop_training = c(vec_pop_training, unique(as.character(x$METRICS_ACROSS_POP_PAIRWISE$pop_training)))
				vec_pop_validation = c(vec_pop_validation, unique(as.character(x$METRICS_ACROSS_POP_PAIRWISE$pop_validation)))
				vec_models = c(vec_models, unique(as.character(x$METRICS_ACROSS_POP_PAIRWISE$model)))
				for (j in 1:ncol(x$METRICS_ACROSS_POP_PAIRWISE)) {
					if (is.numeric(x$METRICS_ACROSS_POP_PAIRWISE[,j])) {
						vec_metrics = c(vec_metrics, colnames(x$METRICS_ACROSS_POP_PAIRWISE)[j])
					}
				}
			}
		}
		vec_list_list_output_names = sort(unique(vec_list_list_output_names))
		vec_traits = sort(unique(vec_traits))
		vec_pop_training = sort(unique(vec_pop_training))
		vec_pop_validation = sort(unique(vec_pop_validation))
		vec_models = sort(unique(vec_models))
		vec_metrics = sort(unique(vec_metrics))
		if (length(vec_traits) > 0) {
			updatePickerInput(session, "across_pair_trait", choices=vec_traits, selected=vec_traits[1])
			updatePickerInput(session, "across_pair_pop_training", choices=vec_pop_training, selected=vec_pop_training[1])
			updatePickerInput(session, "across_pair_vec_pop_validation", choices=vec_pop_validation, selected=vec_pop_validation)
			updatePickerInput(session, "across_pair_vec_models", choices=vec_models, selected=vec_models)
			updatePickerInput(session, "across_pair_metric", choices=vec_metrics, selected=vec_metrics[1])
			updatePickerInput(session, "across_pair_scat_trait", choices=vec_list_list_output_names, selected=vec_list_list_output_names[1])
			updatePickerInput(session, "across_pair_scat_pop_training", choices=vec_pop_training, selected=vec_pop_training[1])
			updatePickerInput(session, "across_pair_scat_pop_validation", choices=vec_pop_validation, selected=vec_pop_validation[2])
			updatePickerInput(session, "across_pair_scat_model", choices=vec_models, selected=vec_models)
		}
		##################################################################
		### Update across population selections (leave-one-population-out)
		vec_list_list_output_names = c()
		vec_traits = c()
		vec_pop_training = c()
		vec_pop_validation = c()
		vec_models = c()
		vec_metrics = c()
		for (i in 1:length(list_list_output)) {
			# i = 2
			x = list_list_output[[i]]
			if (!is.na(x$METRICS_ACROSS_POP_LOPO[1])[1]) {
				vec_list_list_output_names = c(vec_list_list_output_names, names(list_list_output)[i])
				vec_traits = c(vec_traits, as.character(x$TRAIT_NAME))
				vec_pop_training = c(vec_pop_training, unique(as.character(x$METRICS_ACROSS_POP_LOPO$pop_training)))
				vec_pop_validation = c(vec_pop_validation, unique(as.character(x$METRICS_ACROSS_POP_LOPO$pop_validation)))
				vec_models = c(vec_models, unique(as.character(x$METRICS_ACROSS_POP_LOPO$model)))
				for (j in 1:ncol(x$METRICS_ACROSS_POP_LOPO)) {
					if (is.numeric(x$METRICS_ACROSS_POP_LOPO[,j])) {
						vec_metrics = c(vec_metrics, colnames(x$METRICS_ACROSS_POP_LOPO)[j])
					}
				}
			}
		}
		vec_list_list_output_names = sort(unique(vec_list_list_output_names))
		vec_traits = sort(unique(vec_traits))
		vec_pop_training = sort(unique(vec_pop_training))
		vec_pop_validation = sort(unique(vec_pop_validation))
		vec_models = sort(unique(vec_models))
		vec_metrics = sort(unique(vec_metrics))
		if (length(vec_traits) > 0) {
			updatePickerInput(session, "across_lopo_trait", choices=vec_traits, selected=vec_traits[1])
			updatePickerInput(session, "across_lopo_pop_validation", choices=vec_pop_validation, selected=vec_pop_validation[1])
			updatePickerInput(session, "across_lopo_vec_models", choices=vec_models, selected=vec_models)
			updatePickerInput(session, "across_lopo_metric", choices=vec_metrics, selected=vec_metrics[1])
			updatePickerInput(session, "across_lopo_scat_trait", choices=vec_list_list_output_names, selected=vec_list_list_output_names[1])
			updatePickerInput(session, "across_lopo_scat_pop_validation", choices=vec_pop_validation, selected=vec_pop_validation[1])
			updatePickerInput(session, "across_lopo_scat_model", choices=vec_models, selected=vec_models)
		}
		##################################################################
		### Update across population selections (bulk)
		vec_traits = c()
		vec_populations = c()
		vec_trait_x_pop = c()
		vec_models = c()
		vec_metrics = c()
		vec_entry_names = c()
		vec_pop_names = c()
		for (i in 1:length(list_list_output)) {
			list_output = list_list_output[[i]]
			if (is.na(list_output$METRICS_ACROSS_POP_BULK[1])[1]) {next}
			vec_traits = c(vec_traits, as.character(list_output$TRAIT_NAME))
			vec_populations = c(vec_populations, as.character(list_output$POPULATION))
			vec_trait_x_pop = c(vec_trait_x_pop, names(list_list_output)[i])
			vec_models = c(vec_models, unique(as.character(list_output$METRICS_ACROSS_POP_BULK$model)))
			df_tmp = list_output$METRICS_ACROSS_POP_BULK
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
		updatePickerInput(session, "across_bulk_vec_traits", choices=vec_traits, selected=vec_traits[1])
		updatePickerInput(session, "across_bulk_vec_models", choices=vec_models, selected=vec_models)
		updatePickerInput(session, "across_bulk_metric", choices=vec_metrics, selected=vec_metrics[1])
		updatePickerInput(session, "across_bulk_scat_trait_x_pop", choices=vec_trait_x_pop, selected=vec_trait_x_pop[1])
		updatePickerInput(session, "across_bulk_scat_model", choices=vec_models, selected=vec_models[1])
		return(list_list_output)
	})
	############################
	### WITHIN POPULATION CV ###
	############################
	### Within population: violin or box plot
	output$within_plot_violin_or_box = renderPlotly({
		fn_within_violin_server(input=input, list_list_output=data())
	})
	### Within population: download selected raw and aggregated data based on user selections
	output$within_download_metrics = downloadHandler(
		filename = function() {
			paste0("gp_withinPop_table_metrics-", format(Sys.time(),'%Y%m%d%a%H%M%S'), ".zip")
		},
		content = function(file) {
			list_df_stats_raw = fn_within_table_metrics_server(input, list_list_output=data())
			dir_tmp = tempdir()
			vec_fnames_tsv = unlist(lapply(list_df_stats_raw, FUN=function(df){
				if (ncol(df) > ncol(list_df_stats_raw$df_stats)) {
					fname_tsv = "gp_withinPop_raw_selected_data.tsv"
					write.table(df, file=file.path(dir_tmp, fname_tsv), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
				} else {
					fname_tsv = "gp_withinPop_aggregated_stats.tsv"
					write.table(df, file=file.path(dir_tmp, fname_tsv), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
				}
				return(fname_tsv)
			}))
			zip::zip(zipfile=file, files=vec_fnames_tsv, root=dir_tmp)
		}
	)
	### Within population: summary and raw data tables
	output$within_data_tables= renderTable({
		list_df_stats_raw = fn_within_table_metrics_server(input, list_list_output=data())
		head(list_df_stats_raw$df_stats, n=4)
	})
	### Within population: bar plot
	output$within_plot_bar = renderPlotly({
		fn_within_barplot_server(input=input, list_list_output=data())
	})
	### Within population: observed vs prediction scatter plot
	output$within_plot_scatter = renderPlotly({
		fn_within_scatterplot_server(input=input, list_list_output=data())
	})
	### Within population: histogram of observed and predicted trait values
	output$within_plot_histogram = renderPlotly({
		fn_within_histogram_server(input, list_list_output=data())
	})
	### Within population: download selected observed and predicted phenotypes
	output$within_download_predictions = downloadHandler(
		filename = function() {
			paste0("gp_withinPop_table_predictions-", format(Sys.time(),'%Y%m%d%a%H%M%S'), ".tsv")
		},
		content = function(file) {
			df = fn_within_table_predictions_server(input, list_list_output=data())
			write.table(df, file=file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
		}
	)
	#####################################
	### ACROSS POPULATIONS (PAIRWISE) ###
	#####################################
	### Across populations (pairwise): barplot
	output$across_pair_plot_bar = renderPlotly({
		fn_across_pair_barplot_server(input, list_list_output=data())
	})
	### Across populations (pariwise): download selected raw and subset based on user selections
	output$across_pair_download_metrics = downloadHandler(
		filename = function() {
			paste0("gp_acrossPopPair_table_metrics-", format(Sys.time(),'%Y%m%d%a%H%M%S'), ".zip")
		},
		content = function(file) {
			list_df_sub_full = fn_across_pair_table_metrics_server(input, list_list_output=data())
			dir_tmp = tempdir()
			vec_fnames_tsv = unlist(lapply(list_df_sub_full, FUN=function(df){
				if (ncol(df) == 5) {
					fname_tsv = "gp_acrossPopPair_subset_data.tsv"
					write.table(df, file=file.path(dir_tmp, fname_tsv), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
				} else {
					fname_tsv = "gp_acrossPopPair_full_data.tsv"
					write.table(df, file=file.path(dir_tmp, fname_tsv), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
				}
				return(fname_tsv)
			}))
			zip::zip(zipfile=file, files=vec_fnames_tsv, root=dir_tmp)
		}
	)
	### Across populations (pairwise): data tables
	output$across_pair_data_tables = renderTable({
		list_df_sub_full = fn_across_pair_table_metrics_server(input, list_list_output=data())
		head(list_df_sub_full$df_sub, n=3)
	})
	### Across populations (pairwise): observed vs predicted scatter plot
	output$across_pair_plot_scatter = renderPlotly({
		fn_across_pair_scatter_server(input, list_list_output=data())
	})
	### Across populations (pairwise): histogram of observed and predicted phenotypes
	output$across_pair_plot_histogram = renderPlotly({
		fn_across_pair_histogram_server(input, list_list_output=data())
	})
	### Across populations (pairwise): download selected observed and predicted phenotypes
	output$across_pair_download_predictions = downloadHandler(
		filename = function() {
			paste0("gp_acrossPopPair_table_predictions-", format(Sys.time(),'%Y%m%d%a%H%M%S'), ".tsv")
		},
		content = function(file) {
			df = fn_across_pair_table_predictions_server(input, list_list_output=data())
			write.table(df, file=file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
		}
	)
	#####################################################
	### ACROSS POPULATIONS (LEAVE-ONE-POPULATION-OUT) ###
	#####################################################
	### Across populations (leave-one-population-out): barplot
	output$across_lopo_plot_bar = renderPlotly({
		fn_across_lopo_barplot_server(input, list_list_output=data())
	})
	### Across populations (leave-one-population-out): download selected raw and subset based on user selections
	output$across_lopo_download_metrics = downloadHandler(
		filename = function() {
			paste0("gp_acrossPopLOPO_table_metrics-", format(Sys.time(),'%Y%m%d%a%H%M%S'), ".zip")
		},
		content = function(file) {
			list_df_sub_full = fn_across_lopo_table_metrics_server(input, list_list_output=data())
			dir_tmp = tempdir()
			vec_fnames_tsv = unlist(lapply(list_df_sub_full, FUN=function(df){
				if (ncol(df) == 5) {
					fname_tsv = "gp_acrossPopLOPO_subset_data.tsv"
					write.table(df, file=file.path(dir_tmp, fname_tsv), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
				} else {
					fname_tsv = "gp_acrossPopLOPO_full_data.tsv"
					write.table(df, file=file.path(dir_tmp, fname_tsv), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
				}
				return(fname_tsv)
			}))
			zip::zip(zipfile=file, files=vec_fnames_tsv, root=dir_tmp)
		}
	)
	### Across populations (leave-one-population-out): data tables
	output$across_lopo_data_tables = renderTable({
		list_df_sub_full = fn_across_lopo_table_metrics_server(input, list_list_output=data())
		head(list_df_sub_full$df_sub, n=3)
	})
	### Across populations (leave-one-population-out): observed vs predicted scatter plot
	output$across_lopo_plot_scatter = renderPlotly({
		fn_across_lopo_scatter_server(input, list_list_output=data())
	})
	### Across populations (leave-one-population-out): histogram of observed and predicted phenotypes
	output$across_lopo_plot_histogram = renderPlotly({
		fn_across_lopo_histogram_server(input, list_list_output=data())
	})
	### Across populations (leave-one-population-out): download selected observed and predicted phenotypes
	output$across_lopo_download_predictions = downloadHandler(
		filename = function() {
			paste0("gp_acrossPopLOPO_table_predictions-", format(Sys.time(),'%Y%m%d%a%H%M%S'), ".tsv")
		},
		content = function(file) {
			df = fn_across_lopo_table_predictions_server(input, list_list_output=data())
			write.table(df, file=file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
		}
	)
	#################################
	### ACROSS POPULATIONS (BULK) ###
	#################################
	### Across population (bulk): violin or box plot
	output$across_bulk_plot_violin_or_box = renderPlotly({
		fn_across_bulk_violin_server(input=input, list_list_output=data())
	})
	### Across population (bulk): download selected raw and aggregated data based on user selections
	output$across_bulk_download_metrics = downloadHandler(
		filename = function() {
			paste0("gp_acrossPopBulk_table_metrics-", format(Sys.time(),'%Y%m%d%a%H%M%S'), ".zip")
		},
		content = function(file) {
			list_df_stats_raw = fn_across_bulk_table_metrics_server(input, list_list_output=data())
			dir_tmp = tempdir()
			vec_fnames_tsv = unlist(lapply(list_df_stats_raw, FUN=function(df){
				if (ncol(df) > ncol(list_df_stats_raw$df_stats)) {
					fname_tsv = "gp_acrossPopBulk_raw_selected_data.tsv"
					write.table(df, file=file.path(dir_tmp, fname_tsv), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
				} else {
					fname_tsv = "gp_acrossPopBulk_aggregated_stats.tsv"
					write.table(df, file=file.path(dir_tmp, fname_tsv), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
				}
				return(fname_tsv)
			}))
			zip::zip(zipfile=file, files=vec_fnames_tsv, root=dir_tmp)
		}
	)
	### Across population (bulk): summary and raw data tables
	output$across_bulk_data_tables= renderTable({
		list_df_stats_raw = fn_across_bulk_table_metrics_server(input, list_list_output=data())
		head(list_df_stats_raw$df_stats, n=4)
	})
	### Across population (bulk): observed vs prediction scatter plot
	output$across_bulk_plot_scatter = renderPlotly({
		fn_across_bulk_scatterplot_server(input=input, list_list_output=data())
	})
	### Across population (bulk): histogram of observed and predicted trait values
	output$across_bulk_plot_histogram = renderPlotly({
		fn_across_bulk_histogram_server(input, list_list_output=data())
	})
	### Across population (bulk): download selected observed and predicted phenotypes
	output$across_bulk_download_predictions = downloadHandler(
		filename = function() {
			paste0("gp_acrossPopBulk_table_predictions-", format(Sys.time(),'%Y%m%d%a%H%M%S'), ".tsv")
		},
		content = function(file) {
			df = fn_across_bulk_table_predictions_server(input, list_list_output=data())
			write.table(df, file=file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
		}
	)
}
####################################################################################################
### Serve the app
####################################################################################################
options(shiny.maxRequestSize = 50 * 1024^2)
shinyApp(ui = ui, server = server)

