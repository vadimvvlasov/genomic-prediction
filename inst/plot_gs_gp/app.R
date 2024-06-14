library(shiny)
library(shinyWidgets)
library(plotly)
library(bslib)
library(shinyFiles)

####################################################################################################
### Define the light and dark themes
####################################################################################################
light = bslib::bs_theme()
dark = bslib::bs_theme(
  bg = "#202123", fg = "#B8BCC2",
  primary = "#EA80FC", secondary = "#48DAC6",
  base_font = c("Grandstander", "sans-serif"),
  code_font = c("Courier", "monospace"),
  heading_font = "'Helvetica Neue', Helvetica, sans-serif",
  "input-border-color" = "#EA80FC"
)
####################################################################################################
### Front-end
####################################################################################################
ui <- page_fillable(
  titlePanel("PLOT(genomic_selection)"),
  theme = light,
  shinyWidgets::materialSwitch(inputId="dark_mode", label="Dark mode", value=FALSE, status="primary", right=TRUE),
  br(),
  card(
    card_header(h1(strong("Input files and main plot"), style="font-size:21px; text-align:left")),
    min_height="840px",
    layout_sidebar(
      sidebar=sidebar(
         width=500,
        fileInput(
          inputId="input", 
          label="Choose the Rds output files from the genomic_prediction workflow",
          multiple=TRUE,
          accept=".Rds"
        ),
        shinyWidgets::pickerInput(inputId="trait", label="Filter by trait:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="pop", label="Filter by population:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="model", label="Filter by model:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="metric", label="Use the genomic prediction accuracy metric:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="group_by", label="Group by:", choices=c("trait", "pop_training", "model"), selected=c("trait", "pop", "model"), multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="sort_by", label="Sort by:", choices=c("increasing_mean", "decreasing_mean", "increasing_median", "decreasing_median", "alphabetical"), selected="increasing_mean", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::materialSwitch(inputId="within_box_with_labels", label="Mean-labelled boxplot", value=FALSE, status="primary", right=FALSE)
      ),
      mainPanel(
        width=750,
        plotlyOutput(outputId="plot_violin"),
        verbatimTextOutput(outputId="data_sizes1")
      )
    )
  ),
  card(
    card_header(h1(strong("Per trait-population-model combination"), style="font-size:21px; text-align:left")),
    min_height="1200px",
    layout_sidebar(
      sidebar=sidebar(
        width=500,
        shinyWidgets::pickerInput(inputId="scat_trait_x_pop", label="Trait-Pop combination:", choices="", options=list(`live-search`=TRUE)),
        shinyWidgets::pickerInput(inputId="scat_model", label="Model:", choices="", options=list(`live-search`=TRUE)),
        shinyWidgets::materialSwitch(inputId="scat_points", label="Plot points", value=FALSE, status="primary", right=TRUE),
        sliderInput(inputId="bins", label="Number of bins:", min=1, max=100, value=10)
      ),
      mainPanel(
        width=750,
        plotlyOutput(outputId="plot_scatter"),
        br(),
        br(),
        br(),
        br(),
        plotlyOutput(outputId="plot_histogram")
      )
    )
  ),
  card(
    card_header(h1(strong("LOPO (leave-one-population-out) per trait-model combination"), style="font-size:21px; text-align:left")),
    min_height="1200px",
    layout_sidebar(
      sidebar=sidebar(
        width=500,
        shinyWidgets::pickerInput(inputId="scat_lopo_trait", label="Trait (Note: Ignore the population label as the first population contains the LOPO output):", choices="", options=list(`live-search`=TRUE)),
        shinyWidgets::pickerInput(inputId="scat_lopo_model", label="Model:", choices="", options=list(`live-search`=TRUE)),
        sliderInput(inputId="lopo_bins", label="Number of bins:", min=1, max=100, value=10)
      ),
      mainPanel(
        width=750,
        plotlyOutput(outputId="plot_scatter_lopo"),
        br(),
        br(),
        br(),
        br(),
        plotlyOutput(outputId="plot_histogram_lopo")
      )
    )
  ),
  card(
    card_header(h1(strong("LOPO barplot"), style="font-size:21px; text-align:left")),
    min_height="750px",
    layout_sidebar(
      sidebar=sidebar(
        width=500,
        shinyWidgets::pickerInput(inputId="lopo_trait", label="Filter by trait:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="lopo_metric", label="Use the genomic prediction accuracy metric:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::materialSwitch(inputId="lopo_group_model", label="Group by model", value=FALSE, status="primary", right=TRUE),
        shinyWidgets::pickerInput(inputId="lopo_pop", label="Filter by pop:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="lopo_model", label="Filter by model:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE))
      ),
      mainPanel(
        width=750,
        plotlyOutput(outputId="plot_bar_lopo"),
        verbatimTextOutput(outputId="data_sizes2")
      )
    )
  ),
  card(
    card_header(h1(strong("PAIRWISE (pairwise-population cross-validation) per trait-model combination"), style="font-size:21px; text-align:left")),
    min_height="1200px",
    layout_sidebar(
      sidebar=sidebar(
        width=500,
        shinyWidgets::pickerInput(inputId="scat_pairwise_trait", label="Trait (Note: Ignore the population label as the first population contains the PAIRWISE output):", choices="", options=list(`live-search`=TRUE)),
        shinyWidgets::pickerInput(inputId="scat_pairwise_pop_training", label="Training population:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="scat_pairwise_pop_validation", label="Validation population:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="scat_pairwise_model", label="Model:", choices="", options=list(`live-search`=TRUE)),
        sliderInput(inputId="pairwise_bins", label="Number of bins:", min=1, max=100, value=10)
      ),
      mainPanel(
        width=750,
        plotlyOutput(outputId="plot_scatter_pairwise"),
        br(),
        br(),
        br(),
        br(),
        plotlyOutput(outputId="plot_histogram_pairwise")
      )
    )
  ),
  card(
    card_header(h1(strong("PAIRWISE barplot"), style="font-size:21px; text-align:left")),
    min_height="750px",
    layout_sidebar(
      sidebar=sidebar(
        width=500,
        shinyWidgets::pickerInput(inputId="pairwise_trait", label="Filter by trait:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="pairwise_metric", label="Use the genomic prediction accuracy metric:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::materialSwitch(inputId="pairwise_group_model", label="Group by model", value=FALSE, status="primary", right=TRUE),
        shinyWidgets::pickerInput(inputId="pairwise_pop_training", label="Filter by training population:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="pairwise_pop_validation", label="Filter by validation population:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="pairwise_model", label="Filter by model:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE))
      ),
      mainPanel(
        width=750,
        plotlyOutput(outputId="plot_bar_pairwise"),
        verbatimTextOutput(outputId="data_sizes3")

      )
    )
  ),
  card(
    card_header(h1(strong("Genomic predictions per se"), style="font-size:21px; text-align:left")),
    min_height="500px",
    layout_sidebar(
      sidebar=sidebar(
        width=500,
        shinyWidgets::pickerInput(inputId="gp_trait", label="Filter by trait:", choices="", multiple=FALSE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        shinyWidgets::pickerInput(inputId="gp_pop", label="Filter by pop:", choices="", multiple=TRUE, options=list(`live-search`=TRUE, `actions-box`=TRUE)),
        sliderInput(inputId="gp_bins", label="Number of bins:", min=1, max=100, value=10)
      ),
      mainPanel(
        width=750,
        plotlyOutput(outputId="plot_genomic_predictions_per_se")
      )
    )
  ),
  card(
    card_header(h1(strong("Debugging"), style="font-size:21px; text-align:left")),
    min_height="500px",
    layout_sidebar(
      mainPanel(
        verbatimTextOutput(outputId="debug")
      )
    )
  )
)
####################################################################################################
### Back-end
####################################################################################################
server <- function(input, output, session) {
  observe(session$setCurrentTheme(
    if (isTRUE(input$dark_mode)) dark else light
  ))
  #######################################################################
  ### Load input Rds file/s
  #######################################################################
  # roots <- getVolumes()
  # shinyFileChoose(input, "file", roots = roots)
  # file <- reactive(parseFilePaths(roots, input$file))
  data = reactive({
    if (is.null(input$input$datapath)) {
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
    } else {
      ### Load the user-input data
      vec_Rds = input$input$datapath
      # vec_Rds = list.files(".", pattern="*.Rds")
      list_list_output = list()
      for (i in 1:length(vec_Rds)) {
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
        validate(
          need(class(list_output) == "list", paste0("Error: the ", nth, " Rds input is not the correct genomic_selection pipeline output (not a list)!")),
          need(length(list_output) == 12, paste0("Error: the ", nth, " Rds input is not the correct genomic_selection pipeline output (does not have 7 elements)!"))
        )
        validate(
          need(!is.null(list_output$TRAIT_NAME), paste0("Missing filed: TRAIT_NAME")),
          need(!is.null(list_output$POPULATION), paste0("Missing filed: POPULATION")),
          need(!is.null(list_output$METRICS_WITHIN_POP), paste0("Missing filed: METRICS_WITHIN_POP")),
          need(!is.null(list_output$YPRED_WITHIN_POP), paste0("Missing filed: YPRED_WITHIN_POP")),
          need(!is.null(list_output$METRICS_ACROSS_POP_BULK), paste0("Missing filed: METRICS_ACROSS_POP_BULK")),
          need(!is.null(list_output$YPRED_ACROSS_POP_BULK), paste0("Missing filed: YPRED_ACROSS_POP_BULK")),
          need(!is.null(list_output$METRICS_ACROSS_POP_PAIRWISE), paste0("Missing filed: METRICS_ACROSS_POP_PAIRWISE")),
          need(!is.null(list_output$YPRED_ACROSS_POP_PAIRWISE), paste0("Missing filed: YPRED_ACROSS_POP_PAIRWISE")),
          need(!is.null(list_output$METRICS_ACROSS_POP_LOPO), paste0("Missing filed: METRICS_ACROSS_POP_LOPO")),
          need(!is.null(list_output$YPRED_ACROSS_POP_LOPO), paste0("Missing filed: YPRED_ACROSS_POP_LOPO")),
          need(!is.null(list_output$GENOMIC_PREDICTIONS), paste0("Missing filed: GENOMIC_PREDICTIONS")),
          need(!is.null(list_output$ADDITIVE_GENETIC_EFFECTS), paste0("Missing filed: ADDITIVE_GENETIC_EFFECTS"))
        )
        trait = list_output$TRAIT_NAME
        pop = list_output$POPULATION
        eval(parse(text=paste0("list_list_output$`", trait, "_", pop, "` = list_output")))
      }
    }
    ### Within population data
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
    updatePickerInput(session, "trait", choices=vec_traits, selected=vec_traits[1])
    updatePickerInput(session, "pop", choices=vec_populations, selected=vec_populations[1])
    updatePickerInput(session, "model", choices=vec_models, selected=vec_models)
    updatePickerInput(session, "metric", choices=vec_metrics, selected=vec_metrics[1])
    updatePickerInput(session, "scat_trait_x_pop", choices=names(list_list_output), selected=names(list_list_output)[1])
    updatePickerInput(session, "scat_model", choices=vec_models, selected=vec_models[1])
    ### Across population data (leave-one-population-out cross-validation)
    vec_list_list_output_names = c()
    vec_traits_lopo = c()
    vec_populations_lopo = c()
    vec_models_lopo = c()
    vec_metrics_lopo = c()
    for (i in 1:length(list_list_output)) {
      # i = 2
      x = list_list_output[[i]]
      if (!is.na(x$METRICS_ACROSS_POP_LOPO[1])[1]) {
        vec_list_list_output_names = c(vec_list_list_output_names, names(list_list_output)[i])
        vec_traits_lopo = c(vec_traits_lopo, as.character(x$TRAIT_NAME))
        vec_populations_lopo = c(vec_populations_lopo, unique(as.character(x$METRICS_ACROSS_POP_LOPO$pop_validation)))
        vec_models_lopo = c(vec_models_lopo, unique(as.character(x$METRICS_ACROSS_POP_LOPO$model)))
        for (j in 1:ncol(x$METRICS_ACROSS_POP_LOPO)) {
          if (is.numeric(x$METRICS_ACROSS_POP_LOPO[,j])) {
            vec_metrics_lopo = c(vec_metrics_lopo, colnames(x$METRICS_ACROSS_POP_LOPO)[j])
          }
        }
      }
    }
    vec_list_list_output_names = sort(unique(vec_list_list_output_names))
    vec_traits_lopo = sort(unique(vec_traits_lopo))
    vec_populations_lopo = sort(unique(vec_populations_lopo))
    vec_models_lopo = sort(unique(vec_models_lopo))
    vec_metrics_lopo = sort(unique(vec_metrics_lopo))
    if (length(vec_traits_lopo) > 0) {
      updatePickerInput(session, "lopo_trait", choices=vec_traits_lopo, selected=vec_traits_lopo[1])
      updatePickerInput(session, "lopo_pop", choices=vec_populations_lopo, selected=vec_populations_lopo)
      updatePickerInput(session, "lopo_model", choices=vec_models_lopo, selected=vec_models_lopo)
      updatePickerInput(session, "lopo_metric", choices=vec_metrics_lopo, selected=vec_metrics_lopo[1])

      updatePickerInput(session, "scat_lopo_trait", choices=vec_list_list_output_names, selected=vec_list_list_output_names[1])
      updatePickerInput(session, "scat_lopo_model", choices=vec_models_lopo, selected=vec_models_lopo)

    }

















    ### Across population data (leave-one-population-out cross-validation)
    vec_list_list_output_names = c()
    vec_traits_pairwise = c()
    vec_populations_training = c()
    vec_populations_validation = c()
    vec_models_pairwise = c()
    vec_metrics_pairwise = c()
    for (i in 1:length(list_list_output)) {
      # i = 1
      x = list_list_output[[i]]
      if (!is.na(x$METRICS_ACROSS_POP_PAIRWISE[1])[1]) {
        vec_list_list_output_names = c(vec_list_list_output_names, names(list_list_output)[i])
        vec_traits_pairwise = c(vec_traits_pairwise, as.character(x$TRAIT_NAME))
        vec_populations_training = c(vec_populations_training, unique(as.character(x$METRICS_ACROSS_POP_PAIRWISE$pop_training)))
        vec_populations_validation = c(vec_populations_validation, unique(as.character(x$METRICS_ACROSS_POP_PAIRWISE$pop_validation)))
        vec_models_pairwise = c(vec_models_pairwise, unique(as.character(x$METRICS_ACROSS_POP_PAIRWISE$model)))
        for (j in 1:ncol(x$METRICS_ACROSS_POP_PAIRWISE)) {
          if (is.numeric(x$METRICS_ACROSS_POP_PAIRWISE[,j])) {
            vec_metrics_pairwise = c(vec_metrics_pairwise, colnames(x$METRICS_ACROSS_POP_PAIRWISE)[j])
          }
        }
      }
    }
    vec_list_list_output_names = sort(unique(vec_list_list_output_names))
    vec_traits_pairwise = sort(unique(vec_traits_pairwise))
    vec_populations_training = sort(unique(vec_populations_training))
    vec_populations_validation = sort(unique(vec_populations_validation))
    vec_models_pairwise = sort(unique(vec_models_pairwise))
    vec_metrics_pairwise = sort(unique(vec_metrics_pairwise))
    if (length(vec_traits_pairwise) > 0) {
      updatePickerInput(session, "pairwise_trait", choices=vec_traits_pairwise, selected=vec_traits_pairwise[1])
      updatePickerInput(session, "pairwise_model", choices=vec_models_pairwise, selected=vec_models_pairwise)
      updatePickerInput(session, "pairwise_metric", choices=vec_metrics_pairwise, selected=vec_metrics_pairwise[1])
      updatePickerInput(session, "pairwise_pop_training", choices=vec_populations_training, selected=vec_populations_training[1])
      updatePickerInput(session, "pairwise_pop_validation", choices=vec_populations_validation, selected=vec_populations_validation)

      updatePickerInput(session, "scat_pairwise_trait", choices=vec_list_list_output_names, selected=vec_list_list_output_names[1])
      updatePickerInput(session, "scat_pairwise_model", choices=vec_models_pairwise, selected=vec_models_pairwise)
      updatePickerInput(session, "scat_pairwise_pop_training", choices=vec_populations_training, selected=vec_populations_training[1])
      updatePickerInput(session, "scat_pairwise_pop_validation", choices=vec_populations_validation, selected=vec_populations_validation[2])

    }



















    ### Genomic predictions per se
    vec_traits_gp = c()
    vec_populations_gp = c()
    for (x in list_list_output) {
      # x = list_list_output[[4]]
      if (!is.na(x$GENOMIC_PREDICTION[1])[1]) {
        vec_traits_gp = c(vec_traits_gp, as.character(x$TRAIT_NAME))
        vec_populations_gp = c(vec_populations_gp, as.character(x$GENOMIC_PREDICTION$pop))
      }
    }
    vec_traits_gp = sort(unique(vec_traits_gp))
    vec_populations_gp = sort(unique(vec_populations_gp))
    if (length(vec_traits_gp) > 0) {
      updatePickerInput(session, "gp_trait", choices=vec_traits_gp, selected=vec_traits_gp[1])
      updatePickerInput(session, "gp_pop", choices=vec_populations_gp, selected=vec_populations_gp)
    }
    return(list_list_output)
  })
  #######################################################################
  ### WITHIN: violin plot responsive to various aggregation and sorting modes
  #######################################################################
  output$plot_violin = renderPlotly({
    ### Define the grouping, i.e. trait and/or pop and/or model
    list_list_output = data()
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
    idx = (METRICS_WITHIN_POP$trait %in% input$trait) & (METRICS_WITHIN_POP$pop_training %in% input$pop) & (METRICS_WITHIN_POP$model %in% input$model)
    grouping = c()
    x_label_one_level = c()
    x_label_one_level_names = c()
    x_label_multiple_levels = c()
    for (g in input$group_by) {
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
    validate(
      need(grouping > 0, "Error: Please select a grouping factor ('Group by:' selection) with more than two levels selected.")
    )
    df = data.frame(
      grouping=grouping,
      metric=eval(parse(text=paste0("METRICS_WITHIN_POP$", input$metric, "[idx]")))
    )
    if (input$sort_by == "increasing_mean") {
      df_agg = aggregate(metric ~ grouping, data=df, FUN=mean, na.rm=TRUE)
      df$grouping = factor(df$grouping, levels=as.character(df_agg$grouping[order(df_agg$metric, decreasing=FALSE)]))
    }
    if (input$sort_by == "decreasing_mean") {
      df_agg = aggregate(metric ~ grouping, data=df, FUN=mean, na.rm=TRUE)
      df$grouping = factor(df$grouping, levels=as.character(df_agg$grouping[order(df_agg$metric, decreasing=TRUE)]))
    }
    if (input$sort_by == "increasing_median") {
      df_agg = aggregate(metric ~ grouping, data=df, FUN=median, na.rm=TRUE)
      df$grouping = factor(df$grouping, levels=as.character(df_agg$grouping[order(df_agg$metric, decreasing=FALSE)]))
    }
    if (input$sort_by == "decreasing_median") {
      df_agg = aggregate(metric ~ grouping, data=df, FUN=median, na.rm=TRUE)
      df$grouping = factor(df$grouping, levels=as.character(df_agg$grouping[order(df_agg$metric, decreasing=TRUE)]))
    }
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
        data.frame(grouping=list_bp$names, pos_x=seq(0, n_groups+1, length=n_groups+1)[1:n_groups], pos_y=1.05*min(df$metric, na.rm=TRUE)),
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
        showarrow=FALSE
      )
    }
    p = p %>% layout(
      title=title,
      yaxis=list(title=input$metric),
      xaxis=list(title="")
    )
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
  })
  #######################################################################
  ### DATA SIZES 1
  #######################################################################
  output$data_sizes1= renderPrint({ 
    list_list_output = data()
    vec_traits = c()
    vec_populations = c()
    vec_pop_size = c()
    vec_fold_size = c()
    for (list_output in list_list_output) {
      # list_output = list_list_output[[1]]
      vec_traits = c(vec_traits, list_output$TRAIT_NAME)
      vec_populations = c(vec_populations, unique(list_output$METRICS_WITHIN_POP$pop_training))
      vec_pop_size = c(vec_pop_size, length(unique(list_output$YPRED_WITHIN_POP$id)))
      vec_fold_size = c(vec_fold_size, max(list_output$YPRED_WITHIN_POP$fold))
    }
    df = data.frame(trait=vec_traits, population=vec_populations, size=vec_pop_size, fold=vec_fold_size)
    print(df[(df$trait %in% input$trait) & (df$population %in% input$pop), ])
  })
  #######################################################################
  ### WITHIN: Scatter plot per trait-by-population combination and model
  #######################################################################
  output$plot_scatter = renderPlotly({
    list_list_output = data()
    list_output = list_list_output[[input$scat_trait_x_pop]]
    idx = which(list_output$YPRED_WITHIN_POP$model == input$scat_model)
    df_tmp = data.frame(ID=list_output$YPRED_WITHIN_POP$id[idx], Observed=list_output$YPRED_WITHIN_POP$y_true[idx], Predicted=list_output$YPRED_WITHIN_POP$y_pred[idx])
    if (input$scat_points) {
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
        " | Population: ", unique(list_output$SUMMARY$pop)[grep("overall", unique(list_output$SUMMARY$pop), invert=TRUE)], 
        " | Model: ", input$scat_model
      ),
      showlegend=FALSE
    )
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
  })
  #######################################################################
  ### WITHIN: Distribution of observed and predicted trait values per trait-by-population combination and model
  #######################################################################
  output$plot_histogram = renderPlotly({
    list_list_output = data()
    list_output = list_list_output[[input$scat_trait_x_pop]]
    idx = which(list_output$YPRED_WITHIN_POP$model == input$scat_model)
    df = data.frame(ID=list_output$YPRED_WITHIN_POP$id[idx], Observed=list_output$YPRED_WITHIN_POP$y_true[idx], Predicted=list_output$YPRED_WITHIN_POP$y_pred[idx])
    p = plot_ly(alpha=0.5, type="histogram", nbinsx=input$bins)
    p = p %>% add_histogram(data=df, x=~Observed, name="Observed")
    p = p %>% add_histogram(data=df, x=~Predicted, name="Predicted")
    p = p %>% layout(
      title="Phenotype distribution",
      barmode="overlay", 
      xaxis=list(title=list_output$TRAIT_NAME)
    )
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
  })
  #######################################################################
  ### LOPO: Scatter plot per trait-by-population combination and model
  #######################################################################
  output$plot_scatter_lopo = renderPlotly({
    list_list_output = data()
    list_output = list_list_output[[input$scat_lopo_trait]]
    if (!is.na(list_output$YPRED_ACROSS_POP_LOPO[1])[1]) {
      # input = list(scat_lopo_model="Bayes_A")
      idx = which(as.character(list_output$YPRED_ACROSS_POP_LOPO$model) == input$scat_lopo_model)
      df = data.frame(
        ID=list_output$YPRED_ACROSS_POP_LOPO$id[idx], 
        POP=list_output$YPRED_ACROSS_POP_LOPO$pop_validation[idx], 
        Observed=list_output$YPRED_ACROSS_POP_LOPO$y_true[idx], 
        Predicted=list_output$YPRED_ACROSS_POP_LOPO$y_pred[idx])
    } else {
      df = data.frame(ID="No leave-one-population-out data available", POP="", Observed=0, Predicted=0)
    }
    fit = lm(Predicted ~ Observed, data=df); new_x=seq(from=min(df$Observed,na.rm=TRUE), to=max(df$Observed,na.rm=TRUE), length=1000)
    p = plotly::plot_ly(data=df, x=~Observed, y=~Predicted, type="scatter", mode='markers', text=~ID, color=~POP, showlegend=TRUE)
    p = p %>% add_trace(x=new_x, y=predict(fit, newdata=data.frame(Observed=new_x)), 
      color=NULL,
      name="linear fit",
      mode='lines', 
      hoverinfo='text',
      text=~paste(
        "Slope:", round(coef(fit)[["Observed"]], 4), 
        "<br>Intercept:", round(coef(fit)[["(Intercept)"]], 4),
        "<br>R-squared:", round(100*summary(fit)$r.squared), "%",
        "<br>Correaltion:", round(100*cor(df$Observed, df$Predicted, use="na.or.complete")), "%",
        "<br>n = ", nrow(df),
        "<br>folds = ", length(unique(list_output$METRICS_ACROSS_POP_LOPO$fold))),
      error_y=NULL
    )
    p = p %>% layout(
      title=paste0(
        "Trait: ", list_output$TRAIT_NAME, 
        " | Model: ", input$scat_lopo_model
      )
    )
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
  })
  #######################################################################
  ### LOPO: Distribution of observed and predicted trait values per trait-by-population combination and model
  #######################################################################
  output$plot_histogram_lopo = renderPlotly({
    list_list_output = data()
    list_output = list_list_output[[input$scat_lopo_trait]]
    if (!is.na(list_output$YPRED_ACROSS_POP_LOPO[1])[1]) {
      idx = which(as.character(list_output$YPRED_ACROSS_POP_LOPO$model) == input$scat_lopo_model)
      df = data.frame(
        ID=list_output$YPRED_ACROSS_POP_LOPO$id[idx], 
        Observed=list_output$YPRED_ACROSS_POP_LOPO$y_true[idx], 
        Predicted=list_output$YPRED_ACROSS_POP_LOPO$y_pred[idx])
    } else {
      df = data.frame(ID="No leave-one-population-out data available", Observed=0, Predicted=0)
    }
    p = plot_ly(alpha=0.5, type="histogram", nbinsx=input$lopo_bins)
    p = p %>% add_histogram(data=df, x=~Observed, name="Observed")
    p = p %>% add_histogram(data=df, x=~Predicted, name="Predicted")
    p = p %>% layout(
      title="Phenotype distribution",
      barmode="overlay", 
      xaxis=list(title=list_output$TRAIT_NAME)
    )
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
  })
  #######################################################################
  ### LOPO: bar plot (no replication as entire populations are used hence no subsampling performed)
  #######################################################################
  output$plot_bar_lopo = renderPlotly({
    ### Define the grouping, i.e. trait and/or pop and/or model
    list_list_output = data()
    df_metrics_across_pop = NULL
    for (list_output in list_list_output) {
      if (list_output$TRAIT_NAME != input$lopo_trait) {
        next
      }
      if (is.null(list_output$METRICS_ACROSS_POP_LOPO)) {
        next
      }
      if (is.null(df_metrics_across_pop)) {
        df_metrics_across_pop = list_output$METRICS_ACROSS_POP_LOPO
      } else {
        df_metrics_across_pop = rbind(df_metrics_across_pop, list_output$METRICS_ACROSS_POP_LOPO)
      }
    }
    if (is.null(df_metrics_across_pop)) {
      df = data.frame(ID="No leave-one-population-out data available", x=0, y=0)
      p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
    } else {
      ### Prepare pop_validation x model matrix of GP metric data
      vec_lopo_pop = sort(input$lopo_pop)
      n = length(vec_lopo_pop)
      m = length(input$lopo_model)
      M = matrix(NA, nrow=n, ncol=m)
      rownames(M) = vec_lopo_pop
      colnames(M) = input$lopo_model
      for (i in 1:n) {
        for (j in 1:m) {
          idx = which((df_metrics_across_pop$pop_validation == vec_lopo_pop[i]) & (df_metrics_across_pop$model == input$lopo_model[j]))
          M[i, j] = eval(parse(text=paste0("df_metrics_across_pop$", input$lopo_metric, "[idx]")))
        }
      }
      vec_validation_pops = sort(unique(as.character(df_metrics_across_pop$pop_validation)))
      vec_models = sort(unique(as.character(df_metrics_across_pop$model)))
      if(!input$lopo_group_model) {
        df = data.frame(vec_lopo_pop, M)
        colnames(df) = c("pop_validation", colnames(M))
        ### Plot
        eval(parse(text=paste0("p = plotly::plot_ly(data=df, x=~pop_validation, y=~", input$lopo_model[1], ", type='bar', name='", input$lopo_model[1], "')")))
        vec_metric = eval(parse(text=paste0("df$", input$lopo_model[1])))
        if (length(input$lopo_model) > 1) {
          for (i in 2:length(input$lopo_model)) {
            eval(parse(text=paste0("p = p %>% plotly::add_trace(y=~", input$lopo_model[i], ", name='", input$lopo_model[i], "')")))
            vec_metric = c(vec_metric, eval(parse(text=paste0("df$", input$lopo_model[i]))))
          }
        }
      } else {
        df = data.frame(input$lopo_model, t(M))
        colnames(df) = c("model", rownames(M))
        ### Plot
        eval(parse(text=paste0("p = plotly::plot_ly(data=df, x=~model, y=~", vec_lopo_pop[1], ", type='bar', name='", vec_lopo_pop[1], "')")))
        vec_metric = eval(parse(text=paste0("df$", vec_lopo_pop[1])))
        if (length(vec_lopo_pop) > 1) {
          for (i in 2:length(vec_lopo_pop)) {
            eval(parse(text=paste0("p = p %>% plotly::add_trace(y=~", vec_lopo_pop[i], ", name='", vec_lopo_pop[i], "')")))
            vec_metric = c(vec_metric, eval(parse(text=paste0("df$", vec_lopo_pop[i]))))
          }
        }
        
      }
      p = p %>% plotly::layout(
        title=paste0("Trait: ", input$lopo_trait, " | Mean ", input$lopo_metric, " = ", round(mean(vec_metric, na.rm=TRUE), 4)), 
        barmode="group",
        yaxis=list(title=input$lopo_metric)
      )
    }
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
  })
  #######################################################################
  ### DATA SIZES 2
  #######################################################################
  output$data_sizes2= renderPrint({ 
    list_list_output = data()
    vec_traits = c()
    vec_validation_populations = c()
    vec_validation_population_sizes = c()
    for (list_output in list_list_output) {
      # list_output = list_list_output[[1]]
      if (list_output$TRAIT_NAME != input$lopo_trait) {next}
      if (is.na(list_output$YPRED_ACROSS_POP_LOPO[1])[1]) {next}
      idx = which(list_output$YPRED_ACROSS_POP_LOPO$pop_validation == input$lopo_pop)
      vec_pop_entry = unique(paste0(list_output$YPRED_ACROSS_POP_LOPO$pop_validation[idx], "---SEPARATOR---", list_output$YPRED_ACROSS_POP_LOPO$id[idx]))
      vec_pop_entry = unlist(lapply(strsplit(vec_pop_entry, "---SEPARATOR---"), FUN=function(x){x[1]}))
      vec_pop_sizes = table(vec_pop_entry)
      vec_traits = c(vec_traits, rep(list_output$TRAIT_NAME, times=length(vec_pop_sizes)))
      vec_validation_populations = c(vec_validation_populations, names(vec_pop_sizes))
      vec_validation_population_sizes = c(vec_validation_population_sizes, vec_pop_sizes)
    }
    df = data.frame(trait=vec_traits, population=vec_validation_populations, size=vec_validation_population_sizes)
    print("Population sizes:")
    print(df[(df$trait %in% input$lopo_trait) & (df$population %in% input$lopo_pop), ])
  })























  #######################################################################
  ### PAIRWISE: Scatter plot per trait-by-population combination and model
  #######################################################################
  # input = list(scat_pairwise_trait="dummyTrait_dummyPop", scat_pairwise_model="dummyModel1", scat_pairwise_pop_training="dummyPop1", scat_pairwise_pop_validation="dummyPop2")
  # input = list(scat_pairwise_trait="biomass_AUT_2023_20230302", scat_pairwise_model="Bayes_A", scat_pairwise_pop_training="DB-MS-31-22-001", scat_pairwise_pop_validation="DB-MS-31-22-002")
  output$plot_scatter_pairwise = renderPlotly({
    list_list_output = data()
    list_output = list_list_output[[input$scat_pairwise_trait]]
    if (!is.null(list_output$YPRED_ACROSS_POP_PAIRWISE)) {
      vec_idx = which(
        (list_output$YPRED_ACROSS_POP_PAIRWISE$model == input$scat_pairwise_model) & 
        (list_output$YPRED_ACROSS_POP_PAIRWISE$pop_training == input$scat_pairwise_pop_training) &
        (list_output$YPRED_ACROSS_POP_PAIRWISE$pop_validation == input$scat_pairwise_pop_validation))
      POP = paste0("trained: ", list_output$YPRED_ACROSS_POP_PAIRWISE$pop_training[vec_idx], "\nvalidated: ", list_output$YPRED_ACROSS_POP_PAIRWISE$pop_validation[vec_idx])
      df = data.frame(ID=list_output$YPRED_ACROSS_POP_PAIRWISE$id[vec_idx], POP=POP, Observed=list_output$YPRED_ACROSS_POP_PAIRWISE$y_true[vec_idx], Predicted=list_output$YPRED_ACROSS_POP_PAIRWISE$y_pred[vec_idx])
    } else {
      df = data.frame(ID="No leave-one-population-out data available", POP="", Observed=0, Predicted=0)
    }
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
        " | Model: ", input$scat_model
      )
    )
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
  })
  #######################################################################
  ### PAIRWISE: Distribution of observed and predicted trait values per trait-by-population combination and model
  #######################################################################
  output$plot_histogram_pairwise = renderPlotly({
    list_list_output = data()
    list_output = list_list_output[[input$scat_pairwise_trait]]
    if (!is.null(list_output$YPRED_ACROSS_POP_PAIRWISE)) {
      vec_idx = which(
        (list_output$YPRED_ACROSS_POP_PAIRWISE$model == input$scat_pairwise_model) & 
        (list_output$YPRED_ACROSS_POP_PAIRWISE$pop_training == input$scat_pairwise_pop_training) &
        (list_output$YPRED_ACROSS_POP_PAIRWISE$pop_validation == input$scat_pairwise_pop_validation))
      df = data.frame(ID=list_output$YPRED_ACROSS_POP_PAIRWISE$id[vec_idx], Observed=list_output$YPRED_ACROSS_POP_PAIRWISE$y_true[vec_idx], Predicted=list_output$YPRED_ACROSS_POP_PAIRWISE$y_pred[vec_idx])
    } else {
      df = data.frame(ID="No leave-one-population-out data available", Observed=0, Predicted=0)
    }
    p = plot_ly(alpha=0.5, type="histogram", nbinsx=input$pairwise_bins)
    p = p %>% add_histogram(data=df, x=~Observed, name="Observed")
    p = p %>% add_histogram(data=df, x=~Predicted, name="Predicted")
    p = p %>% layout(
      title="Phenotype distribution",
      barmode="overlay", 
      xaxis=list(title=list_output$TRAIT_NAME)
    )
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
  })
  #######################################################################
  ### PAIRWISE: bar plot (no replication as entire populations are used hence no subsampling performed)
  #######################################################################
  # input = list(pairwise_trait="dummyTrait", pairwise_model=paste0("dummyModel", 1:5), pairwise_pop_training="dummyPop1", pairwise_pop_validation=c("dummyPop1", "dummyPop2"))
  # list_list_output = list(
  #     dummyTrait_dummyPop=list(
  #       TRAIT_NAME="dummyTrait", 
  #       SUMMARY=data.frame(pop="dummyPop"),
  #       METRICS_WITHIN_POP=data.frame(trait="dummyTrait", pop="dummyPop", model="dummyModel", corr=rnorm(20)),
  #       METRICS_ACROSS_POP_LOPO=data.frame(trait="dummyTrait", pop_validation=paste0("dummyPop", 1:2), model=paste0("dummyModel", 1:5), corr=rnorm(10)),
  #       METRICS_ACROSS_POP_PAIRWISE=data.frame(trait="dummyTrait", pop_training=paste0("dummyPop", 1:2), pop_validation=paste0("dummyPop", 2:1), model=paste0("dummyModel", 1:5), corr=rnorm(10)),
  #       YPRED_WITHIN_POP=data.frame(id=rep(paste0("id", 1:10), times=5), model=rep(paste0("dummyModel", 1:5), each=10), y_true=rnorm(50), y_pred=rnorm(50)),
  #       YPRED_ACROSS_POP_LOPO=data.frame(id=paste0("id", 1:50), model=paste0("dummyModel", 1:5), y_true=rnorm(50), y_pred=rnorm(50), pop=paste0("dummyPop", 1:2)),
  #       YPRED_ACROSS_POP_PAIRWISE=data.frame(id=paste0("id", 1:50), model=paste0("dummyModel", 1:5), y_true=rnorm(50), y_pred=rnorm(50), pop_training=paste0("dummyPop", 1:2), pop_validation=paste0("dummyPop", 2:1)),
  #       GENOMIC_PREDICTION=data.frame(pop=paste0("dummyPop", 1:2), id=paste0("dummyEntry", 51:60), y_pred=rnorm(20), top_model="dummyModel", corr_from_kfoldcv=runif(20))
  #     )
  # )
  output$plot_bar_pairwise = renderPlotly({
    ### Define the grouping, i.e. trait and/or pop and/or model
    list_list_output = data()
    df_metrics_across_pop = NULL
    for (list_output in list_list_output) {
      if (list_output$TRAIT_NAME != input$pairwise_trait) {
        next
      }
      if (is.null(list_output$METRICS_ACROSS_POP_PAIRWISE)) {
        next
      }
      if (is.null(df_metrics_across_pop)) {
        df_metrics_across_pop = list_output$METRICS_ACROSS_POP_PAIRWISE
      } else {
        df_metrics_across_pop = rbind(df_metrics_across_pop, list_output$METRICS_ACROSS_POP_PAIRWISE)
      }
    }

    ### Revise to account for both training and validation populations
    if (is.null(df_metrics_across_pop)) {
      df = data.frame(ID="No leave-one-population-out data available", x=0, y=0)
      p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
    } else {
      ### Retain only the data set trained with the input population and validated in the defined validation populations
      vec_pairwise_pop_validation = input$pairwise_pop_validation[!(input$pairwise_pop_validation %in% input$pairwise_pop_training)]
      vec_idx = which(
        (df_metrics_across_pop$pop_training == input$pairwise_pop_training) &
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
        m = length(input$pairwise_model)
        M = matrix(NA, nrow=n, ncol=m)
        rownames(M) = vec_pairwise_pop_validation
        colnames(M) = input$pairwise_model
        for (i in 1:n) {
          for (j in 1:m) {
            vec_idx = which((df_metrics_across_pop$pop_validation == vec_pairwise_pop_validation[i]) & (df_metrics_across_pop$model == input$pairwise_model[j]))
            M[i, j] = eval(parse(text=paste0("df_metrics_across_pop$", input$pairwise_metric, "[vec_idx]")))
          }
        }
        vec_validation_pops = sort(unique(as.character(df_metrics_across_pop$pop_validation)))
        vec_models = sort(unique(as.character(df_metrics_across_pop$model)))
        if(!input$pairwise_group_model) {
          df = data.frame(vec_pairwise_pop_validation, M)
          colnames(df) = c("pop_validation", colnames(M))
          ### Plot
          eval(parse(text=paste0("p = plotly::plot_ly(data=df, x=~pop_validation, y=~", input$pairwise_model[1], ", type='bar', name='", input$pairwise_model[1], "')")))
          vec_metric = eval(parse(text=paste0("df$", input$pairwise_model[1])))
          if (length(input$pairwise_model) > 1) {
            for (i in 2:length(input$pairwise_model)) {
              eval(parse(text=paste0("p = p %>% plotly::add_trace(y=~", input$pairwise_model[i], ", name='", input$pairwise_model[i], "')")))
              vec_metric = c(vec_metric, eval(parse(text=paste0("df$", input$pairwise_model[i]))))
            }
          }
        } else {
          df = data.frame(input$pairwise_model, t(M))
          colnames(df) = c("model", rownames(M))
          ### Plot
          eval(parse(text=paste0("p = plotly::plot_ly(data=df, x=~model, y=~", vec_pairwise_pop_validation[1], ", type='bar', name='", vec_pairwise_pop_validation[1], "')")))
          vec_metric = eval(parse(text=paste0("df$", vec_pairwise_pop_validation[1])))
          if (length(vec_pairwise_pop_validation) > 1) {
            for (i in 2:length(vec_pairwise_pop_validation)) {
              eval(parse(text=paste0("p = p %>% plotly::add_trace(y=~", vec_pairwise_pop_validation[i], ", name='", vec_pairwise_pop_validation[i], "')")))
              vec_metric = c(vec_metric, eval(parse(text=paste0("df$", vec_pairwise_pop_validation[i]))))
            }
          }
        }
        p = p %>% plotly::layout(
          title=paste0("Trait: ", input$pairwise_trait, " | Mean ", input$pairwise_metric, " = ", round(mean(vec_metric, na.rm=TRUE), 4)), 
          barmode="group",
          yaxis=list(title=input$pairwise_metric)
        )
      }
    }
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
  })
  #######################################################################
  ### DATA SIZES 3
  #######################################################################
  output$data_sizes3= renderPrint({ 
    list_list_output = data()
    ### Training population size
    vec_training_population_sizes = c()
    for (list_output in list_list_output) {
      # list_output = list_list_output[[1]]
      if (list_output$TRAIT_NAME != input$pairwise_trait) {next}
      if (is.na(list_output$YPRED_ACROSS_POP_PAIRWISE[1])[1]) {next}
      for (pop_validation in input$pairwise_pop_validation) {
        vec_idx = which(
          (list_output$YPRED_ACROSS_POP_PAIRWISE$pop_training == pop_validation) & 
          (list_output$YPRED_ACROSS_POP_PAIRWISE$pop_validation == input$pairwise_pop_training))
        if (length(vec_idx) == 0) {next}
        vec_pop_entry = unique(paste0(list_output$YPRED_ACROSS_POP_PAIRWISE$pop_validation[vec_idx], "---SEPARATOR---", list_output$YPRED_ACROSS_POP_PAIRWISE$id[vec_idx]))
        vec_pop_entry = unlist(lapply(strsplit(vec_pop_entry, "---SEPARATOR---"), FUN=function(x){x[1]}))
        vec_training_population_sizes = c(vec_training_population_sizes, table(vec_pop_entry))
      }
    }
    print("Training population size:")
    print(unique(vec_training_population_sizes))
    ### Validation populations sizes
    vec_traits = c()
    vec_validation_populations = c()
    vec_validation_population_sizes = c()
    for (list_output in list_list_output) {
      # list_output = list_list_output[[1]]
      if (list_output$TRAIT_NAME != input$pairwise_trait) {next}
      if (is.na(list_output$YPRED_ACROSS_POP_PAIRWISE[1])[1]) {next}
      vec_idx = which(list_output$YPRED_ACROSS_POP_PAIRWISE$pop_training == input$pairwise_pop_training)
      vec_pop_entry = unique(paste0(list_output$YPRED_ACROSS_POP_PAIRWISE$pop_validation[vec_idx], "---SEPARATOR---", list_output$YPRED_ACROSS_POP_PAIRWISE$id[vec_idx]))
      vec_pop_entry = unlist(lapply(strsplit(vec_pop_entry, "---SEPARATOR---"), FUN=function(x){x[1]}))
      vec_pop_sizes = table(vec_pop_entry)
      vec_traits = c(vec_traits, rep(list_output$TRAIT_NAME, times=length(vec_pop_sizes)))
      vec_validation_populations = c(vec_validation_populations, names(vec_pop_sizes))
      vec_validation_population_sizes = c(vec_validation_population_sizes, vec_pop_sizes)
    }
    df = data.frame(trait=vec_traits, population=vec_validation_populations, size=vec_validation_population_sizes)
    print("Validation population sizes:")
    print(df[(df$trait %in% input$pairwise_trait) & (df$population %in% input$pairwise_pop_validation), ])
  })


















  #######################################################################
  ### Genomic predictions per se of entries with genotype data bu missing phenotype data
  #######################################################################
  output$plot_genomic_predictions_per_se = renderPlotly({
    list_list_output = data()
    df_genomic_predictions_per_se = NULL
    for (list_output in list_list_output) {
      if (list_output$TRAIT_NAME != input$gp_trait) {
        next
      }
      if (is.null(list_output$GENOMIC_PREDICTION)) {
        next
      }
      if (is.null(df_genomic_predictions_per_se)) {
        df_genomic_predictions_per_se = list_output$GENOMIC_PREDICTION
      } else {
        df_genomic_predictions_per_se = rbind(df_genomic_predictions_per_se, list_output$GENOMIC_PREDICTION)
      }
    }
    if (is.null(df_genomic_predictions_per_se)) {
      df = data.frame(ID="No genomic prediction per se performed", x=0, y=0)
      p = plotly::plot_ly(data=df, x=~x, y=~y, type="scatter", mode='markers', text=~ID)
    } else {
      p = plot_ly(alpha=0.5, type="histogram", nbinsx=input$gp_bins)
      for (pop in input$gp_pop) {
        idx = which(df_genomic_predictions_per_se$pop == pop)
        p = p %>% add_histogram(data=df_genomic_predictions_per_se[idx, ], x=~y_pred, name=pop)
      }
      p = p %>% layout(
        title="Genomic predictions distribution",
        barmode="overlay", 
        xaxis=list(title=input$gp_trait)
      )
    }
    p = p %>% config(toImageButtonOptions = list(format = "svg"))
  })
  #######################################################################
  ### Debugging messages
  #######################################################################
  # output$debug= renderPrint({ str(data()) })
  # output$debug= renderPrint({ 
  # list_list_output = data()
  #   name = grepl(paste0("^", input$scat_lopo_trait, "_"), names(list_list_output))
  #   # list_output = list_list_output[[name]]
  #   print(name)
  # })
  output$debug= renderPrint({
    list_list_output = data()
    print(str(list_list_output))
    # list_output = list_list_output[[input$scat_lopo_trait]]
    # if (!is.null(list_output$YPRED_ACROSS_POP)) {
    #   idx = which(list_output$YPRED_ACROSS_POP$model == input$scat_lopo_model)
    #   df = data.frame(ID=list_output$YPRED_ACROSS_POP$id[idx], POP=list_output$YPRED_ACROSS_POP$pop[idx], Observed=list_output$YPRED_ACROSS_POP$y_true[idx], Predicted=list_output$YPRED_ACROSS_POP$y_pred[idx])
    # } else {
    #   df = data.frame(ID="No leave-one-population-out data available", POP="", Observed=0, Predicted=0)
    # }
    # print(list_output$YPRED_ACROSS_POP)
  })
  # output$debug= renderPrint({ names(data()) })
  # output$debug= renderPrint({ input$scat_lopo_trait })
  # output$debug= renderPrint({ data()[[input$scat_lopo_trait]] })
  # output$debug= renderPrint({ length(data()) })
  # output$debug= renderPrint({ 
  #   list_list_output = data()
  #   df_genomic_predictions_per_se = NULL
  #   for (list_output in list_list_output) {
  #     if (list_output$TRAIT_NAME != input$gp_trait) {
  #       next
  #     }
  #     if (is.null(list_output$GENOMIC_PREDICTION)) {
  #       next
  #     }
  #     if (is.null(df_genomic_predictions_per_se)) {
  #       df_genomic_predictions_per_se = list_output$GENOMIC_PREDICTION
  #     } else {
  #       df_genomic_predictions_per_se = rbind(df_genomic_predictions_per_se, list_output$GENOMIC_PREDICTION)
  #     }
  #   }
  #   str(df_genomic_predictions_per_se)
  # })

}
####################################################################################################
### Serve the app
####################################################################################################
shinyApp(ui = ui, server = server)
