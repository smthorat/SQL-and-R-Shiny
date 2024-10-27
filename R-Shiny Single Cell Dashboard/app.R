# Load libraries
# Load required libraries
library(tools)
library(Seurat)
library(dplyr)
library(shinyjs)
library(shiny)
library(DT)
library(shinydashboard)
library(shinydashboardPlus)
library(ggplot2)
library(shinybusy)
library(glue)
library(markdown)
library(ggthemes)

# UI Definition
ui <- dashboardPage(
  dashboardHeader(title = "scRNAseq Analysis"),
  dashboardSidebar(
    tags$head(tags$style(HTML(".skin-blue .main-header .sidebar-toggle {display: none;}"))),
    sidebarMenu(
      id = 'tab',
      useShinyjs(),
      menuItem("Home Page", tabName = "home", icon = icon("list")),
      menuItem("scRNAseq Analyzer", tabName = "input", icon = icon("edit")),
      conditionalPanel(
        condition = "input.tab == 'input'",
        div(
          style = "padding: 15px;",
          fileInput("file", "Upload File", 
                    multiple = FALSE,
                    accept = c('.rds')),
          actionButton("reset", "Reset", 
                       icon = icon("undo"), 
                       style = "color: #fff; background-color: #dc3545; width: 100%; margin-bottom: 10px;"),
          actionButton("run", "Run", 
                       icon = icon("play"), 
                       style = "color: #fff; background-color: #28a745; width: 100%;")
        )
      )
    )
  ),
  dashboardBody(
    useShinyjs(),
    add_busy_spinner(spin = "fading-circle"),
    tabItems(
      tabItem(
        tabName = "input",
        tabsetPanel(
          id = 'main_tabs',
          tabPanel(
            "Instructions",
            h2("Instructions for Using the scRNAseq Analysis App"),
            tags$div(
              style = "padding: 20px;",
              tags$ol(
                tags$li("Upload your Seurat object (.rds file)"),
                tags$li("Click the 'Run' button to generate visualizations"),
                tags$li("Use the tabs to explore different analyses"),
                tags$li("Use the 'Reset' button to clear and start over")
              )
            )
          )
        )
      ),
      tabItem(
        tabName = "home",
        div(
          style = "padding: 20px;",
          tags$h1(
            style = "text-align: center; color: #2c3e50;",
            HTML("<u>Welcome to The scRNAseq Seurat Analysis RShiny App</u>")
          ),
          tags$p(
            style = "text-align: center; font-size: 18px; margin-top: 20px;",
            "This application allows you to analyze single-cell RNA sequencing data using the Seurat package."
          )
        )
      )
    )
  )
)

# Server Definition
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 300*1024^2)
  
  # Reactive values to store the Seurat object
  values <- reactiveValues(
    seurat_obj = NULL
  )
  
  # Function to process Seurat object
  process_seurat <- function(obj) {
    tryCatch({
      # Basic preprocessing if not already done
      if (!"RNA" %in% names(obj@assays)) {
        stop("RNA assay not found in the Seurat object")
      }
      
      # Check if normalization is needed
      if (!"data" %in% names(obj@assays$RNA)) {
        obj <- NormalizeData(obj)
      }
      
      # Find variable features if not already done
      if (length(VariableFeatures(obj)) == 0) {
        obj <- FindVariableFeatures(obj)
      }
      
      # Scale data if not already done
      if (!"scale.data" %in% names(obj@assays$RNA)) {
        obj <- ScaleData(obj)
      }
      
      # Run PCA if not present
      if (!"pca" %in% names(obj@reductions)) {
        obj <- RunPCA(obj)
      }
      
      # Run UMAP if not present
      if (!"umap" %in% names(obj@reductions)) {
        obj <- RunUMAP(obj, dims = 1:30)
      }
      
      return(obj)
    }, error = function(e) {
      stop(paste("Error in processing:", e$message))
    })
  }
  
  # Function to load Seurat object
  load_seurat_obj <- function(filepath) {
    tryCatch({
      obj <- readRDS(filepath)
      if (!inherits(obj, "Seurat")) {
        return(c("Error: Uploaded file is not a valid Seurat object"))
      }
      # Process the object
      obj <- process_seurat(obj)
      return(obj)
    }, error = function(e) {
      return(c("Error loading/processing file:", e$message))
    })
  }
  
  # Disable run button initially
  observe({
    shinyjs::toggleState("run", !is.null(input$file))
  })
  
  # Handle file upload and analysis
  observeEvent(input$run, {
    shinyjs::disable("run")
    
    # Show loading spinner with progress
    withProgress(message = 'Processing data...', {
      
      # Remove existing tabs if they exist
      for (tab in c("UMAP", "Gene Expression")) {
        if (tab %in% names(input$main_tabs)) {
          removeTab("main_tabs", tab)
        }
      }
      
      # Validate file input
      req(input$file)
      
      incProgress(0.2, detail = "Loading Seurat object...")
      
      # Load and process the Seurat object
      values$seurat_obj <- load_seurat_obj(input$file$datapath)
      
      if (is.character(values$seurat_obj)) {
        # Show error modal if loading failed
        showModal(modalDialog(
          title = "Error with file",
          HTML(paste("<h5>There is an error with the file you uploaded:</h5><br>",
                     paste(values$seurat_obj, collapse = "<br><br>"))),
          easyClose = TRUE
        ))
      } else {
        incProgress(0.4, detail = "Creating visualizations...")
        
        # Add UMAP tab
        insertTab(
          inputId = "main_tabs",
          tabPanel(
            "UMAP",
            div(
              style = "margin: 20px;",
              selectizeInput("metadata_col", "Select Metadata Column", 
                             choices = c("Default" = "seurat_clusters", 
                                         colnames(values$seurat_obj@meta.data))),
              plotOutput("umap", height = "600px")
            )
          ),
          select = TRUE
        )
        
        # Add Gene Expression tab
        insertTab(
          inputId = "main_tabs",
          tabPanel(
            "Gene Expression",
            div(
              style = "margin: 20px;",
              selectizeInput("gene", "Select Gene", 
                             choices = rownames(values$seurat_obj), 
                             selected = NULL),
              plotOutput("featurePlot", height = "600px")
            )
          )
        )
        
        # Generate UMAP plot
        output$umap <- renderPlot({
          req(values$seurat_obj, input$metadata_col)
          DimPlot(values$seurat_obj, 
                  reduction = "umap", 
                  group.by = input$metadata_col,
                  label = TRUE) +
            theme_minimal() +
            ggtitle(paste("UMAP colored by", input$metadata_col))
        })
        
        # Generate Feature plot
        output$featurePlot <- renderPlot({
          req(values$seurat_obj, input$gene)
          FeaturePlot(values$seurat_obj, 
                      features = input$gene,
                      reduction = "umap") +
            theme_minimal() +
            ggtitle(paste("Expression of", input$gene))
        })
      }
      
      incProgress(0.4, detail = "Finishing up...")
    })
    
    shinyjs::enable("run")
  })
  
  # Reset functionality
  observeEvent(input$reset, {
    shinyjs::reset("file")
    values$seurat_obj <- NULL
    
    for (tab in c("UMAP", "Gene Expression")) {
      if (tab %in% names(input$main_tabs)) {
        removeTab("main_tabs", tab)
      }
    }
    
    shinyjs::disable("run")
  })
}

# Run the application
shinyApp(ui = ui, server = server)