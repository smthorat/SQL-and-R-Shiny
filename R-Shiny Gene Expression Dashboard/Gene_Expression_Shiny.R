library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
gene_data <- read.delim("/Volumes/Jagannath/R_Programming/GSE183947_long_format.txt")

# Define UI
ui <- fluidPage(
  titlePanel("Gene Expression Dashboard"),
  sidebarLayout(
    sidebarPanel(
      h4("Gene Expression Analysis"),
      selectInput("selectedGene", "Select Gene:", choices = unique(gene_data$gene))
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Table", DTOutput("geneTable")),
        tabPanel("Barplot", plotOutput("barPlot")), # this is the dashboard of shiny app. 
        tabPanel("Density", plotOutput("densityPlot")),
        tabPanel("Boxplot/Violin", plotOutput("boxPlot")),
        tabPanel("Scatterplot", plotOutput("scatterPlot"))
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  output$geneTable <- renderDT({
    datatable(gene_data)
  })
  
  output$barPlot <- renderPlot({
    filtered_data <- gene_data %>% filter(gene == input$selectedGene)
    
    ggplot(filtered_data, aes(x = samples, y = FPKM, fill = tissue)) +
      geom_col() +
      labs(title = paste("Barplot of FPKM Values for", input$selectedGene),
           x = "Sample",
           y = "FPKM")
  })
  
  output$densityPlot <- renderPlot({
    filtered_data <- gene_data %>% filter(gene == input$selectedGene)
    
    ggplot(filtered_data, aes(x = FPKM, fill = tissue)) +
      geom_density(alpha = 0.3) +
      labs(title = paste("Density Plot of FPKM Values for", input$selectedGene),
           x = "FPKM",
           y = "Density")
  })
  
  output$boxPlot <- renderPlot({
    filtered_data <- gene_data %>% filter(gene == input$selectedGene)
    
    ggplot(filtered_data, aes(x = tissue, y = FPKM)) +
      geom_violin() +
      labs(title = paste("Boxplot/Violin Plot of FPKK Values for", input$selectedGene),
           x = "Tissue",
           y = "FPKM")
  })
  
  output$scatterPlot <- renderPlot({
    filtered_data <- gene_data %>% filter(gene == input$selectedGene)
    
    if (nrow(filtered_data) > 1) { # Ensure there's enough data to plot
      ggplot(filtered_data, aes(x = samples, y = FPKM, color = tissue, label = tissue)) +
        geom_point() +
        geom_text(vjust = -1) +
        labs(title = paste("Scatterplot of FPKM Values for", input$selectedGene),
             x = "Sample",
             y = "FPKM")
    }
  })
}

# Run the app
shinyApp(ui = ui, server = server)
