library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming gene_data is already read from the specified path
gene_data <- read.delim("/Volumes/Jagannath/R_Programming/GSE183947_long_format.txt")

# Check if gene_data is a data frame
if (!is.data.frame(gene_data)) {
  print("Error: gene_data is not a data frame. Please check the file and path.")
  print(str(gene_data)) # This will print the structure of gene_data
} else {
  print("gene_data loaded successfully and is a data frame.")
  print(head(gene_data)) # This prints the first few rows of the data frame
}

# Define UI
ui <- fluidPage(
  titlePanel("Updated Gene Expression Dashboard"),
  sidebarLayout(
    sidebarPanel(
      h4("Filters"),
      selectInput("selectedGene", "Select Gene:", choices = unique(gene_data$gene)),
      sliderInput("fpkmRange", "FPKM Value Range:",
                  min = min(gene_data$FPKM), max = max(gene_data$FPKM),
                  value = c(min(gene_data$FPKM), max(gene_data$FPKM))),
      selectInput("selectedTissue", "Select Tissue:", choices = c("All", unique(gene_data$tissue))),
      actionButton("update", "Update View")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Table View", DTOutput("geneTable")),
        tabPanel("Histogram", plotOutput("histPlot")),
        tabPanel("Line Plot", plotOutput("linePlot")),
        tabPanel("Scatter Plot", plotOutput("scatterPlot")) # New scatter plot tab
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  filtered_data <- reactive({
    fd <- gene_data %>%
      filter(FPKM >= input$fpkmRange[1] & FPKM <= input$fpkmRange[2])
    if (input$selectedTissue != "All") {
      fd <- fd %>% filter(tissue == input$selectedTissue)
    }
    fd
  })
  
  output$geneTable <- renderDT({
    datatable(filtered_data())
  })
  
  output$histPlot <- renderPlot({
    valid_data <- filtered_data() %>% filter(!is.na(FPKM) & FPKM != Inf & FPKM != -Inf)
    valid_data$tissue <- factor(valid_data$tissue, levels = c("normal breast tissue", "breast tumor"))
    ggplot(valid_data, aes(x = FPKM, fill = tissue)) +
      geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
      labs(title = "Distribution of FPKM Values", x = "FPKM", y = "Count") +
      scale_fill_manual(values = c("normal breast tissue" = "skyblue", "breast tumor" = "salmon"))
  })
  
  output$linePlot <- renderPlot({
    gene_specific_data <- filtered_data() %>% filter(gene == input$selectedGene)
    ggplot(gene_specific_data, aes(x = samples, y = FPKM, group = 1)) +
      geom_line() +
      geom_point(aes(color = tissue)) +
      labs(title = paste("Trend of FPKM Values for", input$selectedGene), x = "Sample", y = "FPKM")
  })
  
  # Scatter Plot
  output$scatterPlot <- renderPlot({
    scatter_data <- filtered_data() %>% filter(gene == input$selectedGene)
    ggplot(scatter_data, aes(x = samples, y = FPKM, color = tissue)) +
      geom_point() +
      labs(title = paste("FPKM Values Scatter Plot for", input$selectedGene), x = "Sample", y = "FPKM")
  })
}

# Run the app
shinyApp(ui = ui, server = server)
