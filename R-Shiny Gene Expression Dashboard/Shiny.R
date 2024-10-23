library(shiny)
library(shinydashboard)
library(caret)
library(Biostrings) # Assuming it's necessary elsewhere
library(ggplot2)

# Assuming you've saved the dataTest dataframe and the model externally
# Load them here if necessary

# Define UI
ui <- dashboardPage(
    dashboardHeader(title = "Protein Sequence"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Predict", tabName = "predict", icon = icon("dna")),
            menuItem("Model Performance", tabName = "performance", icon = icon("chart-line"))
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName = "predict",
                    fluidRow(
                        box(title = "Enter Sequence Length", status = "primary", solidHeader = TRUE,
                            numericInput("sequence_length", label = "Sequence Length", value = 500)),
                        box(title = "Prediction", status = "info", solidHeader = TRUE,
                            actionButton("predict", "Predict"),
                            verbatimTextOutput("prediction_result"))
                    )),
            tabItem(tabName = "performance",
                    fluidRow(
                        box(title = "Model Accuracy", status = "warning", solidHeader = TRUE,
                            plotOutput("accuracyPlot"))
                    ))
        )
    )
)

# Define server logic
server <- function(input, output) {
  observeEvent(input$predict, {
    tryCatch({
      sequence_length <- input$sequence_length
      # Ensure model is loaded
      model <- readRDS("protein_model.rds")
      prediction_prob <- predict(model, newdata = data.frame(SequenceLength = sequence_length), type = "response")
      predicted_class <- ifelse(prediction_prob > 0.5, "Polymer", "Monomer")
      
      output$prediction_result <- renderText({
        paste("Based on the sequence length, the protein is predicted to be a:", predicted_class)
      })
    }, error = function(e) {
      output$prediction_result <- renderText({
        paste("Error in prediction:", e$message)
      })
    })
  })

    output$accuracyPlot <- renderPlot({
        # Load or recalculate the necessary data for the confusion matrix
        # For demonstration, let's assume you've loaded or recalculated the dataTest and predictions

        # Pre-calculate or load the confusion matrix and accuracy if needed
        # confMat <- confusionMatrix(as.factor(predictedClasses), dataTest$Label)
        # accuracy <- confMat$overall['Accuracy']
        # For illustration, assuming an accuracy value
        accuracy <- 0.9 # Placeholder for actual accuracy calculation
        df <- data.frame(Metric = "Accuracy", Value = accuracy)
        ggplot(df, aes(x = Metric, y = Value)) + geom_col() +
            ylim(0, 1) + labs(y = "Value", title = "Model Accuracy") +
            theme_minimal()
    })
}

# Run the application
shinyApp(ui = ui, server = server)
