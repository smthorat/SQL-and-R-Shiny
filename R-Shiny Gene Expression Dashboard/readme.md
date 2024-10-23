## Overview: 
This Shiny application is designed to classify protein sequences as either Polygomer or Monomer 
based on the sequence length. It provides a user-friendly interface to input the sequence length 
and view the prediction, as well as to visualize the model's performance in terms of accuracy.

## Installation
To run this application, you need to have R installed on your system.
install.packages("shiny")
install.packages("shinydashboard")
install.packages("caret")
install.packages("Biostrings")
install.packages("ggplot2")

## Running the Application
Model.R
Shiny.R

## Application Structure
Predict: Here, users can input a protein sequence length and click on the "Predict" button to receive 
a prediction stating whether the protein is a Polygomer or Monomer based on the model's inference.

Model Performance: This section displays a graphical representation of the model's accuracy 
in classifying the protein sequences.