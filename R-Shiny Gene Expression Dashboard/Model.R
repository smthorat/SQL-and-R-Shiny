# Load necessary libraries
library(caret)
library(Biostrings) # Included for consistency, though not used here

# Set seed for reproducibility
set.seed(123)

# Simulate dataset
sequence_lengths <- rnorm(100, mean=500, sd=200) # Generate random sequence lengths
labels <- ifelse(sequence_lengths > 500, "Polyomer", "Monomer") # Assign labels based on length
data <- data.frame(SequenceLength = sequence_lengths, Label = as.factor(labels)) # Create dataframe

# Split dataset into training and testing sets
set.seed(123) # Ensure reproducibility
trainIndex <- createDataPartition(data$Label, p = .8, list = FALSE, times = 1) # Partition data
dataTrain <- data[trainIndex,] # Training set
dataTest  <- data[-trainIndex,] # Testing set

# Train logistic regression model
model <- glm(Label ~ SequenceLength, data = dataTrain, family = "binomial")

# Predict on test data
predictions <- predict(model, newdata = dataTest, type = "response") # Get predictions
predictedClasses <- ifelse(predictions > 0.5, "Polymer", "Monomer") # Convert to class labels

# Assess model accuracy
confMat <- confusionMatrix(as.factor(predictedClasses), dataTest$Label) # Generate confusion matrix
print(confMat) # Print confusion matrix

# Save the model
saveRDS(model, "protein_model.rds") # Save model for later use
