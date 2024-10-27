# This file contains backend process
# When we click on Run button we want to read the file by Seurat object
# First we need to check if the file is in correct format or not

# Install all required libraries 
# library(gdata)

# Install all required packages
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

# Function to read and validate Seurat object
load_seurat_obj <- function(path) {
  tryCatch({
    if (tolower(tools::file_ext(path)) != "rds") {
      return("Error: Invalid file format. Only .rds files are supported.")
    }
    
    obj <- readRDS(path)
    
    if (!inherits(obj, "Seurat")) {
      return("Error: The file is not a valid Seurat object.")
    }
    
    # Check if UMAP exists, if not, create it
    if (!"umap" %in% names(obj@reductions)) {
      # Run standard preprocessing if not already done
      if (!"pca" %in% names(obj@reductions)) {
        obj <- NormalizeData(obj)
        obj <- FindVariableFeatures(obj)
        obj <- ScaleData(obj)
        obj <- RunPCA(obj)
      }
      obj <- RunUMAP(obj, dims = 1:30)
    }
    
    return(obj)
    
  }, error = function(e) {
    return(paste("Error:", e$message))
  })
}

# Function to create UMAP plots from metadata
create_metadata_UMAP <- function(obj, col) {
  tryCatch({
    if (!"umap" %in% names(obj@reductions)) {
      return(ggplot() +
               theme_void() +
               geom_text(aes(x = 0.5, y = 0.5, 
                             label = "UMAP reduction not found"), 
                         size = 6, color = "red"))
    }
    
    if (col %in% c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
      # For numeric metadata
      col_df <- data.frame(
        UMAP_1 = obj@reductions$umap@cell.embeddings[,1],
        UMAP_2 = obj@reductions$umap@cell.embeddings[,2],
        data = obj@meta.data[[col]]
      )
      
      p <- ggplot(col_df, aes(x = UMAP_1, y = UMAP_2, color = log10(data))) +
        geom_point(size = 0.1) +
        scale_color_gradientn(colours = rainbow(7)) +
        theme_minimal() +
        labs(title = col,
             color = paste("log10(", col, ")")) +
        theme(legend.position = "right",
              plot.title = element_text(hjust = 0.5))
      
    } else if (col %in% colnames(obj@meta.data)) {
      # For categorical metadata
      p <- DimPlot(obj, 
                   reduction = "umap",
                   group.by = col,
                   pt.size = 0.1,
                   label = TRUE,
                   label.size = 4) +
        theme_minimal() +
        labs(title = col)
      
    } else {
      p <- ggplot() +
        theme_void() +
        geom_text(aes(x = 0.5, y = 0.5, 
                      label = "Column not found in metadata"),
                  size = 6, color = "red")
    }
    
    return(p)
    
  }, error = function(e) {
    return(ggplot() +
             theme_void() +
             geom_text(aes(x = 0.5, y = 0.5, 
                           label = paste("Error:", e$message)), 
                       size = 6, color = "red"))
  })
}

# Function to create feature plots
create_feature_plot <- function(obj, gene) {
  tryCatch({
    if (!"umap" %in% names(obj@reductions)) {
      return(ggplot() +
               theme_void() +
               geom_text(aes(x = 0.5, y = 0.5, 
                             label = "UMAP reduction not found"), 
                         size = 6, color = "red"))
    }
    
    if (gene %in% rownames(obj)) {
      p <- FeaturePlot(obj,
                       features = gene,
                       reduction = "umap",
                       pt.size = 0.1) +
        theme_minimal() +
        scale_color_gradientn(colours = c("lightgrey", "red")) +
        labs(title = paste("Expression of", gene))
      
    } else {
      p <- ggplot() +
        theme_void() +
        geom_text(aes(x = 0.5, y = 0.5, 
                      label = "Gene not found in dataset"), 
                  size = 6, color = "red")
    }
    
    return(p)
    
  }, error = function(e) {
    return(ggplot() +
             theme_void() +
             geom_text(aes(x = 0.5, y = 0.5, 
                           label = paste("Error:", e$message)), 
                       size = 6, color = "red"))
  })
}