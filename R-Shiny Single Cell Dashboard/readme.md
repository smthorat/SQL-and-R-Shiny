
# Single Cell RNA-seq Analysis Dashboard

Welcome to the Single Cell RNA-seq Analysis Dashboard! This interactive RShiny app is designed to make it easier for you to analyze single-cell RNA sequencing (scRNAseq) data using the Seurat package. 

## Overview {#overview}

The dashboard simplifies scRNAseq data analysis by allowing you to upload `.rds` files with Seurat objects, process them, and visualize data through UMAP and gene expression plots. This project helps researchers and students quickly explore single-cell data with minimal coding required.

## Features {#features}

- **Simple Data Upload**: Drag and drop your `.rds` files directly into the app.
- **Automated Data Processing**: Includes common preprocessing steps such as normalization, scaling, PCA, and UMAP.
- **Visual Exploration**: View your data using UMAP and gene expression plots with options for customization.
- **Interactive Dashboard**: Use tabs for easy navigation, get real-time processing updates, and reset to start over when needed.

## Installation {#installation}

### Prerequisites

You'll need **R** and **RStudio** installed on your computer. The dashboard uses various R packages like Seurat and Shiny. To install all required packages, run the following code:

```{r install-packages, eval=FALSE}
install.packages(c("tools", "Seurat", "dplyr", "shinyjs", "shiny", "DT", 
                  "shinydashboard", "shinydashboardPlus", "ggplot2", 
                  "shinybusy", "glue", "markdown", "ggthemes"))
```

### Running the Dashboard

1. Clone the repository to your local machine:

```{bash, eval=FALSE}
git clone https://github.com/yourusername/scRNAseq-analysis-dashboard.git
cd scRNAseq-analysis-dashboard
```

2. Launch the app by opening `app.R` in RStudio and clicking "Run App" or running the following in your R console:

```{r run-app, eval=FALSE}
shiny::runApp("app.R")
```

## Using the Dashboard {#using-the-dashboard}

### 1. Uploading Your Data

To start, prepare your scRNAseq data as a Seurat object saved in an `.rds` file. Once ready:

1. Open the dashboard and navigate to the scRNAseq Analyzer tab on the sidebar.
2. Click Upload File and select your `.rds` file.
3. After the file uploads, the Run button will become clickable.

### 2. Running the Analysis

Once your file is uploaded:

1. Click the Run button to start processing your data.
2. The dashboard will go through the necessary steps to preprocess and analyze the data, including creating UMAP and gene expression visualizations.
3. A progress bar will appear to show the current status of your analysis.

### 3. Exploring Results

After processing, you'll have access to several visualization tabs:

- **UMAP Tab**: This tab shows the UMAP (Uniform Manifold Approximation and Projection) plot of your data, where each point represents a cell. Use the Metadata Column dropdown menu to color the UMAP plot based on different attributes (e.g., cell types or gene expression levels).
- **Gene Expression Tab**: In this tab, you can view expression levels for specific genes across your dataset. Choose a gene from the Select Gene dropdown menu to see how it is expressed across cells in the UMAP plot.

### 4. Resetting the Dashboard

If you want to start over:

1. Press the Reset button, which clears uploaded files, removes generated plots, and resets all fields.
2. This allows you to upload a new file and run a fresh analysis without needing to reload the app.

## Contributing {#contributing}

We'd love for you to contribute to this project! Here's how you can help:

1. Fork the repository on GitHub.
2. Clone your forked repository to make local changes:

```{bash, eval=FALSE}
git clone https://github.com/yourusername/scRNAseq-analysis-dashboard.git
```

3. Create a new branch for your changes:

```{bash, eval=FALSE}
git checkout -b feature/your-feature-name
```

4. Make and test your changes locally. When ready, commit your changes and push them to GitHub.
5. Submit a Pull Request with a short description of what you added or fixed.

## License {#license}

This project is licensed under the MIT License. See the LICENSE file for more details.

```{r session-info, echo=FALSE}
sessionInfo()
```

# Result<img width="1035" alt="Image   1 26 12 PM" src="https://github.com/user-attachments/assets/6b7f99e1-a2aa-4136-a986-97b07865dae8">
<img width="964" alt="Image   1 26 27 PM" src="https://github.com/user-attachments/assets/0edd83b3-732f-4aae-ad10-9893975ac87b">
<img width="966" alt="Image   1 26 44 PM" src="https://github.com/user-attachments/assets/fa01ce2d-707d-495e-84e1-b86fa837a821">
<img width="1009" alt="Image   1 26 59 PM" src="https://github.com/user-attachments/assets/0c0628f9-20e6-4040-9281-735e56bad63e">
<img width="961" alt="Image   1 27 15 PM" src="https://github.com/user-attachments/assets/400960b6-1a60-4837-8d64-3e4b706f9a1f">
<img width="971" alt="Image   1 28 25 PM" src="https://github.com/user-attachments/assets/c3637c5c-2b26-47c9-b40c-9acc8979d3ff">


