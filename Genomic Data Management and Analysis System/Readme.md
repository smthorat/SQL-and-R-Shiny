{\rtf1\ansi\ansicpg1252\cocoartf2818
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Genomic Data Management and Analysis System\
\
This project is about building a database to manage and analyze genomic data, specifically for human genes (*Homo sapiens*). The goal is to take raw data, clean it, store it in a database, run SQL queries, and create visualizations to better understand the data.\
\
## Project Structure\
\
Here is a quick rundown of the files in the project:\
\
- **Data_Cleaning.ipynb**: This is a Jupyter notebook used for cleaning the raw gene data. It includes steps like dealing with missing values, dropping columns that are not needed, and formatting the text properly.\
- **Visualization.ipynb**: This Jupyter notebook is used to create visualizations. It uses Python libraries like `matplotlib` and `seaborn` to make charts and graphs for analyzing the data.\
- **Homo_sapiens_Gene.tsv**: This is the original dataset in TSV format. It contains the raw data about human genes before any cleaning.\
- **cleaned_gene_data.csv**: This is the cleaned version of the dataset, saved as a CSV file. It is used to import the data into the database.\
- **code.sql**: This file has all the SQL commands used in the project, including creating tables, importing data, running queries, and setting up triggers and stored procedures.\
- **10_most_common_genes.png**: A bar chart that shows the 10 most common genes in the dataset.\
- **Gene_Description_Length.png**: A plot that displays the distribution of gene description lengths.\
- **Gene_distribution.png**: A graph that illustrates how genes are distributed across different categories.\
- **Gene_type_proportion.png**: A pie chart that shows the proportion of different gene types in the dataset.\
\
## What I Did\
\
1. **Data Collection**\
   - The starting point was the raw data file (`Homo_sapiens_Gene.tsv`). It is a TSV file downloaded from a reliable source, such as NCBI, that contains information about human genes.\
\
2. **Data Cleaning**\
   - I cleaned the data using `Data_Cleaning.ipynb`. This involved removing rows with missing values, dropping columns that were not needed, formatting the text to be consistent, and trimming extra spaces.\
   - After cleaning, the data was saved as `cleaned_gene_data.csv`, which is ready to be imported into the database.\
\
3. **Database Setup**\
   - I set up the database schema and tables using the `code.sql` file. This file has all the commands needed to create the tables and load the cleaned data.\
   - The cleaned data from `cleaned_gene_data.csv` was imported into the database.\
\
4. **Data Analysis**\
   - I wrote SQL queries to do things like retrieving specific gene information, filtering by gene type, and counting occurrences.\
   - I also added advanced features like stored procedures, triggers, and views to make the database more useful and automated.\
\
5. **Visualization**\
   - The `Visualization.ipynb` notebook was used to make different charts to help understand the dataset better.\
   - The visualizations created include `10_most_common_genes.png`, `Gene_Description_Length.png`, `Gene_distribution.png`, and `Gene_type_proportion.png`.\
\
## Tools and Technologies Used\
\
- **SQL**: Used to create the database, run queries, and manage data.\
- **Python**: Used for data cleaning and making visualizations. I mainly used `pandas`, `matplotlib`, and `seaborn`.\
- **Jupyter Notebooks**: Used for interactive data cleaning and visualization.\
- **Data Formats**: The project used both TSV and CSV files.\
\
## How to Run the Project\
\
1. **Set Up the Database**\
   - Install a SQL database system, like MySQL, PostgreSQL, or SQLite.\
   - Run the commands in `code.sql` to create the database and tables.\
   - Import the `cleaned_gene_data.csv` file into the database.\
\
2. **Run the Jupyter Notebooks**\
   - Start with `Data_Cleaning.ipynb` to see the steps used for cleaning the data.\
   - Then open `Visualization.ipynb` to check out the visualizations.\
\
3. **Query the Database**\
   - Use the queries in `code.sql` to explore and analyze the data.\
\
## Future Plans\
\
- Add data from other species.\
- Build a web interface to make it easier to access the data.\
- Try machine learning techniques to predict gene functions.\
\
## Notes\
\
- Make sure you have all the required Python libraries installed before running the notebooks.\
- The visualizations are optional but help a lot in understanding the data.}