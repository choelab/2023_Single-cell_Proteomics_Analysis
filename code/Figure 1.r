###############################
# title : Figure 1
# author : Jaemyung, Jang (piloter2@kbri.re.kr)
# kenel : R 4.3.0
# Date : Dec. 18, 2023
###############################

source(paste0(path, "/code/proteomic_function_20231218.r"))
# Load necessary libraries
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggVennDiagram)
library(pbapply)
library(stringr)
library(reshape2)

# Define the path to the data directory
path_to_data <- file.path("/", "choelab", "SCP")
raw_proteins_pattern <- "_Proteins.txt$"

# Define Functions
# -----------------------------------------------------------
# Function to list files based on pattern
list_files <- function(path, pattern) {
  list.files(list.dirs(file.path(path, "proteomicDatafromPD")), pattern = pattern, full.names = TRUE)
}

# Function to clean column names
clean_colnames <- function(df) {
  df %>% 
    setNames(str_replace_all(colnames(df), " ", ".") %>% 
             str_replace_all("-", "_"))
}

# Function to read and clean data
read_and_clean_data <- function(filenames) {
  lapply(filenames, function(file) fread(file) %>% clean_colnames)
}

# Function to extract specific data based on column pattern
extract_data_subset <- function(data, pattern) {
  select_cols <- grep(pattern, colnames(data), value = TRUE)
  data[, c("Accession", select_cols), with = FALSE]
}

# Function to get unique non-missing Accession values
get_unique_accessions <- function(data) {
  unique(unlist(lapply(data, function(x) x[!is.na(x[, 2:ncol(x)]), Accession])))
}

# Function to create a Venn Diagram
create_venn_diagram <- function(lists) {
  venn_data <- lapply(lists, get_unique_accessions)
  ggVennDiagram(venn_data) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggtitle("Venn Diagram")
}

# Function to fetch data from UniProt
fetch_uniprot_data <- function(ids) {
  # Fetching data logic here
}

# Function to calculate protein length
calculate_protein_length <- function(data, sample_name) {
  # Protein length calculation logic here
}

# Main Script
# -----------------------------------------------------------
# Load Raw Data
raw_proteins_files <- list_files(path_to_data, raw_proteins_pattern)
protein_data <- read_and_clean_data(raw_proteins_files)

# Process Data
protein_counts <- lapply(protein_data, function(data) {
  extract_data_subset(data, "Abundances.")
})

# Figure 1h. Create Venn Diagram
venn_diagram_data <- list(
  "DDM" = extract_data_subset(protein_counts[[1]], "Sample.DDM$"),
  "SB3-12" = extract_data_subset(protein_counts[[1]], "1_40772$"),
  "Q" = extract_data_subset(protein_counts[[1]], "Point3_Q$")
)
venn_diagram <- create_venn_diagram(venn_diagram_data)

# Fetch Data from UniProt
uniprot_results <- lapply(protein_data[[1]]$Accession, fetch_uniprot_data)
uniprot_df <- data.frame('Accession' = protein_data[[1]]$Accession, 'Gene.names' = unlist(uniprot_results))

# Figure 1i. Protein Length Analysis
protein_length_data <- calculate_protein_length(protein_data[[1]], "SampleName")
