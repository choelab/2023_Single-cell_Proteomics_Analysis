###############################
# title : Figure 3
# author : Jaemyung, Jang (piloter2@kbri.re.kr)
# kenel : R 4.3.0
# Date : Dec. 18, 2023
###############################

library(data.table)
library(dplyr)
library(ggVennDiagram)

# Define the path to the data directory
path <- file.path("/", "choelab", "SCP")

# Source the external R script containing additional functions
source(paste0(path, "/code/proteomic_function_20231218.r"))

# Define file patterns for raw data
rawPD_pattern <- "_Proteins.txt$"

# Get the list of files matching the pattern
rawPD_files <- list.files(path, pattern = rawPD_pattern, full.names = TRUE)

# Function to clean column names
clean_colnames <- function(df) {
  colnames(df) <- gsub(" ", ".", colnames(df))
  colnames(df) <- gsub("-", "_", colnames(df))
  df
}

# Function to read and clean data from files
read_and_clean_data <- function(filenames) {
  lapply(filenames, function(file) {
    fread(file) %>% clean_colnames()
  })
}

# Read protein data
protein_data_list <- read_and_clean_data(rawPD_files)

# Extract counts data for specific cell lines
extract_counts <- function(data_list, pattern) {
  lapply(data_list, function(data) {
    data %>% 
      select(Accession, matches(pattern))
  })
}

# Extract normalized abundances for selected cell lines
abundance_patterns <- c("1_LNCaP3000$", "2_L100$", "3_C4$", "4_U118$")
names(abundance_patterns) <- c("LNCaP3000", "LNCaP", "C4", "U118")

count_scp <- lapply(abundance_patterns, function(pattern) {
  extract_counts(protein_data_list, pattern)
})

# Function to get unique Accession numbers for each cell line
get_unique_accessions <- function(count_data) {
  lapply(count_data, function(data_frame_list) {
    sapply(data_frame_list, function(data_frame) {
      unique(na.omit(data_frame$Accession))
    })
  })
}

# Get unique Accession numbers
unique_accessions <- get_unique_accessions(count_scp)

venn_data <- lapply(unique_accessions, unlist)
# Figure 3c. Create a Venn diagram for the selected cell lines
Figure3c <- ggVennDiagram(venn_data) + scale_fill_gradient(low = "blue", high = "red")

# Helper function to get gene names from Uniprot
fetch_gene_names_from_uniprot <- function(accession_ids) {
  sapply(accession_ids, function(id) {
    response <- httr::GET(url = paste0('https://rest.uniprot.org/uniprotkb/', id))
    content <- httr::content(response, "text")
    gene_name <- str_extract(content, "(?<=<geneName>).*?(?=</geneName>)")
    return(gene_name)
  })
}

# Extract gene names from descriptions and Uniprot
get_gene_names <- function(data_from_PD, index) {
  data <- data_from_PD[[index]]
  accessions <- data$Accession
  descriptions <- data$Description

  # Extract gene names from descriptions
  gene_names_from_descriptions <- sapply(descriptions, function(desc) {
    str_extract(desc, "(?<=GN=)[^ ]+")
  })

  # Fetch gene names from Uniprot
  gene_names_from_uniprot <- fetch_gene_names_from_uniprot(accessions)
  
  # Merge and resolve conflicts
  gene_names <- ifelse(is.na(gene_names_from_uniprot), gene_names_from_descriptions, gene_names_from_uniprot)
  return(data.frame(Accession = accessions, Gene.names = gene_names))
}

# Process for the specified indexes
gene_names_list <- lapply(1:4, function(k) get_gene_names(protein_data_list, k))

# Create SummarizedExperiment objects
create_summarized_experiment <- function(data, design, count_data_index) {
  # Merge count data with gene names
  merged_data <- merge(count_data_index, data, by = "Accession")
}

# Run for all the specified indexes
summarized_experiments <- lapply(c(1:4), function(k) {
  create_summarized_experiment(gene_names_list[[k]], design, protein_data_list[[k]])
})

library(SingleCellExperiment)
library(scran)

# Function to create a SummarizedExperiment object and perform initial processing
create_and_process_se <- function(se_list, index, assay_name = "counts") {
  # Create a SingleCellExperiment object
  se <- SingleCellExperiment(assays = list(counts = se_list[[index]][[assay_name]]))
  
  # Replace NA and zero values with a small number
  counts <- assay(se, assay_name)
  counts[is.na(counts) | counts == 0] <- 1e-300
  
  # Update the assay with modified counts
  assays(se)[[assay_name]] <- counts
  
  # Add logcounts
  logcounts(se) <- log2(counts / librarySizeFactors(counts) + 1)
  
  # Correct NaN values in logcounts if any
  logcounts(se)[is.nan(logcounts(se))] <- max(logcounts(se)[!is.nan(logcounts(se))], na.rm = TRUE)
  
  # Set assay names and sample column data
  assayNames(se) <- c(assay_name, "logcounts")
 
  return(se)
}

# Function to find common genes and combine two SingleCellExperiment objects
combine_sces <- function(se1, se2) {
  common_genes <- intersect(rownames(se1), rownames(se2))
  
  # Subset both SCE objects to have the same genes
  se1_common <- se1[common_genes, ]
  se2_common <- se2[common_genes, ]
  
  # Combine both SCEs
  combined_se <- cbind(se1_common, se2_common)
  return(combined_se)
}

# Assuming summarized_experiments is a list of SingleCellExperiment objects
sce_A <- create_and_process_se(summarized_experiments, 3)
sce_B <- create_and_process_se(summarized_experiments, 4)

colData(sce_A)$sample <- "202301"
colData(sce_B)$sample <- "202303" 


# Find and model the gene variance for both experiments
gene_var_rep1 <- modelGeneVar(sce_A)
gene_var_rep2 <- modelGeneVar(sce_B)

# Combine SCEs and identify highly variable genes
sce_combined <- combine_sces(sce_A, sce_B)
gene_var_combined <- combineVar(gene_var_rep1, gene_var_rep2)

# Identify highly variable genes
hvgs <- !is.na(gene_var_combined$bio) & gene_var_combined$bio > 0
sum(hvgs)  # Count the number of HVGs

# Perform PCA on the combined SCE using the HVGs
sce_combined <- runPCA(sce_combined, subset_row = hvgs)

library(reticulate)

# Define a function to get HVGs per dataset
get_hvgs_per_dataset <- function(sce) {
  var_out <- modelGeneVar(sce, method = "loess")
  hvg_out <- var_out[var_out$FDR <= 0.1 & var_out$bio > 0, ]
  hvg_out <- hvg_out[order(hvg_out$bio, decreasing = TRUE), ]
  rownames(hvg_out)
}

# Split the combined SCE object into a list of SCE objects based on the sample
sce_list <- split(sce_combined, sce_combined$sample)

# Get highly variable genes for each dataset
hvgs_per_dataset <- lapply(sce_list, get_hvgs_per_dataset)
names(hvgs_per_dataset) <- names(sce_list)

# Prepare data for Scanorama integration
scelist <- lapply(sce_list, function(sce) t(logcounts(sce)[hvgs_per_dataset[[sce$sample]], ]))
genelist <- hvgs_per_dataset

# Set up Python environment for Scanorama
use_python("/usr/bin/python3")
scanorama <- import("scanorama")

# Integrate datasets using Scanorama
integrated_data <- scanorama$integrate(datasets_full = scelist, genes_list = genelist)

# Prepare integrated data for downstream analysis
intdimred <- do.call(rbind, integrated_data[[1]])
colnames(intdimred) <- paste0("PC_", 1:ncol(intdimred))
rownames(intdimred) <- rownames(sce_combined)

# Add standard deviations for explained variance
stdevs <- apply(intdimred, 2, sd)
attr(intdimred, "varExplained") <- stdevs

# Add Scanorama PCA results to SCE object
reducedDim(sce_combined, "Scanorama_PCA") <- intdimred

# Run UMAP using Scanorama PCs
sce_combined <- runUMAP(sce_combined, dimred = "Scanorama_PCA", n_dimred = 50, ncomponents = 2, name = "UMAP_on_Scanorama")
                  
# Visualization of UMAP plots
Figure3a <- plotReducedDim(sce_combined, dimred = "UMAP_on_Scanorama", colour_by = "celltype") +
  scale_color_aaas() +
  theme(legend.position = c(0.8, 0.9)) +
  Seurat::NoAxes()

Figure3b <- plotReducedDim(sce_combined, dimred = "UMAP_on_Scanorama", colour_by = "clusters") +
  scale_color_npg() +
  theme(legend.position = c(0.8, 0.9)) +
  Seurat::NoAxes()
                  
library(bluster)
library(scater)
library(ggsci)
library(Seurat)

# Clustering the cells using the UMAP dimensionality reduction
sce_combined$cluster_uncorrected <- clusterCells(sce_combined, use.dimred = "UMAP_on_Scanorama", BLUSPARAM = NNGraphParam(k = 100))

# Calculating the modularity ratio
clust_out <- clusterRows(assays(sce_combined)$counts, BLUSPARAM = NNGraphParam(k = 100), full = TRUE)
g <- clust_out$objects$graph
ratio <- pairwiseModularity(g, sce_combined$louvain_SNNk10, as.ratio = TRUE)
rownames(ratio) <- paste0("Cluster ", rownames(ratio))
colnames(ratio) <- paste0("Cluster ", colnames(ratio))

# Creating a proportional table for different clusters
cabbage_exp <- do.call(rbind, lapply(1:5, function(i) {
  clus_name <- paste0("Cluster ", i)
  data.frame(clus = clus_name, prop.table(table(sce_combined[sce_combined$cluster_uncorrected == i,]$celltype)))
}))

Figure3d <- ggplot(cabbage_exp, aes(x = clus, y = Freq, fill = Var1)) +
  geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_aaas() +
  xlab("") + ylab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust = 1, colour = "black"),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

# Finding marker genes
markers_genes <- scran::findMarkers(sce_combined, groups = sce_combined$cluster_uncorrected, lfc = 0.25, pval.type = "all", direction = "up")


# Dot plot for selected genes
genes_to_plot <- c("C1QBP", "CTSB", "SLC39A7", "KRT2", "ELOVL5", "RTN3", "PDCD5", "RBBP7", "SAFB", "CNPY2", "FXR1", "ITGB1", "CAP1", "S100A6", "VIM")
Figure3e <- scater::plotDots(sce_combined, features = genes_to_plot, group = "clusters", center = TRUE) +
  xlab("Cluster") + ylab("") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust = 1, colour = "black"),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
