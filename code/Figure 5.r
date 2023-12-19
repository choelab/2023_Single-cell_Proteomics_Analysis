###############################
# title : Figure 5
# author : Jaemyung, Jang (piloter2@kbri.re.kr)
# kenel : R 4.3.0
# Date : Dec. 18, 2023
###############################

library(Seurat)
library(scCustomize)
library(harmony)
library(reticulate)
library(ggplot2)

# Function to perform initial processing on the Seurat object
processA <- function(run) {
  # Define mitochondrial, ribosomal, and hemoglobin genes
  MT_symb <- c('MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6','MT-CYB','MT-CO1','MT-CO2','MT-CO3','MT-ATP6','MT-ATP8','MT-RNR2')
  HB <- rownames(run)[grep("^HB", rownames(run))]
  HB <- HB[-na.omit(match(c("HBS1L", "HBEGF", "HBP2"), HB))]

  # Calculate percentages of mitochondrial, ribosomal, and hemoglobin genes
  run[["percent.mt"]] <- PercentageFeatureSet(run, pattern = MT_symb)
  run[["percent.rb"]] <- PercentageFeatureSet(run, pattern = "^RP[LS]")
  run[["percent.mrp"]] <- PercentageFeatureSet(run, pattern = "^MRP[LS]")
  run[["percent.hb"]] <- PercentageFeatureSet(run, features = HB)

  # Log transform of genes per UMI
  run$log10GenesPerUMI <- log10(run$nFeature_RNA) / log10(run$nCount_RNA)

  # Determine thresholds for filtering
  median_nCount <- median(run$nCount_RNA)
  mad_nCount <- mad(run$nCount_RNA)
  median_nFeature <- median(run$nFeature_RNA)
  mad_nFeature <- mad(run$nFeature_RNA)
  median_percent_MT <- median(run$percent.mt)
  mad_percent_MT <- mad(run$percent.mt)
  
  thresholds_nCount <- c(0, median_nCount + 5 * mad_nCount)
  thresholds_nFeature <- c(0, median_nFeature + 5 * mad_nFeature)
  thresholds_percent_MT <- c(0, median_percent_MT + 2 * mad_percent_MT)

  # Subset the data based on thresholds
  run <- subset(run, subset = (nFeature_RNA >= 500) & 
                         (percent.mt <= thresholds_percent_MT[2]) & 
                         (log10GenesPerUMI > 0.80) & 
                         (nCount_RNA <= thresholds_nCount[2]) & 
                         (nFeature_RNA <= thresholds_nFeature[2]))
  return(run)
}

# Reading and preprocessing data
refObj_GSE168669 <- list()
refObj_GSE168669.obj <- list()
data_dir <- file.path("choelab","SCP","GSE168669") #You can download the single-cell RNA-seq analysis files for the study "GSE168669" from the NCBI Gene Expression Omnibus (GEO) database. 

# List of sample names to read
sample_list <- c("GSM5155455_LNCaP-DMSO_", "GSM5155456_LNCaP-ENZ48_", "GSM5155457_LNCaP-RESA_", 
                 "GSM5155458_LNCaP-RESB_", "GSM5161290_D37_1_", "GSM5161291_D37_2_")

# Read and process each sample
for (i in 1:length(sample_list)) {
  refObj_GSE168669[[i]] <- scCustomize::Read10X_GEO(data_dir = data_dir, sample_list = sample_list[i])
  refObj_GSE168669.obj[[i]] <- CreateSeuratObject(counts = refObj_GSE168669[[i]], min.cells = 10)
  refObj_GSE168669.obj[[i]] <- processA(refObj_GSE168669.obj[[i]])
}

# Adding metadata
treatments <- c('DMSO', 'ENZ48', 'RESA', 'RESB', 'DMSO', 'ENZ48')
genetypes <- c('LNCaP', 'LNCaP', 'LNCaP', 'LNCaP', 'VCaP', 'VCaP')

for (i in seq_along(refObj_GSE168669.obj)) {
  refObj_GSE168669.obj[[i]]$treatment <- treatments[i]
  refObj_GSE168669.obj[[i]]$genetype <- genetypes[i]
  refObj_GSE168669.obj[[i]]$platform <- "Drop-seq"
  refObj_GSE168669.obj[[i]]$orig.ident <- paste0(treatments[i], "-", genetypes[i])
}

# Merge and further process
dat.ref <- Reduce(function(x, y) merge(x, y), refObj_GSE168669.obj)

# Normalization and integration using SCTransform and Harmony
dat.combined <- dat.ref %>%
  SCTransform(assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mt', 'nFeature_RNA', 'S.Score', 'G2M.Score'), verbose = T) %>%
  RunPCA(verbose = FALSE, features = SelectIntegrationFeatures(object.list = refObj_GSE168669.obj, nfeatures = 5000), seed.use = 2000) %>%
  RunHarmony(c("genetype"), assay.use = "SCT", dims = 1:50, max.iter.harmony = 100, max.iter.cluster = 200) %>%
  RunUMAP(reduction = "harmony", dims = 1:50, umap.method = "umap-learn") %>%
  FindNeighbors(reduction = "harmony", dims = 1:50) %>%
  FindClusters(resolution = seq(0.2, 2.0, 0.1), algorithm = 2, verbose = T)

Idents(dat.combined) <- "SCT_snn_res.0.7"


library(reshape2)
library(UCell)

# Assigning cell populations based on cluster identities
assign_population <- function(dat.combined) {
    populations <- c("2a", "2d", "2e", "4b", "3a", "5b", "5c", "5a", "1c", "2b", "1a", "2c", "4a", "1b", "3c", "3b", "extra", "3d")
    for (i in 0:17) {
        cluster_id <- as.character(i)
        population_name <- populations[i + 1]
        dat.combined@meta.data$Population[which(Idents(dat.combined) == cluster_id)] <- population_name
    }
    return(dat.combined)
}

# Update dat.combined with population information
dat.combined <- assign_population(dat.combined)



# Figure5a. Visualization

Figure5a_left <- DimPlot(dat.combined, group.by = "Population", label = TRUE, label.size = 5, pt.size = 0.01, raster = FALSE) +
  theme(legend.position = "none", axis.title = element_blank(), plot.title = element_text(size = 12, face = "bold", hjust = 0), 
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  ggtitle("")

# Supplementary Figure. Dot plot for specified genes
gene_list <- c("CD9", "CD63", "CD81", "TSPAN1", "ROM1", "CNPY2", "DNAJA3", "EIF5", "BTF3", "GTPBP4", "GPI", "ATP5ME", "FXR1", "ITGB1", "CAP1", "GTF2E1", "SLC39A7", "DDX39A", "NIPSNAP1", "CTSB", "ATP6V0C", "EIF3C", "PGAM1", "C1QBP", "HSPA5", "ENO1", "GAPDH", "PPIA", "S100A6", "VIM", "EEF2", "EIF5A", "PKM", "SAFB", "RBBP7", "PDCD5", "RTN3", "RO60", "ELOVL5", "RAB14", "KRT2", "ARPC4", "YWHAQ")

DotPlot(dat.combined, idents = names(table(Idents(dat.combined))[which(table(Idents(dat.combined)) > 100)],
        group.by = "Population", features = gene_list, cols = "RdBu", assay = "SCT") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Read marker file and prepare for UCell scoring
SCP.marker <- readRDS("/choelab/SCP/proteomeDatafromPD/SCP_markers_20231011.RDS")

# Extract markers for each cluster
marker <- list(
    cluster1 = SCP.marker %>% filter(FDR < 0.1, summary.logFC > 1, cluster == '1'),
    cluster2 = SCP.marker %>% filter(FDR < 0.1, summary.logFC > 1.5, cluster == '2'),
    cluster3 = SCP.marker %>% filter(FDR < 0.1, summary.logFC > 1.5, cluster == '3'),
    cluster4 = SCP.marker %>% filter(summary.logFC > 0.1, cluster == '4'),  # Custom conditions
    cluster5 = SCP.marker %>% filter(FDR < 0.1, summary.logFC > 1.5, cluster == '5')
)

# Add UCell scores
dat.combined.score <- AddModuleScore_UCell(dat.combined, features = marker, ncores = 4, assay = "RNA", w_neg = 1)
signature.names <- paste0(names(marker), "_UCell")

# Prepare data for plotting
df <- melt(dat.combined.score@meta.data %>% select("Population", signature.names), id.vars = "Population", na.rm = TRUE)

# Plotting module scores for each cluster
plot_module_scores <- function(df, variable_name, threshold) {
    ggplot(df %>% filter(variable == variable_name), aes(x = Population, y = value)) + 
      geom_hline(yintercept = threshold, color = "red") +
      geom_boxplot(notch = TRUE) + theme_bw()
}

# Supplementary Figure. Adjust variable names and thresholds as necessary
plot_module_scores(df, names(table(df$variable))[1], 0.082)
plot_module_scores(df, names(table(df$variable))[2], 0.075)

library(monocle3)
library(SeuratWrappers)
library(igraph)

# Function to identify the correct root state for trajectory analysis
get_correct_root_state <- function(cds, cell_phenotype, root_type) {
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex_matrix <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex_matrix[cell_ids,]))))]
  return(root_pr_nodes)
}

# Convert Seurat object to CellDataSet for Monocle
cds <- as.cell_data_set(dat.combined)

# Clustering cells in the dataset
cds <- cluster_cells(cds, reduction_method = "UMAP")

# Learning the trajectory graph
cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE, learn_graph_control = list(minimal_branch_len = 75))

# Assigning cluster identities from Seurat to Monocle CDS
cds@clusters@listData$UMAP[["clusters"]] <- dat.combined$Population
cds@colData$cell_type <- dat.combined$Population

# Determining the root state for trajectory analysis
DN_node_id <- get_correct_root_state(cds, cell_phenotype = 'cell_type', "5a")

# Ordering cells along the trajectory
cds_ordered <- order_cells(cds, root_pr_nodes = DN_node_id)

#  Plotting trajectory graphs
plot_trajectory_graph <- function(cds_ordered, color_by, file_name) {
  pdf(file = file_name, width = 7.5, height = 7.5)
  plot_cells(cds_ordered,
             color_cells_by = color_by,
             trajectory_graph_color = ifelse(color_by == "pseudotime", "grey50", "black"),
             trajectory_graph_segment_size = ifelse(color_by == "pseudotime", 0.25, 0.5),
             label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = FALSE,
             graph_label_size = 4, cell_size = 0.5) +
    xlim(c(floor(min(cds@int_colData@listData$reducedDims@listData$UMAP[,1])), ceiling(max(cds@int_colData@listData$reducedDims@listData$UMAP[,1])))) +
    ylim(c(floor(min(cds@int_colData@listData$reducedDims@listData$UMAP[,2])), ceiling(max(cds@int_colData@listData$reducedDims@listData$UMAP[,2])))) +
    Seurat::NoLegend() + Seurat::NoAxes()
  dev.off()
}

# Figure 5a. Plotting cells colored by cell type
plot_trajectory_graph(cds_ordered, "cell_type", "/choelab/SCP/results/Figure 5a_Right.pdf")

# Plotting cells colored by pseudotime
# plot_trajectory_graph(cds_ordered, "pseudotime", "/choelab/SCP/results/GSE168669_pseudotime.pdf")

library(magrittr)
library(igraph)
library(gam)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Calculate closest vertices for cells in the trajectory graph
y_to_cells <- principal_graph_aux(cds1)$UMAP$pr_graph_cell_proj_closest_vertex %>%
  as.data.frame()
y_to_cells$cells <- rownames(y_to_cells)
y_to_cells$Y <- y_to_cells$V1

mst <- principal_graph(cds1)$UMAP
root <- cds1@principal_graph_aux$UMAP$root_pr_nodes
endpoints <- names(which(igraph::degree(mst) == 1))
endpoints <- endpoints[!endpoints %in% root]

# Calculate cell weights for each endpoint
cellWeights <- lapply(endpoints, function(endpoint) {
  path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
  path <- as.character(path)
  df <- y_to_cells[y_to_cells$Y %in% path, ]
  df <- data.frame(weights = as.numeric(colnames(cds1) %in% df$cells))
  colnames(df) <- endpoint
  return(df)
}) %>% do.call(what = 'cbind', args = .) %>%
  as.matrix()

rownames(cellWeights) <- colnames(cds1)
pseudotime <- matrix(pseudotime(cds1), ncol = ncol(cellWeights), nrow = nrow(cds1))

# Segmenting graph for each endpoint
cds_sub <- lapply(endpoints, function(ep) {
  choose_graph_segments(cds1, reduction_method = "UMAP", starting_pr_node = root, ending_pr_nodes = ep, return_list = TRUE)
})
cds_subset <- lapply(cds_sub, function(cds_sb) cds1[, rownames(colData(cds1)) %in% cds_sb$cells])

# Analyzing gene expression for each endpoint
topgene <- lapply(seq_along(endpoints), function(nt) {
  t <- pseudotime[, nt]
  Y <- log1p(GetAssayData(dat.combined))
  var100 <- names(apply(Y, 1, var)[apply(Y, 1, var) > 0.1])
  Y <- Y[var100, ]

  gam.pval <- apply(Y, 1, function(z) {
    d <- data.frame(z = z, t = t)
    suppressWarnings(tmp <- gam(z ~ lo(t), data = d))
    p <- summary(tmp)[3][[1]][2,3]
    p
  })

  res <- gam.pval
  return(res)
})

topgenes <- lapply(topgene, function(nt) {
  res <- names(sort(nt[nt < 0.00001], decreasing = FALSE))
  return(res)
})

# Creating heatmaps for gene expression
htkm <- list()
for(x in seq_along(endpoints)) {
  modulated_genes <- graph_test(cds_subset[[x]], neighbor_graph = "principal_graph", cores = 8)
  pt.matrix <- exprs(cds_subset[[x]])[match(intersect(topgenes[[x]], rownames(modulated_genes)), rownames(rowData(cds_subset[[x]]))), order(pseudotime(cds_subset[[x]]))]
  pt.matrix <- t(apply(pt.matrix, 1, function(x) smooth.spline(x, df = 3)$y))
  pt.matrix <- t(apply(pt.matrix, 1, function(x) (x - mean(x)) / sd(x)))
  
  labels <- rep("", nrow(pt.matrix))
  labels[rownames(pt.matrix) %in% insGenes] <- rownames(pt.matrix)[rownames(pt.matrix) %in% insGenes]

  htkm[[x]] <- Heatmap(pt.matrix, name = "z-score", col = colorRamp2(seq(from = -2, to = 2, length = 11), rev(brewer.pal(11, "Spectral"))), show_row_names = TRUE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 6), km = 4, row_title_rot = 0, cluster_rows = TRUE, cluster_row_slices = TRUE, cluster_columns = FALSE) + rowAnnotation(labels = anno_text(labels, which = "row"), width = max(grobWidth(textGrob(labels))))
}

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggVennDiagram)
library(enrichR)

HM <- draw(htkm[[3]])
# List of common genes for scProteomics
scProteomicsGenes <- c( "PSMD12",    "CLIC1",     "PHGDH",     "ACTN4",     "SYNCRIP",   "GGCT",      "CS",        "SNRNP200",  "IDH1",      "GLUD1",
                  "LMNA",      "ALDOA",     "KRT1",      "GAPDH",     "HSPB1",     "KRT8",      "S100A6",    "ENO1",      "NPM1",      "TPM3.1",  
                  "EPHX1",     "DBI",       "LDHB",      "P4HB",      "PFN1",      "HSP90AA1",  "HNRNPC",    "ANXA6",     "VIM"  ,     "ANXA5",  
                  "RPSA",      "HMGB1",     "HNRNPA1",   "PARP1",     "HSPA1B",    "H1-4",      "HSPA5",     "PABPC1",    "XRCC6",     "EEF2",     
                  "KRT10",     "PDIA4",     "ETFA",      "PRKCSH",    "PKM",       "HNRNPL",    "EZR",       "RPS2",      "H2AX" ,     "H1-5" ,    
                  "ATP2A2",    "PGAM1",     "NCL",       "COX5A",     "LMNB1",     "MAOA",      "PAICS",     "UBA1",      "NME2" ,     "HNRNPA2B1",
                  "UQCRC2",    "SFPQ",      "PPIB",      "RPS3",      "CFL1",      "EIF4B",     "ATP5PB",    "ACAT1",     "ATP5F1A",   "RPL13"   , 
                  "PTBP1",     "CALR",      "CANX",      "PRDX3",     "PDIA3",     "PPP2R1A",   "YWHAB",     "STIP1",     "SHMT2" ,    "HSPA4"   , 
                  "PHB1",      "KRT9",      "ETFB",      "HSPA9",     "RPL3",      "ANP32A",    "MDH2",      "HADHA" ,    "GARS1" ,    "USP5",
                  "RPS9",      "CCT5",      "CCT3",      "TUFM" ,     "AARS1",     "ACADVL",    "CCT8",      "ACLY" ,     "YARS1" ,    "USP14",  
                  "RAD23B",    "CSE1L",     "VCP",       "EIF4A1",    "ACTR2",     "RPS3A",     "HSPE1",     "SEC61A1",   "RPS7"  ,    "RPS13",
                  "RPL7A",     "RAB11A",    "RPS4X",     "RPL23A",    "H4C1",      "RAB1A",     "RAN" ,      "RPS15",     "RPL11" ,    "EIF5A",
                  "TPM4",      "H3C1",      "CCT2",      "PRKDC",     "CLTC",      "HNRNPU",    "SET.1",     "OGDH",      "FKBP4" ,    "RPL6", 
                  "EEF1A2",    "PRDX1",     "C1QBP",     "CKAP4",     "DHX9",      "AHNAK" ,    "ILF3",      "TRIM28",    "DYNC1H1",   "GANAB.1",
                  "PLEC",      "NONO",      "HADH",      "HSDL2",     "DCXR",      "ALYREF",    "CAND1",     "DDX17",     "KHSRP" ,    "PARK7", 
                  "ACO2",      "CCT7",      "API5",      "PNPO",      "RRBP1",     "MYO6",      "PA2G4",     "PGK1" ,     "H2BC11",    "STT3A"  ,
                  "DDX39B",    "HNRNPM",    "SUB1",      "RACK1",     "PRDX2",     "RPS27A",    "ATP1A1",    "OXCT1",     "RPN1"  ,    "ANP32B" , 
                  "SF3A3",     "TRAP1",     "XRCC5",     "CTSD",      "PRDX6",     "GDI2",      "SND1",      "CCT6A",     "EPRS1" ,    "FH",  
                  "ACTC1",     "SLC25A6",   "EEF1G",     "TCP1",      "LRPPRC",    "TUBB",      "HSPH1",     "ANXA2",     "ECH1"  ,    "PPIA",
                  "VARS1",     "SSB",       "PPA1" ,     "MYL12B",    "EIF4G1",    "YWHAG",     "RPL18A",    "PRMT1",     "MYO1C" ,    "RPL15",  
                  "PHB2",      "PGD",       "YWHAZ",     "SLC25A11",  "NQO1",      "SLC25A5",   "TALDO1" ,   "PSMD3",     "ESD"   ,    "GCN1" ,  
                  "FBL",       "MAP4",      "ATP5MG",    "HYOU1",     "TARS1",     "TUBB4B",    "SLC25A3",   "COPG1",     "HK1"    ,   "PPA2" ,
                  "LDHA",      "MACROH2A1", "MTHFD1" ,   "H2AZ1" ,    "CAP1",      "CKB",       "PUF60")

# Extract scRNAseq genes from a heatmap object
scRNAseqGenes <- HM@ht_list[["z-score"]]@row_names_param[["labels"]]

commonGenes <- list(  scRNAseq = scRNAseqGenes,
                  scProteomics = scProteomicsGenes)
                       
# Generate heatmap and extract cluster information
drawHeatmapAndGetClusters <- function(htkm) {
  HM <- draw(htkm)
  row_order(HM) 
}

# Create a dataframe of genes and their corresponding clusters
generateClusterDataFrame <- function(clusterlist, pt_matrix) {
  clu_df <- lapply(seq_along(clusterlist), function(i) {
    data.frame(GeneID = rownames(pt_matrix)[clusterlist[[i]]],
               Cluster = paste0("cluster", i),
               stringsAsFactors = FALSE)
  }) %>% 
  do.call(rbind, .)
  return(clu_df)
}

# Combine into a commonGenes list
commonGenes <- list(
  scRNAseq = scRNAseqGenes,
  scProteomics = scProteomicsGenes
)


# Venn Diagram Analysis
performVennDiagramAnalysis <- function(commonGenes) {
  ggVennDiagram(commonGenes) + scale_fill_gradient(low = "blue", high = "red")
}

# Perform Enrichment Analysis
performEnrichmentAnalysis <- function(genes, databases) {
  enrichr(genes = genes, databases = databases)
}

# Generate Enrichment Comparison Plot
generateEnrichmentComparisonPlot <- function(enrich_results_S, enrich_results_P, dbs) {
  enrich_comparsion <- c()
  for (i in seq_along(dbs)) {
    enrich_comparsion <- rbind(enrich_comparsion, 
      data.frame("num" = setdiff(enrich_results_S[[i]] %>% dplyr::filter(P.value < 0.05) %>% pull(Term),
                                 enrich_results_P[[i]] %>% dplyr::filter(P.value < 0.05) %>% pull(Term)) %>% length(), 
                 "type" = "scRNA-seq", 
                 "group" = dbs[i]),
      data.frame("num" = intersect(enrich_results_S[[i]] %>% dplyr::filter(P.value < 0.05) %>% pull(Term),
                                   enrich_results_P[[i]] %>% dplyr::filter(P.value < 0.05) %>% pull(Term)) %>% length(), 
                 "type" = "common", 
                 "group" = dbs[i]),
      data.frame("num" = setdiff(enrich_results_P[[i]] %>% dplyr::filter(P.value < 0.05) %>% pull(Term),
                                 enrich_results_S[[i]] %>% dplyr::filter(P.value < 0.05) %>% pull(Term)) %>% length(), 
                 "type" = "scProteomics", 
                 "group" = dbs[i]))
  }

  ggplot(enrich_comparsion, aes(group, num, fill = type)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    theme_bw() +
    xlab("") +
    ylab("Number of pathways") +
    scale_fill_aaas() +
    theme(legend.position = c(0.8, 0.8), 
          axis.text.x = element_text(size = 16, colour = "black"),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.title = element_blank(),
          legend.background = element_blank()) +
    coord_flip()
}

# Figure 5c. Extract row order from heatmap and generate cluster dataframe
clusterlist <- drawHeatmapAndGetClusters(htkm[[3]])
clu_df <- generateClusterDataFrame(clusterlist, pt.matrix)

# Figure5b. Venn Diagram of common genes between scRNAseq and scProteomics
Figure5b <- performVennDiagramAnalysis(commonGenes)

# Enrichment analysis
databases <- enrichR::listEnrichrDbs()$libraryName
dbs <- c(grep("GO", databases, value = TRUE)[c(20:22)], grep("WikiPathway_2023_Human|Reactome_2022|KEGG_2021_Human", databases, value = TRUE))
enrich_results_S <- performEnrichmentAnalysis(commonGenes[[1]], dbs)
enrich_results_P <- performEnrichmentAnalysis(commonGenes[[2]], dbs)

# Figure5e. Generate comparison plot
Figure5e <- generateEnrichmentComparisonPlot(enrich_results_S, enrich_results_P, dbs)

df_results <- data.frame()

processEnrichResults <- function(enrich_data, dbs, type, filter_function) {
  for(i in seq_along(dbs)) {
    enrich_result <- enrich_data[[i]][enrich_data[[i]]$Term %in% filter_function(enrich_data[[i]]), ]
    selected_terms <- enrich_result %>%
      filter(grepl("Androgen|Prostate|KLK|PKN|RHO GTPase|Protein Kinase|WASPs|WAVEs", Term) & 
             !grepl("Negative|Positive", Term))
    df_results <<- rbind(df_results, data.frame(selected_terms, group = dbs[i], type = type))
  }
}

# Define filter functions for different comparisons
filter_scProteomics <- function(data) setdiff(data %>% filter(P.value < 0.05) %>% pull(Term), enrich_results.S[[i]] %>% filter(P.value < 0.05) %>% pull(Term))
filter_scRNAseq <- function(data) intersect(data %>% filter(P.value < 0.05) %>% pull(Term), enrich_results.P[[i]] %>% filter(P.value < 0.05) %>% pull(Term))

# Apply the function to process enrichment results
processEnrichResults(enrich_results.P, dbs, "scProteomics", filter_scProteomics)
processEnrichResults(enrich_results.S, dbs, "scRNA-seq", filter_scRNAseq)

# Figure 5f. Plotting
Figure5f <- ggplot(df_results, aes(Term, -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", aes(fill = group)) +
  scale_fill_aaas() +
  xlab("") +
  ylab("-Log 10 Adjust P-value") +
  theme_bw() +
  theme(legend.position=c(0.25,0.9),
        axis.text.x = element_text(size = 16, hjust = 0, colour = "black"),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.background = element_blank()) +
  coord_flip() +
  facet_wrap(~type)

# Display the plot
print(pPATHWAYs)
