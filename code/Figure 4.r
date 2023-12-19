library(monocle3)
library(sceasy)
library(gam)
library(forcats)
library(ggplot2)
library(Seurat)

# Function to get the correct root state for trajectory analysis
get_correct_root_state <- function(cds, cell_phenotype, root_type) {
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))]       
  root_pr_nodes
}

# Prepare data for Monocle
pData <- sce_combined@colData
fData <- data.frame(gene_short_name = rownames(sce_combined), row.names = rownames(sce_combined))
expression_matrix <- assays(sce_combined)$counts
expression_matrix <- as.matrix(expression_matrix)

# Create CDS object
cds <- new_cell_data_set(expression_matrix, pData, fData)
rownames(colData(cds)) <- paste0(sce_combined@colData$sample, "_", rownames(colData(sce_combined)))
cds@int_colData@listData$reducedDims@listData$PCA <- reducedDim(sce_combined, "PCA")
cds@int_colData@listData$reducedDims@listData$UMAP <- reducedDim(sce_combined, "UMAP_on_Scanorama")

# Cluster cells and learn the graph
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData$UMAP[["clusters"]] <- sce_combined$cluster_uncorrected
cds@colData$cell_type <- sce_combined$celltype
cds <- learn_graph(cds, use_partition = FALSE, verbose = FALSE, learn_graph_control = list(minimal_branch_len = 5))

# Determine the root node
DN_node_id <- get_correct_root_state(cds, cell_phenotype = 'cluster_uncorrected', root_type = "5")
cds_ordered <- order_cells(cds, root_pr_nodes = DN_node_id)

# Visualization of cells in pseudotime
pseudotime <- pseudotime(cds_ordered) 
cds_ordered@colData$pseudotime <- pseudotime

# Figure4a. Plot cells
plot_cells(cds_ordered, color_cells_by = "cell_type", trajectory_graph_color = "grey50", trajectory_graph_segment_size = 0.25,
           label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = FALSE, graph_label_size = 4, cell_size = 0.5) +
  xlim(c(floor(min(cds@int_colData@listData$reducedDims@listData$UMAP[,1])), ceiling(max(cds@int_colData@listData$reducedDims@listData$UMAP[,1])))) +
  ylim(c(floor(min(cds@int_colData@listData$reducedDims@listData$UMAP[,2])), ceiling(max(cds@int_colData@listData$reducedDims@listData$UMAP[,2])))) +
  Seurat::NoAxes()

# Extract closest vertices for each cell
y_to_cells <- principal_graph_aux(cds_ordered)$UMAP$pr_graph_cell_proj_closest_vertex %>%
  as.data.frame()
y_to_cells$cells <- rownames(y_to_cells)
y_to_cells$Y <- y_to_cells$V1

# Get the minimum spanning tree (MST)
mst <- principal_graph(cds_ordered)$UMAP

# Extract root and endpoints
root <- cds_ordered@principal_graph_aux$UMAP$root_pr_nodes
endpoints <- names(which(igraph::degree(mst) == 1))
endpoints <- endpoints[!endpoints %in% root]

# Calculate cell weights for each endpoint
cellWeights <- lapply(endpoints, function(endpoint) {
  path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
  path <- as.character(path)
  df <- y_to_cells[y_to_cells$Y %in% path, ]
  df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
  colnames(df) <- endpoint
  return(df)
}) %>% do.call(what = 'cbind', args = .) %>%
  as.matrix()

rownames(cellWeights) <- colnames(cds_ordered)
pseudotime <- matrix(pseudotime(cds_ordered), ncol = ncol(cellWeights), nrow = ncol(cds_ordered), byrow = FALSE)

# Visualization of trajectory
cds_sub <- lapply(endpoints, function(ep) {
  choose_graph_segments(cds_ordered, reduction_method = "UMAP", starting_pr_node = root, ending_pr_nodes = ep, return_list = TRUE)
})

cds_subset <- lapply(cds_sub, function(cds_sb) {
  cds_ordered[, rownames(colData(cds1)) %in% cds_sb$cells]
})

# Figure4a. Plot cells for a specific trajectory segment
plot_cells(cds_subset[[3]],
           color_cells_by = "pseudotime",
           show_trajectory_graph = TRUE,
           trajectory_graph_color = "grey50",
           trajectory_graph_segment_size = 0.1,
           label_cell_groups = TRUE,
           label_roots = TRUE,
           label_leaves = TRUE,
           label_branch_points = FALSE,
           graph_label_size = 4,
           cell_size = 0.5) +
  xlim(c(floor(min(cds@int_colData@listData$reducedDims@listData$UMAP[, 1])), ceiling(max(cds@int_colData@listData$reducedDims@listData$UMAP[, 1])))) +
  ylim(c(floor(min(cds@int_colData@listData$reducedDims@listData$UMAP[, 2])), ceiling(max(cds@int_colData@listData$reducedDims@listData$UMAP[, 2])))) +
  Seurat::NoAxes()

library(magrittr)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Calculate correlations for each gene
calculate_correlations <- function(cds, pseudotime_colname, method) {
  cor_results <- lapply(1:nrow(cds), function(n) {
    cor.test(
      t(logcounts(cds))[, n], 
      as.data.frame(colData(cds))[, pseudotime_colname], 
      method = method
    )
  })
  do.call(rbindlist, cor_results)
}

# Add gene names and grouping based on p-values and correlation estimates
format_correlation_results <- function(cor_results, gene_names, thresholds) {
  cor_results$names <- gene_names
  cor_results$group <- "lightgrey"
  
  cor_results$group[cor_results$names %in% (cor_results %>% filter(log10(p.value) < -5) %>% pull(names))] <- "darkgreen"
  cor_results$group[cor_results$names %in% (cor_results %>% filter(log10(p.value) < -5 & estimate > thresholds$positive) %>% pull(names))] <- "red"
  cor_results$group[cor_results$names %in% (cor_results %>% filter(log10(p.value) < -5 & estimate < thresholds$negative) %>% pull(names))] <- "royalblue"
  
  cor_results
}

# Correlation for pseudotime
cor_res_pseudotime <- calculate_correlations(cds2, "pseudotime", "spearman")
cor_res_pseudotime_formatted <- format_correlation_results(
  cor_res_pseudotime, 
  rownames(cds2), 
  list(positive = 0.15, negative = -0.15)
)

# Correlation for AR gene
cor_res_ar <- calculate_correlations(cds2, "AR", "spearman")
cor_res_ar_formatted <- format_correlation_results(
  cor_res_ar, 
  rownames(cds2), 
  list(positive = 0.2, negative = -0.2)
)

# Plot correlation results
plot_correlation_results <- function(cor_results, x_label, y_label) {
  ggplot(cor_results, aes(x = estimate, y = -log10(p.value))) +
    geom_point(color = cor_results$group) +
    geom_label_repel(aes(label = ifelse(names %in% your_genes_of_interest, names, "")),
                     max.overlaps = Inf, min.segment.length = 0.25, direction = "both") +
    geom_vline(xintercept = c(-0.15, 0.15), linetype = 'dotted', color = 'red', size = 1) +
    geom_hline(yintercept = 5, linetype = 'dotted', color = 'red', size = 1) +
    xlab(x_label) + ylab(y_label) +
    theme_bw()
}

# Figure3b. Plot for pseudotime
plot_correlation_results(cor_res_pseudotime_formatted, "Spearman's correlation coeffient to pseudotime", "-Log10 P-value")

library(ComplexHeatmap)
library(circlize)
library(enrichR)
library(ggsci)

# Generate heatmaps for each endpoint
generate_heatmap <- function(cds_subset, topgenes, endpoint_index) {
  modulated_genes <- graph_test(cds_subset[[endpoint_index]], neighbor_graph = "principal_graph", cores = 8)
  pt.matrix <- exprs(cds_subset[[endpoint_index]])[match(intersect(topgenes[[endpoint_index]], modulated_genes$gene_short_name), rownames(rowData(cds_subset[[endpoint_index]]))), order(pseudotime(cds_subset[[endpoint_index]]))]
  pt.matrix <- t(apply(pt.matrix, 1, function(x) { smooth.spline(x, df = 3)$y }))
  pt.matrix <- t(apply(pt.matrix, 1, function(x) { (x - mean(x)) / sd(x) }))
  
  labels <- rep("", nrow(pt.matrix))
  labels[rownames(pt.matrix) %in% unique_genes] <- rownames(pt.matrix)
  
  Heatmap(
    pt.matrix,
    name = "z-score",
    col = colorRamp2(seq(from = -2, to = 2, length = 11), rev(brewer.pal(11, "Spectral"))),
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 6),
    km = 4,
    row_title_rot = 0,
    cluster_rows = TRUE,
    cluster_row_slices = TRUE,
    cluster_columns = FALSE
  ) + rowAnnotation(labels = anno_text(labels, which = "row"), width = max(grobWidth(textGrob(labels))))
}

# Initialize heatmap list
htkm <- lapply(seq_len(length(endpoints)), function(x) generate_heatmap(cds_subset, topgenes, x))

# Enrichment analysis
perform_enrichment_analysis <- function(genes, databases, filters) {
  enrich_results <- enrichr(genes = genes, databases = databases)
  lapply(enrich_results, function(df) {
    df %>% filter(P.value < 0.1) %>% filter(str_detect(Term, filters))
  })
}

databases <- listEnrichrDbs()$libraryName
enrich_results <- perform_enrichment_analysis(insGenes, c(grep("GO", databases, value = TRUE)[20:22], grep("WikiPathway_2023_Human|Reactome_2022|KEGG_2021_Human", databases, value = TRUE)), "Androgen|Prostate|KLK|PKN|RHO GTPase|Protein Kinase|WASPs|WAVEs")

# Combine and format enrichment results for visualization
format_enrichment_results <- function(enrich_results, enrichment_type) {
  rbind(
    data.frame(enrich_results[[1]], group = "GO-BP", enrichment = enrichment_type),
    data.frame(enrich_results[[5]], group = "Reactome", enrichment = enrichment_type),
    data.frame(enrich_results[[6]], group = "Wikipathway", enrichment = enrichment_type)
  )
}

df_results <- rbind(
  format_enrichment_results(enrich_results, "Positive"),
  format_enrichment_results(enrich_results.n, "Negative")
)

# Figure3c. Visualization of enrichment results
ggplot(df_results, aes(Term, -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", aes(fill = group)) +
  scale_fill_aaas() +
  coord_flip() +
  facet_wrap(~enrichment) +
  xlab("") +
  ylab("-Log10 Adjusted P-value") +
  theme_bw()

library(stringr)
pseudo <- list()

for(geneX in unique(unlist(str_split(df_results$Genes, ";")))) {
  
  # Filter the dataset for the current gene
  cds_subset_gene <- cds_subset[[3]][rowData(cds_subset[[3]])$gene_short_name %in% geneX, ]
  
  # Create a dataframe for plotting
  cds_df <- data.frame(
    pseudotime = pseudotime(cds_subset_gene),
    expression = as.numeric(assays(cds_subset_gene)$counts),
    group = str_split_fixed(colnames(cds_subset_gene), "_", n = 2)[, 1],
    celltype = colData(cds_subset_gene)$cell_type
  )
  
  # Filter out very low expressions
  cds_df_filtered <- cds_df[cds_df$expression > 1e-10, ]
  
  # Figure3d. Generate the plot
  pseudo[[geneX]] <- ggplot(cds_df_filtered, aes(x = pseudotime, y = log10(expression))) +
    geom_point() +
    geom_smooth(method = "glm", formula = y ~ splines::bs(x, 3), se = TRUE, fullrange = TRUE) +
    theme_bw() +
    ggtitle(geneX) +
    theme(
      legend.position = c(0.8, 0.8),
      axis.text.x = element_text(size = 16, colour = "black"),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.title = element_blank(),
      legend.background = element_blank()
    )
}
