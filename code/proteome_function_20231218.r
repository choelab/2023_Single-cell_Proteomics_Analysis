# Define the required packages
load_lib <- c("data.table", "dplyr", "DEP", "ggplot2", "ggrepel", "ggpubr", 
              "scater", "scran", "SingleCellExperiment", "rafalib", "remotes",
              "umap", "pbapply", "reshape2", "stringr", "ggvenn", "stringi",
              "BiocManager", "RColorBrewer", "xlsx", "magrittr", "readxl",
              "cowplot", "Matrix", "viridis", "tibble", "forcats", "parallel")

# Install any missing packages
new_packages <- load_lib[!load_lib %in% installed.packages()[,"Package"]]
if(length(new_packages)) install.packages(new_packages)

# Load the required packages
lapply(load_lib, require, character.only = TRUE)

# Set SSL options for HTTP requests (use with caution; disabling SSL verification can be insecure)
httr::set_config(config(ssl_verifypeer = 0L))

# Function to retrieve data from the UniProt database
uniprot_mapping <- function(ids) {
  uri <- 'https://rest.uniprot.org/uniprotkb/search?size=1&query='
  format <- '&fields=accession%2Cgene_names&format=tsv'
  fullUri <- paste0(uri, ids, format)
  dat <- httr::GET(fullUri, httr::accept_json(), httr::add_headers(Accept = 'application/json'))
  return(dat)
}

# Correlation plot function
corrplot2 <- function(data,
                      method = "pearson",
                      sig.level = 0.05,
                      order = "original",
                      diag = FALSE,
                      type = "upper",
                      tl.srt = 90,
                      number.font = 1,
                      number.cex = 1,
                      mar = c(0, 0, 0, 0)) {
  if (!requireNamespace("corrplot", quietly = TRUE)) install.packages("corrplot")
  library(corrplot)
  
  data <- data[complete.cases(data), ]
  mat <- cor(data, method = method)
  p.mat <- cor.mtest(mat, method = method)
  
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  corrplot::corrplot(mat, method = "color", col = col(200), number.font = number.font,
                     mar = mar, number.cex = number.cex, type = type, order = order,
                     addCoef.col = "black", tl.col = "black", tl.srt = tl.srt,
                     p.mat = p.mat, sig.level = sig.level, insig = "blank", diag = diag)
}

# Helper function for correlation testing
cor.mtest <- function(mat, method){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], method = method)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}
