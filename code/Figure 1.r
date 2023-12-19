###############################
# title : Figure 1
# author : Jaemyung, Jang (piloter2@kbri.re.kr)
# kenel : R 4.3.0
# Date : Dec. 18, 2023
###############################

# Define the path to the data directory
path <- file.path("/", "choelab", "SCP")

source(paste0(path, "/code/proteomic_function_20231218.r"))

# Define file patterns for raw data
rawPD_pattern <- "_Proteins.txt$"
# rawPeptide_pattern <- "_PeptideGroups.txt$"

# Get the list of files matching the pattern
rawPD <- list.files(list.dirs(paste0(path,"/proteomicDatafromPD")), pattern = rawPD_pattern, full.names = TRUE)
# rawPeptide <- list.files(list.dirs(paste0(path,"/proteomicDatafromPD")), pattern = rawPeptide_pattern, full.names = TRUE)

# Function to clean column names
clean_colnames <- function(df) {
  colnames(df) <- str_replace_all(colnames(df), " ", ".")
  colnames(df) <- str_replace_all(colnames(df), "-", "_")
  return(df)
}

# Read and clean data from files
read_and_clean <- function(filenames) {
  lapply(filenames, function(file) {
    clean_colnames(fread(file))
  })
}

# Read protein and peptide data
datafromPD <- read_and_clean(rawPD)
#peptidefromPD <- read_and_clean(rawPeptide)

# Process the counts data
counts <- lapply(datafromPD, function(data) {
  data %>%
    select(Accession, starts_with("Abundances."))
})

# Process the peptide counts data
# Peptidecounts <- lapply(peptidefromPD, function(data) {
#  data %>%
#    select(Annotated.Sequence, starts_with("Found.in."))
# })

# Define a function to extract a subset of data based on column names
subset_data <- function(data, pattern) {
  select_cols <- colnames(data)[grep(pattern, colnames(data), value = TRUE)]
  return(subset(data, select = c("Accession", select_cols)))
}

# Define a function to extract unique non-missing values from a data frame
get_unique_non_missing <- function(data) {
  result <- unique(c(unlist(apply(data, 2, function(x) which(!is.na(x)))[-1])))
  return(data$Accession[result])
}

# Define a function to create a Venn diagram
create_venn_diagram <- function(y) {
  require(ggVennDiagram)
  p <- ggVennDiagram(y, label_alpha = 0, 
                     color = c("A" = "yellow", "B" = "steelblue", 'C' = "red"),
                     set_color = c("A" = "yellow", "B" = "steelblue", 'C' = "red")) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    ggtitle("Venn Diagram of none vs Abeta vs AST011AandAbeta")
  return(p)
}

# Define a function to fetch data from UniProt
fetch_data_from_uniprot <- function(ids) {
  data <- uniprot_mapping(ids)
  content(data, as = "text", encoding = 'UTF-8')
  res <- unlist(str_split(unlist(str_split(unlist(str_split(content(data, as = "text", encoding = 'UTF-8'), "\\t")), "\\n"))[4], " "))[1]
  return(res)
}

# Define a function to extract gene names from descriptions
extract_gene_names <- function(desc) {
  str_sub(desc, str_locate(desc, pattern = "GN=")[2] + 1, str_locate(desc, pattern = "PE=")[1] - 2)
}

# Part 1: Modify variables
k = 1
y.count <- list(
  "DDM" = subset_data(counts[[k]], "Sample.DDM$"),
  "SB3-12" = subset_data(counts[[k]], "1_40772$"),
  "Q" = subset_data(counts[[k]], "Point3_Q$")
)

y <- list(
  "DDM" = get_unique_non_missing(y.count[["DDM"]]),
  "SB3-12" = get_unique_non_missing(y.count[["SB3-12"]]),
  "Q" = get_unique_non_missing(y.count[["Q"]])
)

# Create Venn diagram
figure1h <- create_venn_diagram(y)


# Fetch data from UniProt and create data frames
resultUNIPROT <- pbapply::pblapply(datafromPD[[k]]$Accession, fetch_data_from_uniprot)

results_from_uniprot <- data.frame('Accession' = datafromPD[[k]]$Accession, 'Gene.names' = unlist(resultUNIPROT))

# Part 2: Not imputed
data.merged_neo <- subset(data, select = c("Gene.names", "Accession", grep("^Abundance.Ratio.log2.|^Abundance.Ratio.Weight.|^Abundance.Ratio.P_Value", colnames(data), value = TRUE)))

# Define a function to calculate the number of identified proteins
calculate_protein_length <- function(data_unique, sample_name) {
  data_selected <- subset(data_unique, select = c("Accession", colnames(data_unique)[c(2:79)][which(str_sub(colnames(data_unique)[c(2:79)], str_locate(colnames(data_unique)[c(2:79)], "Sample")[, 2] + 2) %in% names(table(str_sub(colnames(data_unique)[c(2:79)], str_locate(colnames(data_unique)[c(2:79)], "Sample")[, 2] + 2)))[sample_name])]))
  return(apply(data_selected, 2, function(x) length(which(!is.na(x))))[-1] %>% melt())
}

# Calculate protein length for different samples
sample_names <- names(table(str_sub(colnames(data_unique)[c(2:79)], str_locate(colnames(data_unique)[c(2:79)], "Sample")[, 2] + 2)))
protein_length <- data.frame()

for (sample_name in sample_names) {
  protein_length <- rbind(protein_length, cbind(calculate_protein_length(data_unique, sample_name), "type" = sample_name))
}

# Create boxplot and jitter plot
require(ggsci)
require(ggplot2)
require(ggpubr)

protein_length_plot <- protein_length %>%
  ggplot(aes(type, value, color = type)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) +
  scale_color_aaas() +
  xlab("") + ylab("Number of identified proteins") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust = 1, colour = "black"),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 16),
    legend.title = element_blank(),
    legend.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# Rename sample types for better readability
protein_length$type[which(protein_length$type %in% "DDM")] <- "1. 0.2% DDM"
protein_length$type[which(protein_length$type %in% "1_40772")] <- "2. 1% SB3-14"
protein_length$type[which(protein_length$type %in% "Point3_Q")] <- "3. 0.3% Q Catalyst"
protein_length$type[which(protein_length$type %in% "1_40772_Point3_Q")] <- "4. 1% SB3-14 with 0.3% Q Catalyst"
protein_length$type[which(protein_length$type %in% "1_40772_Point125_Q")] <- "5. 1% SB3-14 with 0.125% Q Catalyst"

# Create another boxplot and jitter plot with renamed sample types
figure1i <-  <- protein_length_plot
