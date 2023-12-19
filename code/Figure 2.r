###############################
# title : Figure 2
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
rawPD <- list.files(list.dirs(paste0(path,"proteomicDatafromPD")), pattern = rawPD_pattern, full.names = TRUE)
# rawPeptide <- list.files(list.dirs(paste0(path,"proteomicDatafromPD")), pattern = rawPeptide_pattern, full.names = TRUE)

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

# Function to subset data based on cell number
subset_counts_by_number <- function(cnt, pattern) {
  subset(cnt, select = c("Accession", ends_with(pattern)))
}

# Create a list of data frames for different cell numbers
count.number <- lapply(counts[[2]], function(cnt) {
  lapply(c("1_None$", "2_5cells$", "3_10cells$", "4_50cells$", "5_100cells$"), 
         subset_counts_by_number, cnt = cnt)
})

# Extract unique accession numbers based on cell number
extract_unique_accessions <- function(df_list, pattern) {
  df <- df_list[[1]][[pattern]]
  df$Accession[unique(unlist(lapply(df[-1], function(x) which(!is.na(x)))))]

}

count.types <- map(names(count.number[[1]]), extract_unique_accessions, df_list = count.number)

# Figure 2a. Generate a Venn diagram
ggVennDiagram(count.types) + scale_fill_gradient(low = "blue", high = "red")

# Function to count non-NA values for proteins
count_non_na_proteins <- function(x) {
  x %>% na.omit() %>% length()
}

# Generate data for boxplots
num.proteins <- map_df(c("Number = 0", "Number = 5", "Number = 10", "Number = 100"), 
                       function(n) {
                         data.frame(
                           proteins = apply(count.number[[1]][[n]], 2, count_non_na_proteins)[-1],
                           cell_number = str_replace_all(n, "[^0-9]", "")
                         )
                       })

# Figure 2b. Boxplot for proteins
Figure2b <- ggboxplot(num.proteins, x = "cell_number", y = "proteins",
          palette = "aaas", add = "jitter") +
  stat_compare_means(label = "p.signif", method = "t.test", 
                     comparisons = my_comparisons, label.y = c(2500, 2800, 3100)) +
  labs(x = "Cell Number", y = "") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size = 20, colour = "black"),
        axis.text.y = element_text(size = 20), axis.title = element_text(size = 16),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

# Function to calculate log10 abundance and rank
calc_log10_rank <- function(df, number_label) {
  count.avg <- rowMeans(df[-1], na.rm = TRUE)
  df %>% 
    select(Accession) %>%
    mutate(
      Rank = rank(-count.avg, ties.method = "first"),
      Log10Abundance = log10(count.avg),
      Type = number_label
    )
}

# Prepare data frames for abundance-rank plot
df_list <- map(c("Number = 0", "Number = 5", "Number = 10", "Number = 50", "Number = 100"),
               calc_log10_rank, df_list = count.number)

# Combine all data frames
df <- bind_rows(df_list)

# Figure 2c. Abundance-rank plot
Figure2c <- ggplot(df, aes(x = Rank, y = Log10Abundance, color = Type)) +
  geom_point(size = 1) +
  scale_color_brewer(palette = "Pastel2") +
  theme_bw() +
  theme(legend.position = c(0.85, 0.8), axis.text.x = element_text(size = 20, colour = "black"),
        axis.text.y = element_text(size = 20), axis.title = element_text(size = 16),
        legend.title = element_blank(), legend.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

# Function to process count data and create a data frame with rank and log10 abundance
process_counts <- function(data, number_label) {
  count.avg <- rowMeans(data[,-1], na.rm = TRUE)
  data.frame(
    Accession = data[,1],
    Rank = rank(-count.avg, na.last = TRUE),
    Log10Abundance = log10(as.numeric(count.avg))
  ) %>%
    mutate(NumberLabel = number_label)
}

# Process each count number and store the results in a list
df_list <- lapply(count.number, function(cn_list) {
  lapply(names(cn_list), function(number_label) {
    process_counts(cn_list[[number_label]], number_label)
  })
})

# Combine the results into one data frame and reshape for plotting
df_combined <- do.call(rbind, unlist(df_list, recursive = FALSE))

# Convert from wide to long format for plotting
df_long <- pivot_longer(df_combined, cols = -c(Accession, NumberLabel), 
                        names_to = "CellNumber", values_to = "Value")

# Define comparisons for statistical tests
my_comparisons <- list(c("005","010"), c("005","100"))

# Figure 2d. Create violin plot with statistical comparisons
Figure2d <- ggplot(df_long, aes(x = CellNumber, y = Value, fill = CellNumber)) +
  geom_violin(trim = FALSE) +
  geom_jitter(shape = 16, position = position_jitter(0.25), cex = 0.01) +
  geom_boxplot(width = 0.1) +
  scale_fill_brewer(palette = "Pastel2") +
  xlab("Cell Number") +
  ylab("Log10Abundance") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20, hjust = 0, colour = "black"),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  stat_compare_means(comparisons = my_comparisons, aes(label = after_stat(p.signif)), 
                     method = "t.test", label.y = c(11, 12, 13))


# Function to load data and select necessary columns
load_and_select <- function(path, pattern, select_cols) {
  filenames <- list.files(path, pattern = pattern, full.names = TRUE)
  lapply(filenames, function(file) {
    fread(file) %>% 
      select(all_of(select_cols))
  })
}

# Function to combine data and calculate correlations
calculate_correlations <- function(data_list) {
  combined_data <- reduce(data_list, full_join, by = "Accession")
  correlations <- cor(combined_data[,-1], method = "pearson", use = "pairwise.complete.obs")
  setdiff(as.vector(correlations), 1) %>% unique()
}

# Load data for different cell numbers and calculate correlations
cell_numbers <- c("10", "5", "50", "100")
correlation_list <- setNames(vector("list", length(cell_numbers)), cell_numbers)

for (cn in cell_numbers) {
  pattern <- sprintf("^221229_%s_", cn)
  analysis_list <- load_and_select(path, pattern, c("Accession", "Number of PSMs"))
  correlation_list[[cn]] <- calculate_correlations(analysis_list)
}

# Combine correlations into one data frame
correlation_data <- bind_rows(
  lapply(correlation_list, function(corrs, type) {
    data.frame(Corr = na.omit(corrs), Type = type)
  }, type = names(correlation_list))
)

# Calculate mean and standard deviation for the largest cell number
mean_corr_100 <- mean(correlation_list[["100"]], na.rm = TRUE)
sd_corr_100 <- sd(correlation_list[["100"]], na.rm = TRUE)

# Create a boxplot for the correlations
my_comparisons <- list(c("005","010"), c("005","100"))

# Figure2e. 
Figure2e <- ggplot(correlation_data, aes(x = Type, y = Corr)) +
  geom_jitter(shape = 16, position = position_jitter(0.25), cex = 0.01) +
  geom_boxplot(width = 0.1) +
  stat_compare_means(comparisons = my_comparisons,
                     aes(label = after_stat(p.signif)),
                     method = "t.test",
                     label.y = c(1.1, 1.2, 1.3)) +
  xlab("Cell Number") +
  ylab("Correlation Coeff.") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20, hjust = 0, colour = "black"),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        legend.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
