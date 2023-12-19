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
