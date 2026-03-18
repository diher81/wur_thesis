# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#    ___  ___ _____        _____ _               _     
#    |  \/  |/  ___|      |_   _| |             (_)    
#    | .  . |\ `--.  ___    | | | |__   ___  ___ _ ___ 
#    | |\/| | `--. \/ __|   | | | '_ \ / _ \/ __| / __|
#    | |  | |/\__/ / (__    | | | | | |  __/\__ \ \__ \
#    \_|  |_/\____/ \___|   \_/ |_| |_|\___||___/_|___/
#    ______ _       _       _                          
#    | ___ (_)     | |     | |                         
#    | |_/ /_  ___ | | ___ | | __ _ _   _              
#    | ___ \ |/ _ \| |/ _ \| |/ _` | | | |             
#    | |_/ / | (_) | | (_) | | (_| | |_| |             
#    \____/|_|\___/|_|\___/|_|\__, |\__, |             
#                              __/ | __/ |             
#                             |___/ |___/                                       
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Dirk Hermans
# Laboratory of Nematology
# Wageningen University & Research
# 2026
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Prerequisites
# ------------------------------------------------------------------------------

# 1. Download, unzip and add following 43 .txt files to ./data/E-MTAB-11658/
#     https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11658



# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------

# Install and load renv
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
library(renv)

# Define required packages
packages <- list(
  cran = c(
    "tidyverse"
  ),
  bioc = c("limma", 
           "statmod")
)

# Ensure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  renv::install("BiocManager")
}

# Install missing packages
for(type in names(packages)) {
  installer <- switch(type, cran = renv::install, bioc = BiocManager::install)
  for(pkg in packages[[type]]) {
    if (!requireNamespace(pkg, quietly = TRUE)) installer(pkg)
  }
}

# Load all libraries
lapply(unlist(packages), library, character.only = TRUE)


# ------------------------------------------------------------------------------
# Instantiate variables
# ------------------------------------------------------------------------------

# Working directories
dir_home <- "/Users/diher/Repos/wur/thesis_dirk/" 
dir_data <- paste(dir_home, "data/E-MTAB-11658/", sep = "")
dir_dataMA <- paste(dir_home, "data/MA/", sep = "")
dir_functions <- paste(dir_home, "function_scripts/", sep = "")
dir_target <- paste(dir_home, "target/", sep = "")
dir_output <- paste(dir_home, "output/", sep = "")
dir_normalized <- paste(dir_home, "normalized/", sep = "")

setwd(dir_home)


# ------------------------------------------------------------------------------
# Loading data
# ------------------------------------------------------------------------------

# Targets
targets_RIL <- read.delim("./target/Targets_RIL.txt", 
                          stringsAsFactors=FALSE)
agi.id <- read.delim("./target/ArrayID_agilentV2_WS258.txt", 
                     stringsAsFactors = FALSE)


# ------------------------------------------------------------------------------
# Initial data inspection - OPTIONAL
# ------------------------------------------------------------------------------

# Read raw micro-array data
microArrayFileNames <- targets_RIL$FileName

# Read raw red/green intensities
rawDataMicroArrays <- read.maimages(microArrayFileNames, source = "agilent")

# A quick glance at the data loaded
dim(rawDataMicroArrays$R)
dim(rawDataMicroArrays$G)
names(rawDataMicroArrays)
summary(as.vector(rawDataMicroArrays$R))
summary(as.vector(rawDataMicroArrays$G))
head(rawDataMicroArrays$genes)


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

setwd(dir_functions)
for(i in 1:length(dir())){
  source(dir()[i])
}
setwd(dir_home)


# ------------------------------------------------------------------------------
# Plotting themes
# ------------------------------------------------------------------------------

presentation <- 
  theme(axis.text.x = element_text(size=12, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.title.y = element_text(size=14, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=16, face="bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour="lightgrey"),
        legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "lightgrey", fill = NA))


# ------------------------------------------------------------------------------
# Data normalization
# ------------------------------------------------------------------------------

# Generate normalized data (matrix) and save the normalized data files
# R = Red channel = Cy5
# G = Green channel = Cy3
rg.normalized.intensities <- normalize.agilent.transcriptomics(targets = targets_RIL,
                save_dir = dir_normalized)

# Generate a list with "log2.rat.mean", "log2.int", "standard_score"
transformed.intensities <- transcriptomics.transform.norm(rg.normalized.intensities,
                save_dir = dir_dataMA)

# Checks
correlsums <- transcriptomics.check.cor(transformed.intensities, save_dir = dir_output)
transcriptomics.check.genes(transformed.intensities,
                            spot.id = agi.id$gene_public_name,
                            save_dir = dir_dataMA)

# Make and save a list of log2 transformed intensities and log2 ratio means
colnames.names <- c("number","strain","batch","alphasyn","days","sample_number")
list.data <- transcriptomics.list.to.df(trans.int = transformed.intensities, 
                                        spot.id=agi.id$SpotID,
                                        colnames.sep = ":", 
                                        colnames.names = colnames.names)
save(list.data, 
     file = paste(dir_normalized, "/obj_list.data.Rdata", sep = ""))

