# ------------------------------------------------------------------------------
# MSc Thesis 
# Dirk Hermans
# Laboratory of Nematology
# Wageningen University & Research
# 2026
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Prerequisites
# ------------------------------------------------------------------------------

# 1. Download, unzip and add following 43 .txt files to ./data/E-MTAB-11658/
#     https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11658

# 2. Place the targets.txt file in directory ./normalization/


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
# Loading data
# ------------------------------------------------------------------------------

# Working directories
dir_home <- "/Users/diher/Repos/wur/thesis_dirk" 
dir_data <- paste(dir_home, "/data/E-MTAB-11658", sep = "")
dir_dataMA <- paste(dir_home, "/data/MA", sep = "")
dir_functions <- paste(dir_home, "/function_scripts", sep = "")
dir_normalized <- paste(dir_home, "/normalized", sep = "")
setwd(dir_home)

# Targets
targets_RIL <- read.delim("./normalization/Targets_RIL.txt", stringsAsFactors=FALSE)
agi.id <- read.delim("./normalization/ArrayID_agilentV2_WS258.txt", stringsAsFactors = FALSE)

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

presentation <- theme(axis.text.x = element_text(size=12, face="bold", color="black"),
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
                      panel.border = element_rect(colour="lightgrey",fill = NA))



# ------------------------------------------------------------------------------
# Data normalization
# ------------------------------------------------------------------------------

rg.norm <- transcriptomics.agilent.norm(targets = targets_RIL,
                save_dir = dir_normalized)

trans.int <- transcriptomics.transform.norm(rg.norm,
                save_dir = dir_dataMA)

# Checks
correlsums <- transcriptomics.check.cor(trans.int, 
                                        save_dir = dir_dataMA)
transcriptomics.check.genes(trans.int,
                            spot.id = agi.id$gene_public_name,
                            save_dir = dir_dataMA)

###Make a list
colnames.sep <- ":"
colnames.names <- c("number","strain","batch","alphasyn","days","sample_number")

list.data <- transcriptomics.list.to.df(trans.int = trans.int, spot.id=agi.id$SpotID,colnames.sep = colnames.sep, colnames.names = colnames.names)
save(list.data,file=paste(dir_home, "/Normalized/obj_list.data.Rdata", sep = ""))



# # Background correction
# bgCorrected <- backgroundCorrect(rawDataMicroArrays, method="normexp", offset=50) # Bodil method="none"
# 
# # Within-array normalization (LOESS)
# maNorm <- normalizeWithinArrays(bgCorrected, method="loess")
# 
# # Between-array normalization (quantile on A-values)
# maNorm <- normalizeBetweenArrays(maNorm, method="quantile")
# 
# # Convert to M-values (log-ratios) for downstream analysis (*)
# M <- maNorm$M
# A <- maNorm$A
# 
# # Inspect summary to check normalization
# summary(M)
# summary(A)


# transcriptomics.agilent.norm
#transcriptomics.transform.norm
#transcriptomics.check.core
#transcriptomics.check.genes
#transcriptomics.list.to.df (error)



# ------------------------------------------------------------------------------
# (*) Explanations on M- and A-values
# ------------------------------------------------------------------------------

# $M — log-ratio (fold-change)
# 
# $M comes from the original red and green intensities:
#   
#   M=log⁡2(R)−log⁡2(G)
# 
# It represents the log-fold change between the two channels for each spot.
# 
# After within-array normalization, $M is adjusted to remove dye-bias, so comparisons across arrays are more reliable.

# 2. $A — average log-intensity
# 
# $A is the mean log-intensity of the two channels:
#   
#   A = 1/2 (log⁡2(R)+log⁡2(G))
# 
# It gives a measure of overall signal strength.
# 
# $A is used in MA-plots, which plot M vs. A to visualize intensity-dependent biases.

# 3. maNorm holds the MA-transformed data, i.e., $M and $A values after within-array normalization.


