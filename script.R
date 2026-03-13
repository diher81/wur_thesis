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
# Plotting themes
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Loading data
# ------------------------------------------------------------------------------

# Set working directories
workwd <- getwd()
workwd_data <- "./data/E-MTAB-11658"

# Read raw micro-array data
targets <- read.delim("./normalization/Targets.txt", stringsAsFactors=FALSE)
microArrayFileNames <- targets$FileName

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
# Data normalization
# ------------------------------------------------------------------------------

# Background correction
bgCorrected <- backgroundCorrect(rawDataMicroArrays, method="normexp", offset=50)

# Within-array normalization (LOESS)
maNorm <- normalizeWithinArrays(bgCorrected, method="loess")

# Between-array normalization (quantile on A-values)
maNorm <- normalizeBetweenArrays(maNorm, method="Aquantile")

# Convert to M-values (log-ratios) for downstream analysis (*)
M <- maNorm$M
A <- maNorm$A

# Inspect summary to check normalization
summary(M)
summary(A)




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


