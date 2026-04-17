# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#      __  __  _____        _______ _               _     
#     |  \/  |/ ____|      |__   __| |             (_)    
#     | \  / | (___   ___     | |  | |__   ___  ___ _ ___ 
#     | |\/| |\___ \ / __|    | |  | '_ \ / _ \/ __| / __|
#     | |  | |____) | (__     | |  | | | |  __/\__ \ \__ \
#     |_|  |_|_____/ \___|    |_|  |_| |_|\___||___/_|___/
#      _ __  _       _                                    
#     |  _ \(_)     | |                                   
#     | |_) |_  ___ | | ___   __ _ _   _                  
#     |  _ <| |/ _ \| |/ _ \ / _` | | | |                 
#     | |_) | | (_) | | (_) | (_| | |_| |                 
#     |____/|_|\___/|_|\___/ \__, |\__, |                 
#                             __/ | __/ |                 
#                            |___/ |___/                                       
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
#       https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11658
#
# 2. Place following files in ./data/target/ : 
#       "Targets_RIL.txt" and "ArrayID_agilentV2_WS258.txt"
#
# 3. Place following files in ./data/Genetic_map/ : 
#       "asRIL_map_new.txt" and "asRIL_map_new.txt"
#
# 4. Place following files in ./data/proteinAccumulation/ : 
#       "obj_elisa_aS.Rdata" and "obj_qpcr_aS.Rdata"
#
# 5. Place following files in ./data/lifespan/ : 
#       "obj_life_data_stats.Rdata"
#
# 6. Place following files in ./data/QTL/ :
#      "aS.simulation.RData" and "obj_peak.aS.eQTL.Rdata"
#      (unless you want to calculate them again using time-consuming 
#       method calls in this script.)


# ------------------------------------------------------------------------------
# General info
# ------------------------------------------------------------------------------

# ---- Setup reproducible environment ----
# When running this script on a new machine, run these lines first:
###if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
###renv::activate()
###renv::restore()

# ---- Adding packages ----
# To install new packages, run these lines once:
# install.packages("newPackage")
# renv::snapshot()
# Then add library in this script:
# library(newPackage)


# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------

library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(httr)
library(jsonlite)
library(limma)
library(statmod)
library(here)


# ------------------------------------------------------------------------------
# Set working directories
# ------------------------------------------------------------------------------

root <- here::here()
paths <- list(
  R = file.path("R"),
  scripts = file.path("scripts"),
  data = list(
    arrays = file.path("data", "E-MTAB-11658"),
    geneticMap = file.path("data", "Genetic_map"),
    lifespan = file.path("data", "lifespan"),
    protAcc = file.path("data", "proteinAccumulation"),
    qtl = file.path("data", "QTL"),
    target = file.path("data", "target")
  ),
  output = list(
    elisa = file.path("output", "elisa"),
    eqtl = file.path("output", "eqtl"),
    fdr = file.path("output", "FDR"),
    geneticMap = file.path("output", "Genetic_map"),
    genes = file.path("output", "genes"),
    ma = file.path("output", "MA"),
    lifespan = file.path("output", "lifespan"),
    normalized = file.path("output", "normalized"),
    qpcr = file.path("output", "qpcr")
  )
)

setwd(root)

# ---- resolve full paths ----
data_dirs <- file.path(root, unlist(paths$data, recursive = TRUE, use.names = FALSE))
output_dirs <- file.path(root, unlist(paths$output, recursive = TRUE, use.names = FALSE))

# ---- create data directories ----
for (d in data_dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    message("Created data directory: ", d)
  }
}

# ---- create output directories ----
for (d in output_dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    message("Created output directory: ", d)
  }
}

message("Project setup complete.")


# ------------------------------------------------------------------------------
# Loading data and functions
# ------------------------------------------------------------------------------

# ---- Load my data ----

# Targets
targets_RIL <- read.delim(file.path(paths$data$target, "Targets_RIL.txt"),
                          stringsAsFactors = FALSE)
agi.id <- read.delim(file.path(paths$data$target, "ArrayID_agilentV2_WS258.txt"), 
                     stringsAsFactors = FALSE)

# load DNAseq map
populationMap <- data.matrix(read.table(
  file.path(paths$data$geneticMap, "asRIL_map_new.txt"))[,-c(1:3,5,6,8,9,11,13)])
populationMarkers <- read.table(
  file.path(paths$data$geneticMap, "asRIL_map_new.txt"))[,c(1:3)]

# Protein accumulation
load(file = file.path(paths$data$protAcc, "obj_elisa_aS.Rdata"))
load(file = file.path(paths$data$protAcc, "obj_qpcr_aS.Rdata"))

# Lifespan
load(file = file.path(paths$data$lifespan, "obj_life_data_stats.Rdata"))

# ---- Load my functions ----
files <- list.files(paths$R, pattern = "\\.R$", full.names = TRUE)
invisible(lapply(files, source))


# ------------------------------------------------------------------------------
# Initial data inspection - OPTIONAL
# ------------------------------------------------------------------------------

inspectData = FALSE

if(inspectData){
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
}


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

blankTheme <- theme(plot.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks = element_blank())

myColors <- c(brewer.pal(9,"Set1")[c(2,5)],
              brewer.pal(9,"Purples")[c(4,6,6)],
              "black",
              "darkgrey",
              "black",
              brewer.pal(12,"Paired")[c(1,7)],
              rep(brewer.pal(12,"Paired")[c(8,2,7,1)], times = 3))
fillScale <- scale_fill_manual(name = "Treatment", values = myColors)
colScale <- scale_colour_manual(name = "Treatment", values = myColors)


# ------------------------------------------------------------------------------
# Data normalization
# ------------------------------------------------------------------------------

# Generate normalized data (matrix) and save the normalized data files
# R = Red channel = Cy5
# G = Green channel = Cy3
rg.normalized.intensities <- normalize.agilent.transcriptomics(
  targets = targets_RIL, save_dir = file.path(root, paths$output$normalized))

# Generate a list with "log2.rat.mean", "log2.int", "standard_score"
transformed.intensities <- transcriptomics.transform.norm(
  rg.normalized.intensities, save_dir = file.path(root, paths$output$ma))

# Checks
correlsums <- transcriptomics.check.cor(transformed.intensities, 
                                        save_dir = file.path(root, paths$output$normalized))
transcriptomics.check.genes(transformed.intensities,
                            spot.id = agi.id$gene_public_name,
                            save_dir = file.path(root, paths$output$ma))

# Make and save a list of log2 transformed intensities and log2 ratio means
colnames.names <- c("number","strain","batch","alphasyn","days","sample_number")
list.data.RIL <- transcriptomics.list.to.df(trans.int = transformed.intensities, 
                                        spot.id = agi.id$SpotID,
                                        colnames.sep = ":", 
                                        colnames.names = colnames.names)
save(list.data.RIL, 
     file = file.path(root, paths$output$normalized, "/obj_list.data.RIL.Rdata"))


# ------------------------------------------------------------------------------
# Genetic map
# ------------------------------------------------------------------------------

# Figure 1: Genetic map
genetic_map <- populationMap;
genetic_map[is.na(genetic_map)] <- 0

# rewrite the genetic map to a dataframe plottable in ggplot
data.plot <- prep.ggplot.genetic.map(genetic_map, populationMarkers) %>%
    dplyr::group_by(strain) %>%
    dplyr::mutate(genotype = ifelse(genotype == 1, "N2", 
                                    ifelse(genotype == -1, "CB4856", NA))) %>%
    rename(start = position_start, end = position_end, gt = genotype, chrom = chromosome) %>%
    data.frame() %>%
    mutate(strain = as.factor(strain)) %>%
    group_by(strain) %>%
    dplyr::filter(!is.na(gt))

data.plot$index <- dplyr::group_indices(data.plot)
strain_index <- data.plot$strain
names(strain_index) <- data.plot$index + 0.5    

figure_1_gen_map <- ggplot(data.plot, aes(xmin = start, xmax = end, ymin = index, ymax = index + 1, fill = gt)) +
  geom_rect() + 
  scale_fill_manual(values = c("#0080FF","#FF8000",NA)) +
  facet_grid(.~chrom, scales="free", space="free") +
  presentation + 
  labs(x = "Chromosome position (Mb)") + ylab("Strain") +
  ggtitle("Figure 1: Genetic map") +
  scale_x_continuous(labels = function(x) { x/1e6 }, expand = c(0,0)) +
  scale_y_continuous(breaks = unique(data.plot$index) + 0.5, 
                     labels = function(x) {
                       strain_index[as.character(x)] 
                     }, 
                     expand = c(0,0)) +
  theme(strip.background = element_blank(), legend.position = "None", panel.spacing = unit(0,"lines"))
print(figure_1_gen_map)

# Figure 2: centimorgan map
strain.map <- populationMap[,-c(1:4)]
strain.map[strain.map==0] <- NA

# Calculate genetic distance between markers
data.plot <- genetic.distance(strain.map, populationMarkers)

figure_2_cm_map <- ggplot(data.plot, aes(x=position,y=cM)) +
  geom_line() + geom_rug() + facet_grid(~chromosome,space = "free_x", scales = "free_x") + 
  presentation + labs(x = "Physical position (Mb)", y = "Genetic position (cM)") +
  ggtitle("Figure 2: centimorgan map") +
  scale_x_continuous(breaks=c(5,10,15,20)*10^6, labels = c(5,10,15,20))
print(figure_2_cm_map)

# Figure 3: Genotype distribution
data.plot <- cbind(populationMarkers, populationMap) %>%
  pivot_longer(
    cols = -c(Name, Chromosome, Position),
    names_to = "strain",
    values_to = "genotype"
  ) %>%
  group_by(Chromosome, Position, genotype) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    genotype = case_when(
      genotype == -1 ~ "CB4856",
      genotype == 1  ~ "N2",
      TRUE           ~ "NA"
    ),
    genotype = factor(genotype, levels = c("NA", "CB4856", "N2"))
  )

figure_3_gen_distr <- ggplot(data.plot, aes(x = Position, y = n, fill = genotype)) +
  geom_area() + facet_grid(~Chromosome, space = "free_x", scales="free_x") + fillScale +
  presentation + labs(x = "Position (Mb)", y = "Genotype (n strains)") +
  ggtitle("Figure 3: Genotype distribution") +
  scale_x_continuous(breaks = c(5,10,15,20)*10^6, 
                     labels = c(5,10,15,20)) + scale_y_continuous(expand = c(0,0)) + guides(fill = guide_legend(title = "Genotype"))
print(figure_3_gen_distr)

# Figure 4: Genetic map
blank.plot <- ggplot() + geom_blank(aes(1,1)) + blankTheme
layout <- rbind(c(1,2),c(1,3),c(1,4))

annotation.grobA <- title.grob <- grid::textGrob(label = "A",
                                           x = grid::unit(0, "lines"),
                                           y = grid::unit(0, "lines"),
                                           hjust = 0, 
                                           vjust = 0,
                                           gp = grid::gpar(fontsize = 20,fontface="bold"))
annotation.grobB <- title.grob <- grid::textGrob(label = "B",
                                           x = grid::unit(0, "lines"),
                                           y = grid::unit(0, "lines"),
                                           hjust = 0, 
                                           vjust = 0,
                                           gp = grid::gpar(fontsize = 20,fontface="bold"))
annotation.grobC <- title.grob <- grid::textGrob(label = "C",
                                           x = grid::unit(0, "lines"),
                                           y = grid::unit(0, "lines"),
                                           hjust = 0, 
                                           vjust = 0,
                                           gp = grid::gpar(fontsize = 20,fontface="bold"))

graph.oject.gen.map <- arrangeGrob(figure_1_gen_map, top = annotation.grobA)
graph.oject.cm.map  <- arrangeGrob(figure_2_cm_map, top = annotation.grobB)
graph.oject.gen.distr  <- arrangeGrob(figure_3_gen_distr, top = annotation.grobC)

# Open pdf device
pdf(file = file.path(root, paths$output$geneticMap, "figure4-Genetic_map.pdf"), width = 10, height = 14)
# Draw to pdf
grid.arrange(graph.oject.gen.map, graph.oject.cm.map, graph.oject.gen.distr, 
             blank.plot, layout_matrix = layout, heights = c(3,3,8))
# Close pdf
dev.off()


# ------------------------------------------------------------------------------
# eQTL
# ------------------------------------------------------------------------------

# Figure 5: RILs used for eQTL mapping
popmap_in <- populationMap
popmap_in[is.na(popmap_in)] <- 0
data.plot.RIL <- prep.ggplot.genetic.map(popmap_in, populationMarkers) %>%
  dplyr::filter(strain %in% list.data.RIL$strain, 
                !strain %in% c("N2", "CB4856", "NL5901", "SCH4856")) %>%
  dplyr::mutate(index = as.numeric(as.factor(strain)),
                genotype_name = ifelse(genotype == -1, 
                                       "CB4856", 
                                       ifelse(genotype == 1, "N2", NA))) %>%
  dplyr::filter(!is.na(genotype_name))

fig5_RILs.for.eQTL.mapping <- ggplot(
  data.plot.RIL,
  aes(xmin = position_start,
      xmax = position_end,
      ymin = index,
      ymax = index + 1,
      fill = genotype_name)) +
  geom_rect(aes(alpha = 0.8)) +
  facet_grid(. ~ chromosome, scales = "free_x", space = "free_x") +
  fillScale +
  presentation +
  theme(
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  scale_y_continuous(
    breaks = unique(data.plot.RIL$index) + 0.5,
    labels = unique(data.plot.RIL$strain),
    expand = c(0, 0.5)) + 
  scale_x_continuous(
    breaks = c(5, 10, 15, 20) * 10^6,
    labels = c(5, 10, 15, 20),
    expand = c(0, 0.5)
  ) +
  xlab("Genomic location (Mbp)") +
  ylab("Strain") +
  ggtitle("Figure 5: RILs used for eQTL mapping") 

fig5_RILs.for.eQTL.mapping

# eQTL mapping
data.eQTL <- filter(list.data.RIL, !strain %in% c("CB4856", "SCH4856", "N2", "NL5901")) %>%
  dplyr::select(SpotID, strain, log2_intensities) %>%
  spread(key = strain, value = log2_intensities) %>%
  tibble::column_to_rownames("SpotID")

data.eQTL <- data.matrix(data.eQTL)

# Data preparation for QTL mapping
# Returns a list with the entries Trait (log2 intensities), Map, and Marker, 
# where the Traits are aligned with the map
# AGIWURxxxx = SpotIDs
data.eQTL <- QTL.data.prep(data.eQTL, 
                           colnames(data.eQTL), 
                           populationMap, 
                           populationMarkers)
lapply(data.eQTL, head)  

# function for single marker mapping
# Generate a list with names: LOD, Effect, Trait, Map, and Marker.
aS.eQTL <- QTL.map.1(data.eQTL[[1]], data.eQTL[[2]], data.eQTL[[3]])
save(aS.eQTL, file = file.path(root, paths$output$eqtl, "/obj_aS.eQTL.Rdata"))
str(aS.eQTL)

# Build eQTL list file with calculated threshold value 4.3
# These lines are not executable on workstations due to memory limits
executeLongMethod = FALSE
if(executeLongMethod){
  peak.aS.eQTL <- QTL.map1.dataframe(map1.output = aS.eQTL) %>%
    QTL.peak.finder(threshold = 4.3)
  save(peak.aS.eQTL, file = file.path(root, paths$output$eqtl, "obj_peak.aS.eQTL.Rdata"))
  load(file.path(root, paths$output$eqtl, "obj_peak.aS.eQTL.Rdata"))
} else {
  # Instead of generating: load previously generated peak.aS.eQTL file here:
  load(file.path(root, paths$data$qtl, "obj_peak.aS.eQTL.Rdata"))
}

# Create eQTL table
aS.eQTL.table <- QTL.eQTL.table(
  QTL.peak.dataframe = peak.aS.eQTL,
  cis.trans.window = "LOD.drop",
  trait.annotation = cbind(trait = agi.id[,1],
                           dplyr::select(agi.id,
                                        chromosome, 
                                        gene_bp_start, 
                                        gene_WBID, 
                                        gene_sequence_name, 
                                        gene_public_name))) %>%
  QTL.table.addR2(aS.eQTL) %>%
  group_by(gene_public_name,qtl_chromosome) %>%
  mutate(qtl_representative=ifelse(qtl_R2_sm == max(qtl_R2_sm),"yes","no")) %>%
  data.frame() %>%
  filter(!is.na(gene_bp), qtl_representative == "yes") %>%
  arrange(trait); head(aS.eQTL.table)

# Add transband
call.transbands.file <- QTL.eQTL.call.TBs(aS.eQTL.table, 
                                          window = 0.5e6,
                                          chromosomes = "III",
                                          chromosome_size = 13783801)

# This function finds the number of successes that corresponds to a certain percentile 
# based on an average rate of success. 
# (The corresponding Poisson quantiles for a set of probabilities are obtained.)
qpois(0.0001, 18.36, lower.tail = F)

# For every QTL on chromosome III, indicate whether it is located in a trans-band or cis-enriched area
aS.eQTL.table <-  QTL.eQTL.table.addTB(aS.eQTL.table,
                                       call.transbands.file,
                                       merge_condition = 1)

# Save
save(aS.eQTL.table, file = file.path(root, paths$output$eqtl, "obj_aS.eQTL.table.Rdata"))         


# ------------------------------------------------------------------------------
# power analysis (R only step)
# ------------------------------------------------------------------------------

# Conduct mapping simulation to investigate power
# +/- 35 min
executeLongMethod = FALSE
if(executeLongMethod){
  aS.simulation <- QTL.map.1.sim(strain.map = data.eQTL$Map,
                                 strain.marker = data.eQTL$Marker,
                                 threshold = 4.3,
                                 n_per_marker = 10)
  
  aS.simulation
  save(aS.simulation, file = file.path(root, paths$output$eqtl, "aS.simulation.RData"))
  load(file.path(root, paths$output$eqtl, "aS.simulation.RData"))
} else {
  # Instead of generating: load previously generated aS.simulation file here:
  load(file.path(root, paths$data$qtl, "aS.simulation.RData"))
}

writexl::write_xlsx(aS.eQTL.table,
                    path = file.path(root, paths$output$eqtl, "Supplementary_table5-eQTL.xlsx"))

# Make a table with eQTL peaks for plotting.
# This function adds points to faithfully indicate the chromosome boundaries (standard is C. elegans))
plot.data <- prep.ggplot.eQTL.table(aS.eQTL.table) 

figure_6_locations <- ggplot(plot.data, aes(x = qtl_bp,y = gene_bp)) +
  geom_segment(aes(x = qtl_bp_left,
                   y = gene_bp,
                   xend = qtl_bp_right,
                   yend = gene_bp),
               alpha = 0.25,colour = "darkgrey") + 
  geom_point(aes(colour  =  qtl_type),
             size = 0.6, alpha = 0.5)  +
  facet_grid(gene_chromosome ~ qtl_chromosome, 
             space  =  "free",
             scales = "free") + 
  presentation + colScale + 
  theme(legend.position  =  "none") +
  labs(x = "eQTL peak position (Mb)",
       y = "Gene position (Mb)") + 
  ggtitle("Figure 6: eQTL locations") +
  theme(panel.spacing = unit(0.1,"lines")) +
  scale_x_continuous(breaks = c(5,10,15,20)*10^6, labels = c(5,10,15,20)) +
  scale_y_continuous(breaks = c(5,10,15,20)*10^6, labels = c(5,10,15,20))

figure_6_locations

plot.data <- filter(aS.eQTL.table, qtl_type == "trans") %>%
  prep.ggplot.eQTL.table() 

figure_7_eQTL_counts <- ggplot(plot.data, aes(x = qtl_bp,fill = qtl_type)) +
  geom_histogram(binwidth = 500000) + geom_hline(yintercept = 36,lty = 2,lwd = 1.2) +
  facet_grid(c("")~qtl_chromosome, space  =  "free",scales = "free") + presentation + 
  fillScale + theme(legend.position  =  "none",
                    panel.spacing = unit(0.1,"lines")) +
  labs(x = "eQTL peak position (Mb)", y = "eQTL counts") +
  ggtitle("Figure 7: eQTL") +
  scale_x_continuous(breaks = c(5,10,15,20)*10^6,
                     labels = c(5,10,15,20))

figure_7_eQTL_counts


# ------------------------------------------------------------------------------
# Filter out genes of interest
# ------------------------------------------------------------------------------

# Filter out genes between 1M bp and 2.5M bp.
aS.eQTL.table <- aS.eQTL.table %>%
  dplyr::filter(qtl_chromosome == "III",
                dplyr::between(qtl_bp, 1000000, 2500000))

# Get a list of selected gene names and 
genes <- unique(aS.eQTL.table$gene_public_name)
genes <- genes[!is.na(genes)]
cat(genes, sep = "\n")
# copy & paste the result in https://wormbase.org/tools/mine/simplemine.cgi to
# generate an overview of their functions. Save resulting html file in ./output

# Alternative: Use MyGene api to store gene name in var geneInfo:
res <- POST(
  "https://mygene.info/v3/query",
  body = list(
    q = genes,                 
    scopes = "symbol",
    fields = "symbol,name,summary",
    species = "6239"
  ),
  encode = "json"
)

geneInfo <- fromJSON(content(res, "text", encoding = "UTF-8"))
geneInfo <- geneInfo[, c("query", "name")]
geneInfo

# Write these genes to a .pdf file
df <- geneInfo
colnames(df) <- c("Gene", "Description")
df <- geneInfo
colnames(df) <- c("Gene", "Description")
df$Description <- str_wrap(df$Description, width = 50)
pdf(file.path(paths$output$genes, "geneInfo.pdf"), width = 8.5, height = 11)
rows_per_page <- 35
n <- nrow(df)
tt <- ttheme_default(base_size = 8)
for (i in seq(1, n, by = rows_per_page)) {
  grid.table(
    df[i:min(i + rows_per_page - 1, n), ],
    theme = tt
  )
  if (i + rows_per_page <= n) grid.newpage()
}
dev.off()

# ------------------------------------------------------------------------------
# Draw eQTL profiles and boxplots for genes of interest
# ------------------------------------------------------------------------------

# Gene names and spotIds of the selected genes
spotIds <- aS.eQTL.table[[1]]
geneNames <- aS.eQTL.table[[16]]

# Open pdf device
pdf(file = file.path(root, paths$output$genes, "boxplotsForGenes.pdf"), width = 10, height = 14)

# Create an aQTL profile and a boxplot for every selected gene
for (i in seq_along(spotIds)) {
  message(paste(i, "of", length(spotIds)))
  
  data.plot <- prep.ggplot.QTL.profile(peak.aS.eQTL, aS.eQTL, spotIds[i])
  data.plot[[2]] <- mutate(data.plot[[2]], geno_strain = ifelse(genotype == -1, "CB4856", "N2"))

  plotEqtlProfile <- ggplot(data.plot$QTL_profile, aes(x = qtl_bp,y = qtl_significance,alpha = 0.2)) +
    geom_line(size = 1.5,colour = brewer.pal(9,"Set1")[4]) + 
    facet_grid(~qtl_chromosome,scales = "free",
               space = "free_x") + 
    presentation + 
    theme(legend.position = "none",
          plot.margin = margin(10, 30, 10, 30)) +
    geom_abline(intercept = 4.3, slope = 0, linetype = 2, size = 1) + 
    labs(x = "QTL position (Mb)",
         y = expression(bold("significance"~(-log[10](p)))),
         parse = TRUE) +
    scale_x_continuous(breaks = c(0,10,20)*10^6,labels = c(0,10,20)) + ylim(0,5.5)
  
  plotBoxplot <- ggplot(data.plot[[2]], aes(x = geno_strain, y = trait_value)) +
    geom_jitter(height = 0,
                width = 0.25,
                aes(colour = geno_strain),
                alpha = 0.2) + 
    geom_boxplot(outlier.shape = NA, alpha = 0.2, aes(fill = geno_strain)) +
    labs(x = "Genotype at marker", 
         y = paste(geneNames[i], " expression"),
         parse = TRUE) + 
    facet_grid(~Chromosome+Position) +
    presentation + 
    colScale + 
    fillScale + 
    theme(legend.position = "none",
          plot.margin = margin(10, 30, 10, 30)) +
    annotate("text",
             x=1.5,
             y=max(data.plot[[2]]$trait_value,na.rm=T),
             label = paste0("italic(R)^{2}==",
                         round(data.plot[[2]]$R_squared[1], digits=2)),
             parse = TRUE)
  
  globalTitle <- grid::textGrob(
    label = paste("Gene", geneNames[i], ": eQTL profile and genotype split-out."),
    gp = grid::gpar(fontsize = 22, fontface = "bold"),
    hjust = 0.5
  )

  # Draw to pdf
  grid.arrange(globalTitle,
               plotEqtlProfile, 
               grid::nullGrob(),
               plotBoxplot,
               grid::nullGrob(),
               ncol = 1,
               heights = c(1, 5, 1.5, 5, 1))
}

# Close pdf
dev.off()


# ------------------------------------------------------------------------------
# Protein accumulation - ELISA
# ------------------------------------------------------------------------------

# pQTL mapping
data.pQTL <- filter(elis_data, !Strain %in% c("SCH4856", "NL5901")) %>%
  dplyr::select(Strain, Norm_NLSCH) %>%
  spread(key = Strain, value = Norm_NLSCH)

data.pQTL <- data.matrix(data.pQTL)
rownames(data.pQTL) <- "elisa"

# Data preparation for QTL mapping
data.pQTL <- QTL.data.prep(data.pQTL,
                           colnames(data.pQTL),
                           populationMap,
                           populationMarkers)
lapply(data.pQTL, head)  

# function for single marker mapping
aS.pQTL <- QTL.map.1(data.pQTL[[1]], data.pQTL[[2]], data.pQTL[[3]])
save(aS.pQTL, file = file.path(paths$output$elisa, "/obj_aS.pQTL.Rdata"))

# Permutation for single marker mapping
elisa.QTL.perm <- QTL.map.1.perm(aS.pQTL[[1]], aS.pQTL[[2]], aS.pQTL[[3]], 1000)
save(elisa.QTL.perm, file = file.path(paths$output$elisa, "obj_elisa.QTL.perm.Rdata"))

# FDR calculation for single marker mapping
elisa.FDR <- QTL.map.1.FDR(map1.output = aS.pQTL,
                           filenames.perm = "/obj_elisa.QTL.perm.Rdata",
                           FDR_dir = file.path(root, paths$output$elisa),
                           q.value = 0.025,
                           small = TRUE)
thresholdElisa <- elisa.FDR[[1]] / 10 # 3.21
save(elisa.FDR, file = file.path(paths$output$elisa, "obj_RIL.elis.FDR.Rdata"))

peak.aS.pQTL <- QTL.map1.dataframe(map1.output = aS.pQTL) %>%
  QTL.peak.finder(threshold = thresholdElisa)
save(peak.aS.pQTL, file = file.path(paths$output$elisa, "obj_peak.aS.QTL.Rdata"))

# Plot peak.aS.pQTL
plotPqtlProfile <- ggplot(peak.aS.pQTL, aes(x = qtl_bp, y = qtl_significance, alpha = 0.2)) +
  geom_line(size = 1.5, colour = brewer.pal(9,"Set1")[4]) + 
  facet_grid(~qtl_chromosome, scales = "free",
             space = "free_x") + 
  presentation + 
  theme(legend.position = "none",
        plot.margin = margin(10, 30, 10, 30)) +
  geom_abline(intercept = thresholdElisa, slope = 0, linetype = 2, size = 1) + 
  labs(x = "QTL position (Mb)",
       y = expression(bold("QTL significance")),
       parse = TRUE) +
  scale_x_continuous(breaks = c(0,10,20)*10^6, labels = c(0,10,20)) + ylim(0,5.5) +
  ggtitle("Figure 8: pQTL mapping - ELISA")
print(plotPqtlProfile)


# ------------------------------------------------------------------------------
# Protein accumulation - QPCR
# ------------------------------------------------------------------------------

# pQTL mapping
data.pQTL <- qpcr_aS %>%
  dplyr::select(Strain, mean) %>%
  spread(key = Strain, value = mean)

data.pQTL <- data.matrix(data.pQTL)
rownames(data.pQTL) <- "qpcr"

# Data preparation for QTL mapping
data.pQTL <- QTL.data.prep(data.pQTL,
                           colnames(data.pQTL),
                           populationMap,
                           populationMarkers)
lapply(data.pQTL, head)  

# function for single marker mapping
qpcr.pQTL <- QTL.map.1(data.pQTL[[1]], data.pQTL[[2]], data.pQTL[[3]])
save(qpcr.pQTL, file = file.path(paths$output$qpcr, "obj_qpcr.pQTL.Rdata"))

# Permutation for single marker mapping
qpcr.QTL.perm <- QTL.map.1.perm(qpcr.pQTL[[1]], qpcr.pQTL[[2]], qpcr.pQTL[[3]], 1000)
save(qpcr.QTL.perm, file = file.path(paths$output$qpcr, "obj_qpcr.pQTL.perm.Rdata"))

# FDR calculation for single marker mapping
qpcr.FDR <- QTL.map.1.FDR(map1.output = qpcr.pQTL,
                          filenames.perm = "/obj_qpcr.pQTL.perm.Rdata",
                          FDR_dir = file.path(root, paths$output$qpcr),
                          q.value = 0.025,
                          small = TRUE)
# thresholdQpcr <- qpcr.FDR[[1]] / 10 # 1.67 -> this can't be right
thresholdQpcr <- 3.1
save(qpcr.FDR, file = file.path(paths$output$qpcr, "obj_RIL.qpcr.FDR.Rdata"))

peak.qpcr.pQTL <- QTL.map1.dataframe(map1.output = qpcr.pQTL) %>%
  QTL.peak.finder(threshold = thresholdQpcr)
save(peak.qpcr.pQTL, file = file.path(paths$output$qpcr, "obj_peak.qpcr.pQTL.Rdata"))

# Plot peak.qpcr.pQTL
plotQpcrProfile <- ggplot(peak.qpcr.pQTL, aes(x = qtl_bp, y = qtl_significance, alpha = 0.2)) +
  geom_line(size = 1.5,colour = brewer.pal(9,"Set1")[4]) + 
  facet_grid(~qtl_chromosome, scales = "free",
             space = "free_x") + 
  presentation + 
  theme(legend.position = "none",
        plot.margin = margin(10, 30, 10, 30)) +
  geom_abline(intercept = thresholdQpcr, slope = 0, linetype = 2, size = 1) + 
  labs(x = "QTL position (Mb)",
       y = expression(bold("QTL significance")),
       parse = TRUE) +
  scale_x_continuous(breaks = c(0,10,20)*10^6, labels = c(0,10,20)) + ylim(0,5.5) +
  ggtitle("Figure 9: QTL mapping for qPCR")
print(plotQpcrProfile)


# ------------------------------------------------------------------------------
# Lifespan
# ------------------------------------------------------------------------------

# QTL mapping
traits <- unique(life_data_stats$trait)
blank.plot <- ggplot() + geom_blank(aes(1,1)) + blankTheme

# Open pdf device
pdf(file = file.path(root, paths$output$lifespan, "lifespanQTLmapping.pdf"), width = 10, height = 14)

for (tr in traits) {
  data.QTL <- life_data_stats %>%
    dplyr::filter(!Strain %in% c("N2", "CB4856", "SCH4856", "NL5901"),
                  trait == tr) %>%
    dplyr::select(Strain, Treatment, value) %>%
    tidyr::pivot_wider(id_cols = Treatment,
                       names_from = Strain,
                       values_from = value)
  
  # Split data per treatment
  data_DR <- data.QTL %>%
    dplyr::filter(Treatment == "DR") %>%
    dplyr::select(-Treatment)
  data_DR <- data.matrix(data_DR)
  rownames(data_DR) <- tr
  
  data_NGM <- data.QTL %>%
    dplyr::filter(Treatment == "NGM") %>%
    dplyr::select(-Treatment)
  data_NGM <- data.matrix(data_NGM)
  rownames(data_NGM) <- tr
  
  data_pl <- data.QTL %>%
    dplyr::filter(Treatment == "Plasticity") %>%
    dplyr::select(-Treatment)
  data_pl <- data.matrix(data_pl)
  rownames(data_pl) <- tr
  
  # Data preparation for QTL mapping
  data_DR <- QTL.data.prep(data_DR,
                           colnames(data_DR),
                           populationMap,
                           populationMarkers)
  data_NGM <- QTL.data.prep(data_NGM,
                            colnames(data_NGM),
                            populationMap,
                            populationMarkers)
  data_pl <- QTL.data.prep(data_pl,
                           colnames(data_pl),
                           populationMap,
                           populationMarkers)
  
  # function for single marker mapping
  dr.QTL <- QTL.map.1(data_DR[[1]], data_DR[[2]], data_DR[[3]])
  save(dr.QTL, file = file.path(paths$output$lifespan, "obj_dr.QTL.Rdata"))
  ngm.QTL <- QTL.map.1(data_NGM[[1]], data_NGM[[2]], data_NGM[[3]])
  save(ngm.QTL, file = file.path(paths$output$lifespan, "obj_ngm.QTL.Rdata"))
  pl.QTL <- QTL.map.1(data_pl[[1]], data_pl[[2]], data_pl[[3]])
  save(pl.QTL, file = file.path(paths$output$lifespan, "obj_pl.QTL.Rdata"))
  
  # Permutation for single marker mapping
  dr.QTL.perm <- QTL.map.1.perm(dr.QTL[[1]], dr.QTL[[2]], dr.QTL[[3]], 1000)
  save(dr.QTL.perm,
       file = file.path(paths$output$lifespan, "obj_dr.QTL.perm.Rdata"))
  ngm.QTL.perm <- QTL.map.1.perm(ngm.QTL[[1]], ngm.QTL[[2]], ngm.QTL[[3]], 1000)
  save(ngm.QTL.perm,
       file = file.path(paths$output$lifespan, "obj_ngm.QTL.perm.Rdata"))
  pl.QTL.perm <- QTL.map.1.perm(pl.QTL[[1]], pl.QTL[[2]], pl.QTL[[3]], 1000)
  save(pl.QTL.perm,
       file = file.path(paths$output$lifespan, "obj_pl.QTL.perm.Rdata"))
  
  # FDR calculation for single marker mapping
  dr.FDR <- QTL.map.1.FDR(
    map1.output = dr.QTL,
    filenames.perm = "/obj_dr.QTL.perm.Rdata",
    FDR_dir = file.path(root, paths$output$lifespan),
    q.value = 0.025,
    small = TRUE
  )
  thresholdDr <- dr.FDR[[1]] / 10 # 2.5
  save(dr.FDR,
       file = file.path(paths$output$lifespan, "obj_RIL.DR.FDR.Rdata"))
  peak.dr.QTL <- QTL.map1.dataframe(map1.output = dr.QTL) %>%
    QTL.peak.finder(threshold = thresholdDr)
  save(peak.dr.QTL,
       file = file.path(paths$output$lifespan, "obj_peak.dr.QTL.Rdata"))
  
  ngm.FDR <- QTL.map.1.FDR(
    map1.output = ngm.QTL,
    filenames.perm = "/obj_ngm.QTL.perm.Rdata",
    FDR_dir = file.path(root, paths$output$lifespan),
    q.value = 0.025,
    small = TRUE
  )
  thresholdNgm <- ngm.FDR[[1]] / 10 # 3.84
  save(ngm.FDR,
       file = file.path(paths$output$lifespan, "obj_RIL.ngm.FDR.Rdata"))
  peak.ngm.QTL <- QTL.map1.dataframe(map1.output = ngm.QTL) %>%
    QTL.peak.finder(threshold = thresholdNgm)
  save(peak.ngm.QTL,
       file = file.path(paths$output$lifespan, "obj_peak.ngm.QTL.Rdata"))
  
  pl.FDR <- QTL.map.1.FDR(
    map1.output = pl.QTL,
    filenames.perm = "/obj_pl.QTL.perm.Rdata",
    FDR_dir = file.path(root, paths$output$lifespan),
    q.value = 0.025,
    small = TRUE
  )
  thresholdPl <- pl.FDR[[1]] / 10 # 2.66
  save(pl.FDR,
       file = file.path(paths$output$lifespan, "obj_RIL.pl.FDR.Rdata"))
  peak.pl.QTL <- QTL.map1.dataframe(map1.output = pl.QTL) %>%
    QTL.peak.finder(threshold = thresholdPl)
  save(peak.pl.QTL,
       file = file.path(paths$output$lifespan, "obj_peak.pl.QTL.Rdata"))
  
  # Plot peak.dr.pQTL
  plotLifespanDr <- ggplot(peak.dr.QTL, aes(x = qtl_bp, y = qtl_significance, alpha = 0.2)) +
    geom_line(size = 1.5, colour = brewer.pal(9, "Set1")[4]) +
    facet_grid( ~ qtl_chromosome, scales = "free", space = "free_x") +
    presentation +
    theme(legend.position = "none",
          plot.margin = margin(10, 30, 10, 30)) +
    geom_abline(
      intercept = thresholdDr,
      slope = 0,
      linetype = 2,
      size = 1
    ) +
    labs(x = "QTL position (Mb)",
         y = expression(bold("QTL significance")),
         parse = TRUE) +
    scale_x_continuous(breaks = c(0, 10, 20) * 10^6,
                       labels = c(0, 10, 20)) + ylim(0, 5.5) + 
    ggtitle("Figure 10: Dietary restriction") 
  
  # Plot peak.ngm.pQTL
  plotLifespanNgm <- ggplot(peak.ngm.QTL, aes(x = qtl_bp, y = qtl_significance, alpha = 0.2)) +
    geom_line(size = 1.5, colour = brewer.pal(9, "Set1")[4]) +
    facet_grid( ~ qtl_chromosome, scales = "free", space = "free_x") +
    presentation +
    theme(legend.position = "none",
          plot.margin = margin(10, 30, 10, 30)) +
    geom_abline(
      intercept = thresholdNgm,
      slope = 0,
      linetype = 2,
      size = 1
    ) +
    labs(x = "QTL position (Mb)",
         y = expression(bold("QTL significance")),
         parse = TRUE) +
    scale_x_continuous(breaks = c(0, 10, 20) * 10^6,
                       labels = c(0, 10, 20)) + ylim(0, 5.5) + 
    ggtitle("Figure 11: NGM") 
  
  # Plot peak.pl.pQTL
  plotLifespanPl <- ggplot(peak.pl.QTL, aes(x = qtl_bp, y = qtl_significance, alpha = 0.2)) +
    geom_line(size = 1.5, colour = brewer.pal(9, "Set1")[4]) +
    facet_grid( ~ qtl_chromosome, scales = "free", space = "free_x") +
    presentation +
    theme(legend.position = "none",
          plot.margin = margin(10, 30, 10, 30)) +
    geom_abline(
      intercept = thresholdPl,
      slope = 0,
      linetype = 2,
      size = 1
    ) +
    labs(x = "QTL position (Mb)",
         y = expression(bold("QTL significance")),
         parse = TRUE) +
    scale_x_continuous(breaks = c(0, 10, 20) * 10^6,
                       labels = c(0, 10, 20)) + ylim(0, 5.5) + 
    ggtitle("Figure 12: Plasticity") 
  
  globalTitle <- grid::textGrob(
    label = paste("QTL mapping for trait", tr),
    gp = grid::gpar(fontsize = 22, fontface = "bold"),
    hjust = 0.5
  )
  
  # Draw to pdf
  left_col <- arrangeGrob(
    plotLifespanNgm,
    plotLifespanDr,
    ncol = 1
  )
  
  right_col <- arrangeGrob(
    plotLifespanPl, 
    blank.plot,
    ncol = 1
  )
  
  # Combine everything
  grid.arrange(
    globalTitle,
    arrangeGrob(left_col, right_col, ncol = 2),
    ncol = 1,
    heights = c(1, 10)
  )
}

# Close pdf
dev.off()







################################################################################
################################################################################
# Lifespan code above should be replaced with this:
################################################################################
################################################################################

popmap_in <- populationMap
popmap_in[is.na(popmap_in)] <- 0

# eQTL mapping
data.QTL.lifespan <- life_data_stats %>%
  dplyr::filter(!Strain %in% c("CB4856", "SCH4856", "N2", "NL5901")) %>%
  dplyr::mutate(treatment_trait = paste0(Treatment, "_", trait)) %>%
  dplyr::select(Strain, treatment_trait, value) %>%
  tidyr::pivot_wider(names_from = Strain, values_from = value) %>%
  tibble::column_to_rownames("treatment_trait")

data.QTL.lifespan <- data.matrix(data.QTL.lifespan)
# CHECK V




# Data preparation for QTL mapping
# Returns a list with the entries Trait (log2 intensities), Map, and Marker, 
# where the Traits are aligned with the map
# AGIWURxxxx = SpotIDs
data.eQTL <- QTL.data.prep(data.eQTL, 
                           colnames(data.eQTL), 
                           populationMap, 
                           populationMarkers)
lapply(data.eQTL, head)  

# function for single marker mapping
# Generate a list with names: LOD, Effect, Trait, Map, and Marker.
aS.eQTL <- QTL.map.1(data.eQTL[[1]], data.eQTL[[2]], data.eQTL[[3]])
save(aS.eQTL, file = file.path(root, paths$output$eqtl, "/obj_aS.eQTL.Rdata"))
str(aS.eQTL)

# Build eQTL list file with calculated threshold value 4.3
# These lines are not executable on workstations due to memory limits
executeLongMethod = FALSE
if(executeLongMethod){
  peak.aS.eQTL <- QTL.map1.dataframe(map1.output = aS.eQTL) %>%
    QTL.peak.finder(threshold = 4.3)
  save(peak.aS.eQTL, file = file.path(root, paths$output$eqtl, "obj_peak.aS.eQTL.Rdata"))
  load(file.path(root, paths$output$eqtl, "obj_peak.aS.eQTL.Rdata"))
} else {
  # Instead of generating: load previously generated peak.aS.eQTL file here:
  load(file.path(root, paths$data$qtl, "obj_peak.aS.eQTL.Rdata"))
}

# Create eQTL table
aS.eQTL.table <- QTL.eQTL.table(
  QTL.peak.dataframe = peak.aS.eQTL,
  cis.trans.window = "LOD.drop",
  trait.annotation = cbind(trait = agi.id[,1],
                           dplyr::select(agi.id,
                                         chromosome, 
                                         gene_bp_start, 
                                         gene_WBID, 
                                         gene_sequence_name, 
                                         gene_public_name))) %>%
  QTL.table.addR2(aS.eQTL) %>%
  group_by(gene_public_name,qtl_chromosome) %>%
  mutate(qtl_representative=ifelse(qtl_R2_sm == max(qtl_R2_sm),"yes","no")) %>%
  data.frame() %>%
  filter(!is.na(gene_bp), qtl_representative == "yes") %>%
  arrange(trait); head(aS.eQTL.table)

# Add transband
call.transbands.file <- QTL.eQTL.call.TBs(aS.eQTL.table, 
                                          window = 0.5e6,
                                          chromosomes = "III",
                                          chromosome_size = 13783801)

# This function finds the number of successes that corresponds to a certain percentile 
# based on an average rate of success. 
# (The corresponding Poisson quantiles for a set of probabilities are obtained.)
qpois(0.0001, 18.36, lower.tail = F)

# For every QTL on chromosome III, indicate whether it is located in a trans-band or cis-enriched area
aS.eQTL.table <-  QTL.eQTL.table.addTB(aS.eQTL.table,
                                       call.transbands.file,
                                       merge_condition = 1)

# Save
save(aS.eQTL.table, file = file.path(root, paths$output$eqtl, "obj_aS.eQTL.table.Rdata"))         


# ------------------------------------------------------------------------------
# power analysis (R only step)
# ------------------------------------------------------------------------------

# Conduct mapping simulation to investigate power
# +/- 35 min
executeLongMethod = FALSE
if(executeLongMethod){
  aS.simulation <- QTL.map.1.sim(strain.map = data.eQTL$Map,
                                 strain.marker = data.eQTL$Marker,
                                 threshold = 4.3,
                                 n_per_marker = 10)
  
  aS.simulation
  save(aS.simulation, file = file.path(root, paths$output$eqtl, "aS.simulation.RData"))
  load(file.path(root, paths$output$eqtl, "aS.simulation.RData"))
} else {
  # Instead of generating: load previously generated aS.simulation file here:
  load(file.path(root, paths$data$qtl, "aS.simulation.RData"))
}

writexl::write_xlsx(aS.eQTL.table,
                    path = file.path(root, paths$output$eqtl, "Supplementary_table5-eQTL.xlsx"))

# Make a table with eQTL peaks for plotting.
# This function adds points to faithfully indicate the chromosome boundaries (standard is C. elegans))
plot.data <- prep.ggplot.eQTL.table(aS.eQTL.table) 

figure_6_locations <- ggplot(plot.data, aes(x = qtl_bp,y = gene_bp)) +
  geom_segment(aes(x = qtl_bp_left,
                   y = gene_bp,
                   xend = qtl_bp_right,
                   yend = gene_bp),
               alpha = 0.25,colour = "darkgrey") + 
  geom_point(aes(colour  =  qtl_type),
             size = 0.6, alpha = 0.5)  +
  facet_grid(gene_chromosome ~ qtl_chromosome, 
             space  =  "free",
             scales = "free") + 
  presentation + colScale + 
  theme(legend.position  =  "none") +
  labs(x = "eQTL peak position (Mb)",
       y = "Gene position (Mb)") + 
  ggtitle("Figure 6: eQTL locations") +
  theme(panel.spacing = unit(0.1,"lines")) +
  scale_x_continuous(breaks = c(5,10,15,20)*10^6, labels = c(5,10,15,20)) +
  scale_y_continuous(breaks = c(5,10,15,20)*10^6, labels = c(5,10,15,20))

figure_6_locations

plot.data <- filter(aS.eQTL.table, qtl_type == "trans") %>%
  prep.ggplot.eQTL.table() 

figure_7_eQTL_counts <- ggplot(plot.data, aes(x = qtl_bp,fill = qtl_type)) +
  geom_histogram(binwidth = 500000) + geom_hline(yintercept = 36,lty = 2,lwd = 1.2) +
  facet_grid(c("")~qtl_chromosome, space  =  "free",scales = "free") + presentation + 
  fillScale + theme(legend.position  =  "none",
                    panel.spacing = unit(0.1,"lines")) +
  labs(x = "eQTL peak position (Mb)", y = "eQTL counts") +
  ggtitle("Figure 7: eQTL") +
  scale_x_continuous(breaks = c(5,10,15,20)*10^6,
                     labels = c(5,10,15,20))

figure_7_eQTL_counts
################################################################################
################################################################################
# Lifespan code above should be replaced with this
################################################################################
################################################################################


