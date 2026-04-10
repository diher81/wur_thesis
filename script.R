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
    "tidyverse",
    "RColorBrewer",
    "grid",
    "gridExtra",
    "httr",
    "jsonlite"
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
dirHome <- "/Users/diher/Repos/wur/thesis_dirk/" 
dirData <- paste0(dirHome, "data/")
dirDataMA <- paste0(dirHome, "data/MA/")
dirFunctions <- paste0(dirHome, "function_scripts/")
dirTarget <- paste0(dirHome, "target/")
dirOutput <- paste0(dirHome, "output/")
dirOutputFdr <- paste0(dirHome, "output/FDR/")
dirNormalized <- paste0(dirHome, "normalized/")

setwd(dirHome)


# ------------------------------------------------------------------------------
# Loading data
# ------------------------------------------------------------------------------

# Targets
targets_RIL <- read.delim(paste0(dirTarget, "Targets_RIL.txt"), 
                          stringsAsFactors = FALSE)
agi.id <- read.delim(paste0(dirTarget, "ArrayID_agilentV2_WS258.txt"), 
                     stringsAsFactors = FALSE)

# load DNAseq map
populationMap <- data.matrix(read.table(paste0(dirData, 
  "Genetic_map/asRIL_map_new.txt"))[,-c(1:3,5,6,8,9,11,13)])
populationMarkers <- read.table(paste0(dirData, 
  "Genetic_map/asRIL_map_new.txt"))[,c(1:3)]

# eQTL
load(file = paste0(dirOutput, "obj_peak.aS.eQTL.Rdata"))

# Protein accumulation
load(file = paste0(dirData, "proteinAccumulation/obj_elisa_aS.Rdata"))
load(file = paste0(dirData, "proteinAccumulation/obj_qpcr_aS.Rdata"))

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
# Functions
# ------------------------------------------------------------------------------

setwd(dirFunctions)
for(i in 1:length(dir())){
  source(dir()[i])
}
setwd(dirHome)


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
rg.normalized.intensities <- normalize.agilent.transcriptomics(targets = targets_RIL,
                save_dir = dirNormalized)

# Generate a list with "log2.rat.mean", "log2.int", "standard_score"
transformed.intensities <- transcriptomics.transform.norm(rg.normalized.intensities,
                save_dir = dirDataMA)

# Checks
correlsums <- transcriptomics.check.cor(transformed.intensities, save_dir = dirOutput)
transcriptomics.check.genes(transformed.intensities,
                            spot.id = agi.id$gene_public_name,
                            save_dir = dirDataMA)

# Make and save a list of log2 transformed intensities and log2 ratio means
colnames.names <- c("number","strain","batch","alphasyn","days","sample_number")
list.data.RIL <- transcriptomics.list.to.df(trans.int = transformed.intensities, 
                                        spot.id = agi.id$SpotID,
                                        colnames.sep = ":", 
                                        colnames.names = colnames.names)
save(list.data.RIL, 
     file = paste0(dirNormalized, "/obj_list.data.RIL.Rdata"))


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
  geom_rect() + scale_fill_manual(values = c("#0080FF","#FF8000",NA)) +
  facet_grid(.~chrom, scales="free", space="free") +
  presentation + xlab("Chromosome position (Mb)") + ylab("Strain") +
  ggtitle("Figure 1: Genetic map") +
  scale_x_continuous(labels = function(x) { x/1e6 }, expand = c(0,0)) +
  scale_y_continuous(breaks = unique(data.plot$index) + 0.5, labels = function(x) { strain_index[as.character(x)] }, expand = c(0,0)) +
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

annotation.grobA <- title.grob <- textGrob(label = "A",
                                           x = unit(0, "lines"),
                                           y = unit(0, "lines"),
                                           hjust = 0, 
                                           vjust = 0,
                                           gp = gpar(fontsize = 20,fontface="bold"))
annotation.grobB <- title.grob <- textGrob(label = "B",
                                           x = unit(0, "lines"),
                                           y = unit(0, "lines"),
                                           hjust = 0, 
                                           vjust = 0,
                                           gp = gpar(fontsize = 20,fontface="bold"))
annotation.grobC <- title.grob <- textGrob(label = "C",
                                           x = unit(0, "lines"),
                                           y = unit(0, "lines"),
                                           hjust = 0, 
                                           vjust = 0,
                                           gp = gpar(fontsize = 20,fontface="bold"))

graph.oject.gen.map <- arrangeGrob(figure_1_gen_map, top = annotation.grobA)
graph.oject.cm.map  <- arrangeGrob(figure_2_cm_map, top = annotation.grobB)
graph.oject.gen.distr  <- arrangeGrob(figure_3_gen_distr, top = annotation.grobC)

# Open pdf device
pdf(file = paste0(dirOutput, "figure1-Genetic_map.pdf"), width = 10, height = 14)
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
# Returns a list with the entries Trait (log2 intensities), Map, and Marker, where the Traits are aligned with the map
# AGIWURxxxx = SpotIDs
data.eQTL <- QTL.data.prep(data.eQTL, 
                           colnames(data.eQTL), 
                           populationMap, 
                           populationMarkers)
lapply(data.eQTL, head)  

# function for single marker mapping
# Generate a list with names: LOD, Effect, Trait, Map, and Marker.
# LOD score: Logarithm Of Odds. A statistical estimate of the relative probability
# of a QTL at this locus.
aS.eQTL <- QTL.map.1(data.eQTL[[1]], data.eQTL[[2]], data.eQTL[[3]])
save(aS.eQTL, file = paste0(dirOutput, "/obj_aS.eQTL.Rdata"))
str(aS.eQTL)

# Build eQTL list file with calculated threshold value 4.3
# These lines are not executable on workstations due to memory limits
# Existing obj_peak.aS.eQTL.Rdata file was loaded at the start of this script
# Genetisch kaart niet onafhankelijk want linkage, Dus ik weet niet hoeveel onafhankelijke testen er zijn
# hoe groot is de kans dat ik deze piek vind op basis van random data (toevallig dus)
executeLongMethod = FALSE
if(executeLongMethod){
  peak.aS.eQTL <- QTL.map1.dataframe(map1.output = aS.eQTL) %>%
    QTL.peak.finder(threshold = 4.3)
  save(peak.aS.eQTL,
       file = paste0(dirOutput, "obj_peak.aS.eQTL.Rdata"))
}

# Create eQTL table
# lod drop is natte vinger werk
# niet enkel de piek, maar de hele regio moet boven de top van de piek min 1,5 liggen
# breedte van die regio is afhank van linkage en grootte van populatie
# veel recombinaties belangrijk voor het opbreken van de linkage
aS.eQTL.table <- QTL.eQTL.table(QTL.peak.dataframe = peak.aS.eQTL, 
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

# Adds whether a QTL is located in a trans-band (or cis-enriched area)
aS.eQTL.table <-  QTL.eQTL.table.addTB(aS.eQTL.table,
                                       call.transbands.file,
                                       merge_condition = 1)

# Save
save(aS.eQTL.table, file = paste0(dirOutput, "obj_aS.eQTL.table.Rdata"))         


# BELANGRIJK Voor transband. Niet voor mij. Mag dus verwijderd worden.

###Stats for text
###genes in transband
sum(table(aS.eQTL.table$trans_band)) - table(aS.eQTL.table$trans_band)[names(table(aS.eQTL.table$trans_band)) == "none"]

###number of transband
length(table(aS.eQTL.table$trans_band)) - 1

###cis / trans
table(aS.eQTL.table$qtl_type)

###total
length(unique(aS.eQTL.table$gene_WBID))

###TBs on chromosome III
table(aS.eQTL.table$trans_band)


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
  save(aS.simulation, file = file.path(dirOutput, "aS.simulation.RData"))
} else {
  # Instead of generating: load previously generated aS.simulation here:
  load(file.path(dirOutput, "aS.simulation.RData"))
}

writexl::write_xlsx(aS.eQTL.table,
                    path = paste0(dirOutput, "Supplementary_table5-eQTL.xlsx"))

# Make a table with eQTL peaks for plotting.
# This function adds points to faithfully indicate the chromosome boundaries (standard is C. elegans))
plot.data <- prep.ggplot.eQTL.table(aS.eQTL.table) 

figure_6_positions <- ggplot(plot.data, aes(x = qtl_bp,y = gene_bp)) +
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
  ggtitle("Figure 6: positions") +
  theme(panel.spacing = unit(0.1,"lines")) +
  scale_x_continuous(breaks = c(5,10,15,20)*10^6, labels = c(5,10,15,20)) +
  scale_y_continuous(breaks = c(5,10,15,20)*10^6, labels = c(5,10,15,20))

figure_6_positions

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


# ------------------------------------------------------------------------------
# Draw eQTL profiles and boxplots for genes of interest
# ------------------------------------------------------------------------------

# Gene names and spotIds of the selected genes
spotIds <- aS.eQTL.table[[1]]
geneNames <- aS.eQTL.table[[16]]

# Open pdf device
pdf(file = paste0(dirOutput, "boxplotsForGenes.pdf"), width = 10, height = 14)

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
    geom_jitter(height=0,
                width=0.25,
                aes(colour=geno_strain),
                alpha=0.2) + 
    geom_boxplot(outlier.shape=NA, alpha=0.2, aes(fill=geno_strain)) +
    labs(x = "Genotype at marker", 
         y = paste(geneNames[i], " expression"),
         parse=TRUE) + 
    facet_grid(~Chromosome+Position) +
    presentation + 
    colScale + 
    fillScale + 
    theme(legend.position = "none",
          plot.margin = margin(10, 30, 10, 30)) +
    annotate("text",
             x=1.5,
             y=max(data.plot[[2]]$trait_value,na.rm=T),
             label=paste("italic(R)^{2}==",
                         round(data.plot[[2]]$R_squared[1],
                               digits=2),
                         sep=""),
             parse=TRUE)
  
  globalTitle <- textGrob(
    label = paste("Gene", geneNames[i], ": eQTL profile and genotype split-out."),
    gp = gpar(fontsize = 22, fontface = "bold"),
    hjust = 0.5
  )

  # Draw to pdf
  grid.arrange(globalTitle,
               plotEqtlProfile, 
               nullGrob(),
               plotBoxplot,
               nullGrob(),
               ncol = 1,
               heights = c(1, 5, 1.5, 5, 1))
}

# Close pdf
dev.off()


# ------------------------------------------------------------------------------
# Protein accumulation
# ------------------------------------------------------------------------------

# eQTL mapping
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
# TODO I get this error: 
# Error in cor.test.default(data.input[i, !is.na(variable)], variable[!is.na(variable)],  : not enough finite observations
aS.pQTL <- QTL.map.1(data.pQTL[[1]], data.pQTL[[2]], data.pQTL[[3]])
save(aS.pQTL, file = paste0(dirOutput, "/obj_aS.pQTL.Rdata"))

# Calculate threshold
# TODO same error here
# Error in cor.test.default(data.input[i, !is.na(variable)], variable[!is.na(variable)],  : not enough finite observations

# berekening om de theshold te kunnen berekenen
aS.pQTL.perm <- QTL.map.1.perm(data.pQTL[[1]], data.pQTL[[2]], data.pQTL[[3]], 1000) # hier 1000 permutaties
save(aS.pQTL.perm,
     file = paste0(dirOutputFdr, "obj_aS.pQTL.perm", i, ".Rdata"))

# hier de berekende value invullen die ik opgeslagen had. Zie mail Mark (ROND DE 3)
# deze berekent de threshold
aS.FDR <- QTL.map.1.FDR(map1.output = aS.pQTL,
                        filenames.perm = filenames.perm, 
                        FDR_dir = dirOutputFdr,
                        q.value =0.025); aS.FDR[[1]]
save(aS.FDR, file = paste0(dirOutputFdr, "obj_RIL.aS.FDR.Rdata"))

# Build QTL list file with calculated threshold value 4.3 
# TODO QUESTION: Update in order to use calculated threshold value
peak.aS.pQTL <- QTL.map1.dataframe(map1.output = aS.pQTL) %>%
  QTL.peak.finder(threshold = 4.3)
save(peak.aS.pQTL, file = paste0(dirOutput, "obj_peak.aS.QTL.Rdata"))

# Klaar. Figuur plotten van piek peak.aS.pQTL
# call peak en plot

# TODO QUESTION: Should I be using this "...eQTL..." method?
# niet nodig
aS.pQTL.table <- QTL.eQTL.table(QTL.peak.dataframe = peak.aS.pQTL, 
                                cis.trans.window = "LOD.drop",
                                trait.annotation = cbind(trait = agi.id[,1],
                                                         dplyr::select(agi.id,
                                                                       chromosome, 
                                                                       gene_bp_start, 
                                                                       gene_WBID, 
                                                                       gene_sequence_name, 
                                                                       gene_public_name))) %>%
  QTL.table.addR2(aS.pQTL) %>%
  group_by(gene_public_name, qtl_chromosome) %>%
  mutate(qtl_representative = ifelse(qtl_R2_sm == max(qtl_R2_sm),"yes","no")) %>%
  data.frame() %>%
  filter(!is.na(gene_bp), qtl_representative == "yes") %>%
  arrange(trait); head(aS.pQTL.table)

# Add transband
# niet meer
call.transbands.file <- QTL.eQTL.call.TBs(aS.pQTL.table, 
                                          window = 0.5e6,
                                          chromosomes = "III",
                                          chromosome_size = 13783801)











