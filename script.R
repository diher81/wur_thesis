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

# load DNAseq map
population.map <- data.matrix(read.table(
  "./Data/Genetic_map/asRIL_map_new.txt")[,-c(1:3,5,6,8,9,11,13)])
population.markers <- read.table(
  "./Data/Genetic_map/asRIL_map_new.txt")[,c(1:3)]
IL_gen <- read.delim("./Data/Genetic_map/asIL_map_new.txt") 
str(IL_gen) # print its structure
IL.map <- load(file = "./Data/Genetic_map/obj_IL_map.Rdata")
str(get(IL.map)) # print its structure


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


# ------------------------------------------------------------------------------
# Genetic map
# ------------------------------------------------------------------------------

# Figure 1: Genetic map --------------------------------------------------------

genetic_map <- population.map;
# replace missing values with 0
genetic_map[is.na(genetic_map)] <- 0
# rewrite the genetic map to a dataframe plottable in ggplot
data.plot <- prep.ggplot.genetic.map(genetic_map, population.markers) %>%
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


# Figure 2: centimorgan map-----------------------------------------------------

strain.map <- popmap[,-c(1:4)]
strain.map[strain.map==0] <- NA

data.plot <- genetic.distance(strain.map,popmrk)

sf1b <- ggplot(data.plot,aes(x=position,y=cM)) +
  geom_line() + geom_rug() + facet_grid(~chromosome,space = "free_x", scales = "free_x") + 
  presentation + labs(x="Physical position (Mb)",y="Genetic position (cM)")+
  ggtitle("Suppl figure 1B centimorgan map") +
  scale_x_continuous(breaks=c(5,10,15,20)*10^6,labels=c(5,10,15,20))
print(sf1b)


# Figure 3: Genotype distribution ----------------------------------------------

data.plot <- tidyr::gather(data.frame(cbind(popmrk,popmap[,-c(1:4)])),key=strain,value=genotype,-Name,-Chromosome,-Position) %>%
  dplyr::group_by(Chromosome,Position,genotype) %>%
  dplyr::summarise(n=length(genotype)) %>%
  dplyr::mutate(genotype=ifelse(genotype==-1,"CB4856",
                                ifelse(genotype==1,"N2","NA"))) %>%
  data.frame() %>%
  dplyr::mutate(genotype=factor(genotype,levels = c("NA","CB4856","N2")))

sf1c <- ggplot(data.plot,aes(x=Position,y=n,fill=genotype)) +
  geom_area() + facet_grid(~Chromosome,space = "free_x",scales="free_x") + fillScale +
  presentation + labs(x="Position (Mb)",y="Genotype (n strains)") +
  ggtitle("Suppl figure 1C Genotype distribution ") +
  scale_x_continuous(breaks=c(5,10,15,20)*10^6,labels=c(5,10,15,20)) + scale_y_continuous(expand=c(0,0)) +
  guides(fill=guide_legend(title="Genotype"))
print(sf1c)


# Figure 4: Genetic map --------------------------------------------------------

BlPl <- ggplot()+geom_blank(aes(1,1))+ blank_theme

lay <- rbind(c(1,2),c(1,3),c(1,4))

annotation.grobA <- title.grob <- textGrob(label = "A",x = unit(0, "lines"),y = unit(0, "lines"),hjust = 0, vjust = 0,gp = gpar(fontsize = 20,fontface="bold"))
annotation.grobB <- title.grob <- textGrob(label = "B",x = unit(0, "lines"),y = unit(0, "lines"),hjust = 0, vjust = 0,gp = gpar(fontsize = 20,fontface="bold"))
annotation.grobC <- title.grob <- textGrob(label = "C",x = unit(0, "lines"),y = unit(0, "lines"),hjust = 0, vjust = 0,gp = gpar(fontsize = 20,fontface="bold"))

sf1a <- arrangeGrob(sf1a,top=annotation.grobA)
sf1a
sf1b <- arrangeGrob(sf1b,top=annotation.grobB)
sf1b
sf1c <- arrangeGrob(sf1c,top=annotation.grobC)
sf1c

pdf(file="Supplementary_figure1-Genetic_map.pdf",width=10,height=14)
grid.arrange(sf1a,sf1b,sf1c,BlPl,layout_matrix=lay,heights=c(3,3,8))
dev.off()


