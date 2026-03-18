################################################################################################################################################################
### genes.check(trans.int,spot.id,genes.to.check,filename) v1.0
###   returns a plot of genes that should be checked (saved as pdf)
###   input
###     spot.id, vector of identifiers of reference genes (gene names), only if genes.check=TRUE
###     genes.to.check, identifiers that are also present in the spot.id vector



transcriptomics.check.genes <- 
  function(trans.int,spot.id,genes.to.check,filename,save_dir){
                if(missing(spot.id)){                                   
                  stop("Require spot IDs (gene names)")}
                if(missing(genes.to.check)){                            
                  genes.to.check <- c("gst-27","vit-1","vit-2","hsp-16.2","hsp-16.41")}
                if(missing(filename)){                                  
                  filename <- ""}
                if(missing(save_dir)){                                  
                  save_dir <- getwd()}

                checkgenes <- as.data.frame(cbind(id=as.character(spot.id[spot.id %in% genes.to.check]),round(trans.int$log2_intensities[spot.id %in% genes.to.check,],digits=5))) %>%
                              gather(Sample,Expression,-1)
                checkgenes$Expression <- as.numeric(checkgenes$Expression)

                    ggplot(checkgenes) + aes(x=Sample,y=Expression) +
                    facet_grid(id ~.) + geom_jitter(position = position_jitter(width = .25)) +
                    geom_boxplot(outlier.size = 0, aes(alpha=0.5)) + presentation +
                    theme(axis.text.x = element_text(angle = 90,vjust = 0.5))

                if(filename==""){ggsave(filename=paste(save_dir,date.now(),"_genes.check.pdf",sep=""),width=3+length(unique(checkgenes$Sample))*0.25,height=15,limitsize=FALSE)}
                if(filename!=""){ggsave(filename=paste(save_dir,date.now(),"_",filename,"_genes.check.pdf",sep=""),width=3+length(unique(checkgenes$Sample))*0.25,height=15,limitsize=FALSE)}
               }
