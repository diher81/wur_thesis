    ###Makes a table with eQTL peaks for plotting
    ###this function adds points to faithfully indicate the chromosome boundaries (standard is C. elegans))
    prep.ggplot.eQTL.table <- function(eQTL.table.file,condition,chromosomes,chromosome_size){
                                       if(missing(eQTL.table.file)){                        stop("need a eQTL.table.file")}
                                       if(missing(condition)){                              condition <- "ct"}
                                       if(missing(chromosomes)){                            chromosomes <- c("I","II","III","IV","V","X")}
                                       if(missing(chromosome_size)){                        chromosome_size <- c(15072434,15279421,13783801,17493829,20924180,17718942)}
                                       if(length(chromosomes) != length(chromosome_size)){  stop("chromosome names and sizes should have the same length")}

                                       output <- dplyr::filter(eQTL.table.file, gene_chromosome %in% chromosomes) %>%
                                                 dplyr::mutate(gene_chromosome = factor(gene_chromosome,levels=rev(chromosomes))) %>%
                                                 dplyr::mutate(qtl_chromosome = factor(qtl_chromosome,levels=chromosomes),condition=condition) %>%
                                                 dplyr::mutate(qtl_type = factor(qtl_type,levels=rev(sort(unique(qtl_type)))))

                                       ###Chromosome sizes cheat
                                       cheatpoints <- output[rep(1,times=length(chromosomes)*2),]
                                       for(i in which(unlist(lapply(output,class)) == "factor")){
                                           cheatpoints[,i] <- as.character(unlist(cheatpoints[,i]))
                                       }
                                       cheatpoints[1:nrow(cheatpoints),unlist(lapply(output,class)) != "character"] <- 0
                                       cheatpoints[1:nrow(cheatpoints),unlist(lapply(output,class)) != "numeric"] <- ""
                                       cheatpoints$qtl_chromosome <- as.character(rep(chromosomes,times=2))
                                       cheatpoints$gene_chromosome <- as.character(rep(chromosomes,times=2))
                                       cheatpoints$qtl_bp <- as.numeric(as.character(c(chromosome_size,rep(1,times=length(chromosomes)))))
                                       cheatpoints$gene_bp <- as.numeric(as.character(c(chromosome_size,rep(1,times=length(chromosomes)))))
                                       cheatpoints$trait <- paste(as.character(rep(chromosomes,times=2)),rep(c(1,2),each=length(chromosomes)))


                                       ###combine
                                       output <- rbind(output,cheatpoints) %>%
                                                 dplyr::mutate(qtl_type = factor(qtl_type,levels = c("cis","trans","")))
                                       
                                       ###return
                                       return(output)
                                      }

