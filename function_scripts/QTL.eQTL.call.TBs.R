###Calls transbands based on a poisson distribution
#Brem, Rachel B., and Leonid Kruglyak. "The landscape of genetic complexity across 5,700 gene expression traits in yeast." Proceedings of the National Academy of Sciences 102.5 (2005): 1572-1577.

###input
#eQTL.table.file (output from QTL.eQTL.table)
#window (interval size in which QTL are counted)
#chromosomes (names of the chromosomes; standard is C. elegans)
#chromosome_size (sizes of the chromosomes; standard is C. elegans)

###output
#a dataframe with sifnificances for cis- and trans-eQTL enrichments

###Description



QTL.eQTL.call.TBs <- function(eQTL.table.file,window,chromosomes,chromosome_size){
                               if(missing(eQTL.table.file)){                        stop("need a eQTL.table.file")}
                               if(missing(window)){                                 window <- 1.0e6; warning("set default window-size, might not be optimal")}
                               if(missing(chromosomes)){                            chromosomes <- c("I","II","III","IV","V","X")}
                               if(missing(chromosome_size)){                        chromosome_size <- c(15072434,15279421,13783801,17493829,20924180,17718942)}
                               if(length(chromosomes) != length(chromosome_size)){  stop("chromosome names and sizes should have the same length")}

                               maxsize <- max(chromosome_size+window)  
                               chrnum <- length(chromosomes)
    
                               ###call locations with many eQTL
                               transband.id <- dplyr::mutate(eQTL.table.file,Interval=findInterval(qtl_bp,seq(1,maxsize,by=window))) %>%
                                               dplyr::filter(!is.na(qtl_type)) %>%
                                               dplyr::group_by(qtl_chromosome,Interval,qtl_type) %>%
                                               dplyr::summarise(n.ct=length(unique(trait))) %>%
                                               data.frame() %>%
                                               dplyr::group_by(qtl_type) %>%
                                               dplyr::mutate(exp.ct=mean(as.numeric(unlist(n.ct)))) %>%
                                               data.frame() %>%
                                               dplyr::mutate(transband_significance=ppois(n.ct,lambda=exp.ct,lower.tail=F))
                                
                               ###Add locations that have no cis/trans eQTL
                               transband.id <- rbind(transband.id,data.frame(cbind(qtl_chromosome=rep(chromosomes,each=ceiling(maxsize/window)*2),Interval=rep(1:ceiling(maxsize/window),times=chrnum*2),
                                                                                   qtl_type=rep(rep(c("cis","trans"),each=ceiling(maxsize/window)),times=chrnum),
                                                                                   n.ct=0,
                                                                                   exp.ct=rep(rep(c(filter(transband.id,!duplicated(exp.ct),qtl_type=="cis")$exp.ct,filter(transband.id,!duplicated(exp.ct),qtl_type=="trans")$exp.ct),each=ceiling(maxsize/window)),times=chrnum),
                                                                                   transband_significance=1))) %>%
                                               dplyr::mutate(tester=paste(qtl_chromosome,Interval,qtl_type)) %>%
                                               dplyr::filter(!duplicated(tester)) %>%
                                               dplyr::mutate(Interval=as.numeric(as.character(unlist(Interval))),transband_significance=as.numeric(transband_significance)) %>%
                                               dplyr::mutate(Interval_left=(Interval-1)*window,Interval_right=1+Interval*window) %>%
                                               dplyr::arrange(qtl_type,qtl_chromosome,Interval)
                               
                               ###resize to chromosome sizes
                               output <- NULL
                               for(i in 1:chrnum){
                                   output <- rbind(output,filter(transband.id,qtl_chromosome==chromosomes[i],Interval_right<(chromosome_size[i]+window)))
                               }
                               for(i in c(2,4,5,6,8,9)){
                                   output[,i] <- as.numeric(as.character(unlist(output[,i])))
                               }
                               
                               return(output)
                              }