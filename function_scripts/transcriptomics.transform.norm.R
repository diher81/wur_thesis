################################################################################################################################################################
### transform.norm(rg.norm) v1.0
###   returns a list with "log2.rat.mean","log2.int","standard_score"
###   obligatory input
###     rg.norm, output of the limma.norm function, or other any other non-transformed data matrix
###   optional input
###     filename, output file name, default is "obj_trans.int.[date].out"


transcriptomics.transform.norm <- function(rg.norm,filename,save_dir){
                                           if(missing(rg.norm)){                                stop("Require output of normalization function")}
                                           if(missing(filename)){                               filename <- ""}
                                           if(missing(save_dir)){                               save_dir <- getwd()}
                
                                           trans.int <- NULL; trans.int <- as.list(trans.int)
                                           trans.int[[1]] <- log2(rg.norm/apply(rg.norm,1,mean))
                                           trans.int[[2]] <- log2(rg.norm)
                                           trans.int[[3]] <- (rg.norm-apply(rg.norm,1,mean))/apply(rg.norm,1,sd)
                                           names(trans.int) <- c("log2_ratio_mean","log2_intensities","z_score")
                
                                           if(filename==""){save(trans.int,file=paste(save_dir,"/obj_trans.norm.",date.now(),".Rdata",sep=""))}
                                           if(filename!=""){save(trans.int,file=paste(save_dir,"/obj_trans.norm.",filename,".",date.now(),".Rdata",sep=""))}
                
                                           return(trans.int)
                                          }
