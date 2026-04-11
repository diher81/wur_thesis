################################################################################################################################################################
### list.to.dataframe(trans.int,spot.id,colnames.sep,colnames.names) v1.1
###   returns a data frame in long format for easy processing with dplyer and ggplot2
###   obligatory input
###     trans.int, output from transform.int function
###     spot.id, for example the agilentID numbers
###   optional input
###     colnames.sep, a regular expression that separates the column descriptors (e.g. : or ; or _)
###     colnames.names, the proposed names of the columns that exist after separation

transcriptomics.list.to.df <- function(trans.int,spot.id,colnames.sep,colnames.names){
                                      if(missing(trans.int)){                           stop("Require output of transform.norm")}
                                      if(missing(spot.id)){                             stop("Require spot IDs (spot numbers)")}
                                      if(missing(colnames.sep) == FALSE){               mut.colnames <- TRUE}
                                      if(missing(colnames.sep)){                        mut.colnames <- FALSE}
                                      if(missing(colnames.sep) == FALSE &
                                         missing(colnames.names)){                      stop("Give the column names for the samples")}
        
                                      ###Depends on tidyr, transform to list
                                      tmp1 <- cbind(spot.id,as.data.frame(trans.int$log2_ratio_mean)) %>%
                                              tidyr::gather(Identifier,log2.rat,-spot.id)
        
                                      tmp2 <- cbind(spot.id,as.data.frame(trans.int$log2_intensities)) %>%
                                              tidyr::gather(Identifier,log2.int,-spot.id)
        
                                      tmp3 <- cbind(spot.id,as.data.frame(trans.int$z_score)) %>%
                                              tidyr::gather(Identifier,z.score,-spot.id)
        
                                      ###Merge together
                                      output <- cbind(tmp1,tmp2$log2.int,tmp3$z.score)
                                          colnames(output) <- c("SpotID","Identifier","log2_ratio_mean","log2_intensities","z_score")
        
                                      ###Depends on tidier, split column names
                                      if(mut.colnames){
                                          output <- tidyr::separate(output,Identifier,into=colnames.names, sep=colnames.sep)
                                          output <- cbind(tmp1$Identifier,output)
                                          colnames(output) <- c("SampleID",colnames(output)[2:length(colnames(output))])
                                      }
                                      return(output)
                                     }