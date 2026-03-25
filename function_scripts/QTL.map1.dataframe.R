###function for transforming output of QTL.map.1 to a dataframe

###input
#QTL.map.1 output
#small (logic) if a single trait was mapped set to TRUE

###output
#dataframe with the columns trait, qtl_chromosome, qtl_bp, qtl_marker, qtl_significance, and qtl_effect

###Description
#This transforms the relevant output to a dataframe, making handling in dplyr and plotting using ggplot easier
#This takes output of map.1 functions (list file with LOD, effect, Trait, and Map and converts it to a list.
#Optimized (a bit) by generating an empty matrix and filling the columns one-by-one.
#added an option for single traits

    QTL.map1.dataframe <- function(map1.output,small){
                                if(missing(map1.output)){        stop("Missing map1.output ([[LOD]],[[Effect]],[[Trait]],[[Map]],[[Marker]])")}
                                if(missing(small)){              small <- FALSE}
      
                                if(is.null(dim(map1.output[[1]]))){
                                    small <- TRUE
                                }
 
                                if(small){
                                    ###Split up to long format
                                    output <- matrix(NA,length(map1.output[[1]]),6)
                                    output[,1] <- rep(rownames(map1.output[[3]])[1],each=length(map1.output[[1]]))
                                    output[,2] <- as.character(unlist(map1.output[[5]][,2]))
                                    output[,3] <- as.numeric(as.character(unlist(map1.output[[5]][,3])))
                                    output[,4] <- c(1:length(map1.output[[1]]))
                                    output[,5] <- c(map1.output[[1]])
                                    output[,6] <- c(map1.output[[2]])
                                    output <- as.data.frame(output)
                                }     
                                if(!small){
                                    ###Split up to long format
                                    output <- matrix(NA,ncol(map1.output[[1]])*nrow(map1.output[[1]]),6)
                                    output[,1] <- rep(rownames(map1.output[[1]]),each=ncol(map1.output[[1]]))
                                    output[,2] <- rep(as.character(unlist(map1.output[[5]][,2])),times=nrow(map1.output[[1]]))
                                    output[,3] <- rep(as.numeric(as.character(unlist(map1.output[[5]][,3]))),times=nrow(map1.output[[1]]))
                                    output[,4] <- rep(c(1:ncol(map1.output[[1]])),times=nrow(map1.output[[1]]))
                                    output[,5] <- c(t(map1.output[[1]]))
                                    output[,6] <- c(t(map1.output[[2]]))
                                    output <- as.data.frame(output)
                                }

                                ###overwrite dataframe formats
                                output[,1] <- as.character(unlist(output[,1]))
                                output[,2] <- as.character(unlist(output[,2]))
                                output[,3] <- as.numeric(as.character(unlist(output[,3])))
                                output[,4] <- as.numeric(as.character(unlist(output[,4])))
                                output[,5] <- as.numeric(as.character(unlist(output[,5])))
                                output[,6] <- as.numeric(as.character(unlist(output[,6])))
                                colnames(output) <- c("trait","qtl_chromosome","qtl_bp","qtl_marker","qtl_significance","qtl_effect")

                                ###Return
                                return(output)
                               }
