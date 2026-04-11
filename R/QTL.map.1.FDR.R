###FDR calculation function for single marker mapping
#Snoek & Sterken, 2017; Sterken, 2017

###input
#map1.output (output from the QTL.map.1 function)
#filenames.perm (files containing permutated data)
#FDR_dir (location of the files with permutated data)
#q.value (desired q-value)
#small (logic; set to TRUE for a small set)

###output
#a list with names: LOD, Effect, Trait_perm, Map, and Marker.

###Description
# Based on Benjamini and Yakuteili multiple testing under dependency



QTL.map.1.FDR <- function(map1.output,filenames.perm,FDR_dir,q.value,small){
                         if(missing(map1.output)){        stop("Missing map1.output ([[LOD]],[[Effect]],[[Trait]],[[Map]],[[Marker]])")}
                         if(missing(filenames.perm)){     stop("Need a vector with filenames for permutation")}
                         if(missing(FDR_dir)){            FDR_dir <- getwd()}
                         if(missing(q.value)){            q.value <- 0.025}
                         if(missing(small)){              small <- FALSE}

                         output <- as.list(NULL)

                         output[[2]] <- matrix(0,ncol=length(filenames.perm),nrow=5000)
                         rownames(output[[2]]) <- seq(0.1,500,by=0.1)
                         colnames(output[[2]]) <- filenames.perm

                         ###in case of one trait
                         if(small){
                             if(length(filenames.perm) > 1){
                                 stop("needs one permutatation file, with n permutations")
                             }
                             tmp <- load(file=paste(FDR_dir,filenames.perm,sep=""))
                             ###No NA's
                             tmp <- get(tmp)[[1]]
                             tmp[is.na(tmp)] <- 0.1
                             
                             ###no Non-traits
                             tmp <- tmp[rownames(tmp) != "",]
                             output[[2]] <- apply(tmp,1,max,na.rm=T)
                             
                             if((length(output[[2]]) * q.value) < 10){
                                  warning(paste("inadvisable to calculate q=",q.value,"based on ",length(output[[2]])," permutations. Increase to at least ", ceiling(10/q.value)))
                             }
                             
                             ###take the q-cutoff, *10 is to be consistent with eQTL methodology (is per 0.1)
                             output[[1]] <- rev(sort(output[[2]]))[ceiling(q.value*length(output[[2]]))]*10
                             
                             ###again, to be consistent; 
                             output[[3]] <- mean(output[[2]],na.rm = TRUE)
                             
                             ###real data
                             output[[4]] <- max(map1.output$LOD,na.rm=TRUE)
                             
                             names(output) <- c("Significance_threshold","False_discoveries","False_discoveries_average","Real_discoveries")
                         }
                         if(!small){
                             for(i in 1:length(filenames.perm)){
                                 tmp <- load(file=paste(FDR_dir,filenames.perm[i],sep=""))
                                 ###No NA's
                                 tmp <- get(tmp)[[1]]
                                 tmp[is.na(tmp)] <- 0.1

                                 tmp <- table(apply(round(tmp,digits=1),1,max,na.rm=T))

                                 output[[2]][(as.numeric(names(tmp))*10),i] <- tmp
                                 output[[2]][,i] <- rev(cumsum(rev(output[[2]][,i])))
                             }
                             ###Take the average over permutations
                             output[[3]] <- apply(output[[2]],1,mean)

                             ###Calculate the real discoveries, no NA's
                             tmp <- map1.output$LOD
                             tmp[is.na(tmp)] <- 0.1

                             RDS.tmp <- table(round(apply(tmp,1,max,na.rm=T),digits=1))
                             RDS <- rep(0,times=5000); RDS[(as.numeric(names(RDS.tmp))*10)] <- RDS.tmp
                             RDS <- rev(cumsum(rev(RDS)))

                             output[[4]] <- RDS

                             output[[1]] <- which(round(((dim(map1.output$LOD)[1]-RDS)/dim(map1.output$LOD)[1])*q.value*log10(dim(map1.output$LOD)[1]),digits=3)-round(output[[3]]/output[[4]],digits=3)>0)[1]
                             names(output) <- c("Significance_threshold","False_discoveries","False_discoveries_average","Real_discoveries")
                         }
                         
                         ###Return
                         return(output)
                        }
