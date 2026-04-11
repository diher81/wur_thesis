###Permutation function for single marker mapping
#Snoek & Sterken, 2017; Sterken, 2017

###input
#trait.matrix (traits in rows, genotypes in columns)
#strain.map (markers in rows, genotypes in columns)
#strain.marker (markers in rows, with columns name, chromosome and position)
#n.perm (the number of times each trait should be permutated)

###output
#a list with names: LOD, Effect, Trait_perm, Map, and Marker.

###Description
# map.1 is optimized for dataset in which many markers are un-informative
# these markers are not mapped but the p-vals and effect are copied from the previous marker.

QTL.map.1.perm <- function(trait.matrix,strain.map,strain.marker,n.perm){
                           if(missing(trait.matrix)){                                      stop("specify trait.matrix, matrix with traits in rows and genotypes in columns")}
                           if(missing(strain.map)){                                        stop("specify strain.map (genetic map)")}
                           if(ncol(strain.map) !=  dim(as.matrix(trait.matrix))[2] &
                              ncol(strain.map) !=  prod(dim(as.matrix(trait.matrix)))){    stop("number of strains is not equal to number of genotypes in map")}
                           if(missing(strain.marker)){                                     strain.marker <- data.frame(cbind(name=1:dim(strain.map)[1],chromosome=NA,bp=1:dim(strain.map)[1]))
                                                                                           rownames(strain.marker) <- 1:dim(strain.map)[1]}
                           if(missing(n.perm)){                                            n.perm <- 1000}
                           if(n.perm == 1 & dim(as.matrix(trait.matrix))[2] ==1){          stop("cannot conduct 1 permutation on 1 trait, set n.perm > 1")}


                           NAs <- sum(is.na(trait.matrix))>0 | sum(is.na(strain.map))>0

                           small <- FALSE
                           if(dim(as.matrix(trait.matrix))[1] ==1){
                               traitnames <- rownames(trait.matrix)
                               trait.matrix <- rbind(trait.matrix,trait.matrix)
                               rownames(trait.matrix) <- c(traitnames,"B")
                               small <- TRUE
                           }

                           traitnames <- rownames(trait.matrix)
                           trait.matrix <- trait.matrix[rep(1:nrow(trait.matrix),each=n.perm),]
                           trait.matrix <- permutate.traits(trait.matrix)
                           rownames(trait.matrix) <- paste(rep(traitnames,each=n.perm),1:n.perm,sep="_")

                           if(small){
                               trait.matrix <- trait.matrix[1:n.perm,]
                           }

                           eff.out <- matrix(NA,nrow(trait.matrix),nrow(strain.map))
                           pval.out <- matrix(NA,nrow(trait.matrix),nrow(strain.map))
                          
                           for(i in 1:nrow(strain.map)){
                               noseg <- length(unique(strain.map[i,])) ==1
                               if(noseg){
                                   output.tmp <- matrix(c(0,0),byrow=TRUE,ncol=2,nrow=nrow(trait.matrix))
                               }
                               if( i == 1 & !noseg){
                                   if(!NAs){output.tmp <- fast.lm.1(trait.matrix,strain.map[i,])}
                                   if(NAs){output.tmp <- lm.1(trait.matrix,strain.map[i,])}
                               }
                               if( i != 1 & sum(abs(strain.map[i-1,]-strain.map[i,]),na.rm=T) != 0 & !noseg){
                                   if(!NAs){output.tmp <- fast.lm.1(trait.matrix,strain.map[i,])}
                                   if(NAs){output.tmp <- lm.1(trait.matrix,strain.map[i,])}
                               }
                               if( i != 1 & sum(abs(strain.map[i-1,]-strain.map[i,]),na.rm=T) == 0 & !noseg){
                                   output.tmp <- output.tmp
                               }
                               eff.out[,i] <- output.tmp[,2]
                               pval.out[,i] <- output.tmp[,1]
                           }

                       colnames(eff.out) <- rownames(strain.marker)
                       rownames(eff.out) <- rownames(trait.matrix)
                       colnames(pval.out) <- rownames(strain.marker)
                       rownames(pval.out) <- rownames(trait.matrix)


                       output <- NULL; output <- as.list(output)
                       output[[1]] <- round(pval.out,digits=2)
                       output[[2]] <- round(eff.out,digits=3)
                       output[[3]] <- trait.matrix
                       output[[4]] <- strain.map
                       output[[5]] <- strain.marker
                       names(output) <- c("LOD","Effect","Trait_perm","Map","Marker")
                       
                       return(output)                           
                      }
