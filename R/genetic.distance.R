###calculate genetic distance between markers


###input
#strain.map (of the investigated strains: markers in rows, genotypes in columns)
#strain.marker (of the investigated strains: markers in rows, with columns name, chromosome and position)
#chromosomes (character vector with chromosome names; standard is C. elegans)
#chromosome_size (numeric vector with chromosome sizes; standard is C. elegans)
#ignore.na (logic, tells if NAs should be incorporated in the analysis)


###output
#a dataframe with chromosome, physical position, and genetic position (cM)

###Description
#ignore.na makes the algorithm find te nearest non-NA marker


genetic.distance <- function(strain.map, strain.marker, chromosomes, chromosome_size, ignore.na){
                            if(missing(strain.map)){                             
                              stop("Requires genetic map with markers per row and strains per column")}
                            if(missing(strain.marker)){                          
                              stop("Requires marker file with 3 columns: name, chromosome, position")}
                            if(missing(chromosomes)){                            
                              chromosomes <- c("I","II","III","IV","V","X")}
                            if(missing(chromosome_size)){                        
                              chromosome_size <- c(15072434,15279421,13783801,17493829,20924180,17718942)}
                            if(length(chromosomes) != length(chromosome_size)){  
                              stop("chromosome names and sizes should have the same length")}
                            if(missing(ignore.na)){                              
                              ignore.na <- FALSE}

                            ###rolling mean function
                            rollingmean <- function(x){
                              out <- NULL
                              for(i in 2:length(x)){
                                out <- c(out,mean(x[(i-1):i]))
                                }
                              return(out)
                              }

                            ###make colnames comply to format
                            colnames(strain.marker) <- c("name", "chromosome", "position")

                            ###output
                            output <- NULL

                            ###per chromosome
                            for(i in 1:length(chromosomes)){
                                currchr <- strain.map[strain.marker[,2] == chromosomes[i],]
                                currmrk <- strain.marker[strain.marker[,2] == chromosomes[i],]
                                ###fills in gaps based on the genotype of the closest non-NA marker
                                if(!ignore.na){
                                    for(j in 1:ncol(currchr)){
                                        for(k in which(is.na(currchr[,j]))){
                                            currori <- currchr[,j]
                                            closest <- currmrk[,3]-currmrk[k,3]
                                            closest[k] <- NA
                                            closest <- order(abs(closest))
                                            ii <- 1
                                            while(is.na(currori[closest[ii]])){
                                                ii <- ii + 1
                                            }
                                            currchr[k,j] <- currori[closest[ii]]
                                        }
                                    }
                                }

                                ###count recombinations
                                diffmap <- apply(currchr, 2, diff) != 0

                                notrecomb <- apply(!diffmap & !is.na(diffmap), 1, sum)
                                yesrecomb <- apply(diffmap & !is.na(diffmap), 1, sum)

                                fractions <- cumsum(yesrecomb/(notrecomb+yesrecomb))
                                fractions <- c(fractions[1], fractions, rev(fractions)[1])

                                locations <- rollingmean(currmrk[,3])
                                locations <- c(1, locations, chromosome_size[i])

                                ###into output
                                output <- rbind(output,
                                                cbind(chromosome = chromosomes[i], 
                                                      position = locations, 
                                                      cM = fractions))
                            }
                            output <- data.frame(output)
                                output[,1] <- as.character(unlist(output[,1]))
                                output[,2] <- as.numeric(as.character(unlist(output[,2])))
                                output[,3] <- as.numeric(as.character(unlist(output[,3])))*100

                            return(output)
                           }

