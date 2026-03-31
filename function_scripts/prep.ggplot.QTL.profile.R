    prep.ggplot.QTL.profile <- function(QTL.peak.output,map1.output,trait){
                                        output <- list()

                                        output[[1]] <- QTL.peak.output[as.character(unlist(QTL.peak.output$trait))==trait,]
                                            names(output)[1] <- "QTL_profile"
                                        tmp <- output[[1]]$qtl_peak[!is.na(output[[1]]$qtl_peak)]

                                        if(length(tmp)>0){
                                            for(i in 1:length(tmp)){
                                                output[[i+1]] <- cbind(trait=trait,
                                                                       map1.output$Marker[tmp[i],],
                                                                       strain=colnames(map1.output$Trait),
                                                                       genotype=map1.output$Map[tmp[i],],
                                                                       trait_value=map1.output$Trait[rownames(map1.output$Trait)==trait,],
                                                                       R_squared=R.sq(x=map1.output$Map[tmp[i],],y=map1.output$Trait[rownames(map1.output$Trait)==trait,]))
                                                names(output)[i+1] <- paste("QTLsplit_",as.character(unlist(map1.output$Marker[tmp[i],1])),sep="_")
                                            }
                                        }
                                        return(output)
                                       }
