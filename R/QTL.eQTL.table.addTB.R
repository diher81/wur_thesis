###Adds whether a QTL is located in a trans-band (or cis-enriched area)

###input
#eQTL.table.file (output from QTL.eQTL.table)
#eQTL.call.TBs.file (output from QTL.eQTL.call.TBs)
#transband_thr (statistical threshold for calling a trans-band)
#qtl_type_select (type of qtl to call enrichment for)
#merge_condition (how adjacent chromosomal intervals should be for merging trans-bands)
#digits_Mb (amount of digits behind the Mb should appear in TB-name)

###output
#the eQTL.table.file with a trans_band column added

###Description

###version
#Fixed issue with merging transbands

QTL.eQTL.table.addTB <- function(eQTL.table.file,call.transbands.file,transband_thr,qtl_type_select,merge_condition,digits_Mb){                               
                                   if(missing(eQTL.table.file)){                        stop("need a eQTL.table.file")}
                                   if(missing(call.transbands.file)){                   stop("need a call.transbands.file")}
                                   if(missing(transband_thr)){                          transband_thr <- 0.0001; warning("set default significance, might not be optimal")}
                                   if(missing(qtl_type_select)){                        qtl_type_select <- "trans"}
                                   if(missing(merge_condition)){                        merge_condition <- 2}
                                   if(missing(digits_Mb)){                              digits_Mb <- 1}                
                                    
                                   transbands <- dplyr::filter(call.transbands.file,transband_significance < transband_thr,qtl_type==qtl_type_select) %>%
                                                 dplyr::group_by(qtl_chromosome) %>%
                                                 dplyr::mutate(mergeornot=c(0,diff(Interval))) %>%
                                                 data.frame()
                                   
                                   transbands.merged <- NULL
                                   for(i in 1:length(unique(transbands$qtl_chromosome))){
                                       tmp <- dplyr::filter(transbands,qtl_chromosome == unique(transbands$qtl_chromosome)[i])
                                       tmp2 <- NULL
                                       if(nrow(tmp) > 1){
                                           for(j in 2:nrow(tmp)){
                                               if(tmp$mergeornot[j] <= merge_condition){
                                                   tmp2 <- data.frame(rbind(tmp2,cbind(qtl_chromosome=tmp[j,1], qtl_type=tmp[j,3], n.ct=sum(tmp[c((j-1),j),4]), exp.ct=sum(tmp[c((j-1),j),5]), 
                                                                                       transband_significance=max(tmp[j,6]), Interval_left=tmp[j-1,8], Interval_right=tmp[j,9])))
                                                   for(k in c(3,4,5,6,7)){tmp2[,k] <- as.numeric(as.character(unlist(tmp2[,k])))}  
                                               } else {
                                                   if(j==2){
                                                       tmp3 <- tmp[1:2,-c(2,7,10)]
                                                       for(k in c(3,4,5,6,7)){tmp3[,k] <- as.numeric(as.character(unlist(tmp3[,k])))}  
                                                       tmp2 <- rbind(tmp2,tmp3)
                                                       for(k in c(3,4,5,6,7)){tmp2[,k] <- as.numeric(as.character(unlist(tmp2[,k])))}  
                                                   } else {
                                                       tmp3 <- tmp[j,-c(2,7,10)]
                                                       for(k in c(3,4,5,6,7)){tmp3[,k] <- as.numeric(as.character(unlist(tmp3[,k])))}  
                                                       tmp2 <- rbind(tmp2,tmp3)
                                                       for(k in c(3,4,5,6,7)){tmp2[,k] <- as.numeric(as.character(unlist(tmp2[,k])))}  
                                                   }
                                               }
                                           }
                                       } else {
                                           tmp2 <- tmp[,-c(2,7,10)]
                                           for(k in c(3,4,5,6,7)){tmp2[,k] <- as.numeric(as.character(unlist(tmp2[,k])))}  
                                       }
                                       transbands.merged <- data.frame(rbind(transbands.merged,tmp2))
                                       for(k in c(3,4,5,6,7)){transbands.merged[,k] <- as.numeric(as.character(unlist(transbands.merged[,k])))}  
                                   }
                                   for(k in c(1,2)){transbands.merged[,k] <- as.character(unlist(transbands.merged[,k]))}  
                                   transbands.merged <- mutate(transbands.merged,Name=paste(qtl_chromosome,":",round(Interval_left/1e6,digits=digits_Mb),"-",round(Interval_right/1e6,digits=digits_Mb),"Mb",sep=""))
                                   
                                   output <- dplyr::mutate(eQTL.table.file,trans_band = "none")
                                   
                                   for(i in 1:nrow(transbands.merged)){
                                        output[output$qtl_type == transbands.merged$qtl_type[i] &
                                               output$qtl_chromosome == transbands.merged$qtl_chromosome[i] &
                                               output$qtl_bp > transbands.merged$Interval_left[i] &
                                               output$qtl_bp <= transbands.merged$Interval_right[i] &
                                               !is.na(output$qtl_type) & !is.na(output$qtl_chromosome) & !is.na(output$qtl_bp),]$trans_band <- transbands.merged$Name[i]
                                   }
                                   
                                   return(output)
                                  }
