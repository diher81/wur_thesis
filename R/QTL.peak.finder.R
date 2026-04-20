###function for single marker mapping
#Snoek & Sterken, 2017; Sterken, 2017

###input
#map1.dataframe; output of the QTL.map1.dataframe function
#threshold; the significance treshold for calling QTL (in -log10(p))
#LOD.drop; the decrease in LOD-score compared to the peak to call the boundaries of the QTL. Standard is 1.5

###output
#the original dataframe with the added columns: qtl_peak, qtl_marker_left, qtl_bp_left, qtl_marker_right, qtl_bp_right
#these columns are only non-NA where a QTL peak is located.

###Description
#This takes output of mapping.to.list function and identifies peaks and confidence intervals v1.3
#Based on a given threshold and a 1.5 LOD-drop.
#v1.2 made the function more efficient, the peak borders are selected in one processing
#v1.3 fixed the problem calling the peak if the peak borders the edge of the chromosome
#v1.4 fixed 'larger then, smaller then' bug
#v1.5 fixed 'incompatible types' bug; NA is not recognized as numeric
#v1.6 rewritten the code, so it is fully in dplyr framework

QTL.peak.finder <- function(map1.dataframe,threshold,LOD.drop){
                            if(missing(map1.dataframe)){            stop("Missing map1.dataframe (trait, qtl_chromosome, qtl_bp, qtl_marker, qtl_significance, qtl_effect)")}
                            if(missing(threshold)){            stop("Missing threshold")}
                            if(missing(LOD.drop)){             LOD.drop <- 1.5}

                            ###locate the peak
                            output <- dplyr::group_by(map1.dataframe,trait,qtl_chromosome) %>%
                                      dplyr::mutate(qtl_peak_tmp=ifelse(qtl_significance == max(qtl_significance) & max(qtl_significance)>=threshold,1,NA)) %>%
                                      dplyr::group_by(trait,qtl_chromosome,qtl_peak_tmp) %>%
                                      dplyr::mutate(qtl_peak=ifelse(is.na(qtl_peak_tmp),NA,
                                                               ifelse(abs(qtl_marker - round(mean(qtl_marker),digits = 0)) == min(abs(qtl_marker - round(mean(qtl_marker),digits = 0))),qtl_marker,NA))) %>%
                                      dplyr::group_by(trait,qtl_chromosome) %>%
                                      dplyr::mutate(Border=ifelse((qtl_significance < max(qtl_significance)-LOD.drop & max(qtl_significance)>=threshold),T,
                                                             ifelse((qtl_bp == min(qtl_bp) & min(qtl_bp) %in% qtl_bp[qtl_significance >= max(qtl_significance)-LOD.drop & max(qtl_significance)>=threshold]),T,
                                                               ifelse((qtl_bp == max(qtl_bp) & max(qtl_bp) %in% qtl_bp[qtl_significance >= max(qtl_significance)-LOD.drop & max(qtl_significance)>=threshold]),T,F))),
                                                    peakloc=ifelse(sum(!is.na(qtl_peak)) == 1,qtl_peak[!is.na(qtl_peak)],NA)) %>%
                                      dplyr::group_by(trait,qtl_chromosome,Border) %>%
                                      dplyr::mutate(qtl_marker_left=as.numeric(ifelse(Border,ifelse(qtl_marker == max(qtl_marker[qtl_marker <= peakloc]),qtl_marker,NA),NA)),
                                                    qtl_bp_left=as.numeric(ifelse(Border,ifelse(qtl_bp == max(qtl_bp[qtl_marker <= peakloc]),qtl_bp,NA),NA)),
                                                    qtl_marker_right=as.numeric(ifelse(Border,ifelse(qtl_marker == min(qtl_marker[qtl_marker >= peakloc]),qtl_marker,NA),NA)),
                                                    qtl_bp_right=as.numeric(ifelse(Border,ifelse(qtl_bp == min(qtl_bp[qtl_marker >= peakloc]),qtl_bp,NA),NA))) %>% 
                                      dplyr::group_by(trait,qtl_chromosome) %>%
                                      dplyr::mutate(qtl_peak_tmp=ifelse(sum(!is.na(qtl_peak))==1,1,0)) %>%
                                      dplyr::group_by(trait,qtl_chromosome,qtl_peak_tmp) %>%
                                      dplyr::mutate(qtl_marker_left=ifelse(qtl_peak_tmp==1,qtl_marker_left[!is.na(qtl_marker_left)],NA),
                                                    qtl_bp_left=ifelse(qtl_peak_tmp==1,qtl_bp_left[!is.na(qtl_bp_left)],NA),
                                                    qtl_marker_right=ifelse(qtl_peak_tmp==1,qtl_marker_right[!is.na(qtl_marker_right)],NA),
                                                    qtl_bp_right=ifelse(qtl_peak_tmp==1,qtl_bp_right[!is.na(qtl_bp_right)],NA)) %>%
                                      dplyr::ungroup() %>%
                                      dplyr::select(trait,qtl_chromosome,qtl_bp,qtl_marker,qtl_significance,qtl_effect,qtl_peak,qtl_marker_left,qtl_bp_left,qtl_marker_right,qtl_bp_right) %>%
                                      dplyr::mutate(qtl_marker_left=ifelse(!is.na(qtl_peak),qtl_marker_left,NA),
                                                    qtl_bp_left=ifelse(!is.na(qtl_peak),qtl_bp_left,NA),
                                                    qtl_marker_right=ifelse(!is.na(qtl_peak),qtl_marker_right,NA),
                                                    qtl_bp_right=ifelse(!is.na(qtl_peak),qtl_bp_right,NA))
                        
                            ###return called peaks
                            return(output)
                           }

