###function for making an eQTL table

###input
#QTL.peak.dataframe (output from QTL.peak.finder function)
#trait.annotation (dataframe with columns: trait, gene_chromosome, gene_bp, gene_database_ID, gene_name_1, gene_name_2)
#cis.trans.window (how cis and trans eQTL are distinguised: physical, LOD.drop, or both)
#cis.trans.window.size (physical window size, standard 1000000)
#colnames.trait.annotations (names for the trait annotations added)


###output
#a dataframe listing the eQTL per trait.

###Description    
# Makes a table with eQTL peaks v1.1
# QTL.peak.dataframe: output from peak.finder
# trait.annotation: gene ID's of the mapped traits
# cis.trans.window: Cis-trans windowd can be defined based on physical distance "physical", based on LOD-drop interval "LOD.drop", or on both "both"
    
QTL.eQTL.table <- function(QTL.peak.dataframe, trait.annotation,cis.trans.window,cis.trans.window.size,colnames.trait.annotations){
                           if(missing(QTL.peak.dataframe)){                                  stop("need an eQTL list file")}
                           if(missing(trait.annotation)){                                    stop("need annotation file with order: trait, gene_chromosome, gene_bp, gene_WBID, gene_sequence_name, gene_public_name")}
                           if(missing(cis.trans.window)){                                    cis.trans.window <- "physical"} 
                           if(!cis.trans.window %in% c("physical","LOD.drop","both")){       stop("cis.trans.window should be in physical, LOD.drop, or both")}
                           if(missing(cis.trans.window.size)){                               cis.trans.window.size <- 1000000}
                           if(missing(colnames.trait.annotations)){                          colnames.trait.annotations <- c("trait","gene_chromosome","gene_bp","gene_WBID","gene_sequence_name","gene_public_name")}

                           ###Rename colnames traits
                           colnames(trait.annotation) <- colnames.trait.annotations

                           ###Take only the peaks
                           if(cis.trans.window =="physical"){
                               output <- dplyr::filter(QTL.peak.dataframe,!is.na(qtl_peak)) %>%
                                         merge(trait.annotation,by.x=1,by.y=1) %>%
                                         dplyr::mutate(qtl_type = ifelse((qtl_chromosome==gene_chromosome & abs(qtl_bp-gene_bp) < cis.trans.window.size),"cis","trans"))
                           }
                      
                           if(cis.trans.window =="LOD.drop"){
                               output <- dplyr::filter(QTL.peak.dataframe,!is.na(qtl_peak)) %>%
                                         merge(trait.annotation,by.x=1,by.y=1) %>%
                                         dplyr::mutate(qtl_type = ifelse((qtl_chromosome==gene_chromosome & (gene_bp > qtl_bp_left & gene_bp < qtl_bp_right)),"cis","trans"))
                           }
                           
                           if(cis.trans.window =="both"){
                               output <- dplyr::filter(QTL.peak.dataframe,!is.na(qtl_peak)) %>%
                                         merge(trait.annotation,by.x=1,by.y=1) %>%
                                         dplyr::mutate(qtl_type = ifelse((qtl_chromosome==gene_chromosome & abs(qtl_bp-gene_bp) < cis.trans.window.size) |
                                                                         ((qtl_chromosome==gene_chromosome & abs(qtl_bp-gene_bp) < cis.trans.window.size) & (gene_bp > qtl_bp_left & gene_bp < qtl_bp_right)),"cis","trans"))
                           }
                           
                           ###Return file
                           return(output)
                          }
