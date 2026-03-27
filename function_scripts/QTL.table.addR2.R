###function for determining variance explained by QTL

###input
#QTL.table.file (should contain marker locations: qtl_marker, qtl_marker_1, qtl_marker_2 and trait names)
#QTL.list.file (list with QTL data, should contain entries Map and Trait)

###output
#the R-squared is added to the QTL.table.file

###Description

# This function will add the sum of squares of the used model to a QTL table output, it also calculates the r-squared of the model


    QTL.table.addR2 <- function(QTL.table.file,QTL.list.file){
                                 if(missing(QTL.table.file)){                          stop("Give a QTL summary table, should contain marker locations (qtl_marker, qtl_marker_1, qtl_marker_2) and trait names")}
                                 if(missing(QTL.list.file)){                           stop("Give the input on which was mapped (output of QTL.data.prep)")}

                                 ###R-squared function
                                 R.sq <- function(x,y){return(cor(x,y,use="na.or.complete")^2)}
                                 
                                 ###True if single marker mapping
                                 if("qtl_marker" %in% colnames(QTL.table.file)){
                                     output <- NULL
                                     for(i in 1:nrow(QTL.table.file)){
                                         ###get marker
                                         mrk <- QTL.list.file$Map[as.numeric(as.character(unlist(QTL.table.file$qtl_marker[i]))),]
                                         ###get expression
                                         y <- unlist(QTL.list.file$Trait[rownames(QTL.list.file$Trait)==as.character(unlist(QTL.table.file$trait[i])),])
                                         ###run model
                                         output <- c(output,R.sq(y,mrk))
                                     }
                                     output <- data.frame(cbind(QTL.table.file,qtl_R2_sm=output))
                                 }

                                 ###True if multiple marker mapping
                                 if("marker_1" %in% colnames(QTL.table.file) & "marker_2" %in% colnames(QTL.table.file)){
                                     output <- cbind(QTL.table.file,SS_int_mrk1=NA,SS_int_mrk2=NA,SS_int=NA,SS_int_res=NA,SS_add_mrk1=NA,SS_add_mrk2=NA,SS_add_res=NA)
                                     for(i in 1:nrow(QTL.table.file)){
                                         ###get marker
                                         mrk1 <- QTL.list.file$Map[QTL.table.file$marker_1[i],]
                                         mrk2 <- QTL.list.file$Map[QTL.table.file$marker_2[i],]
                                         ###get expression
                                         y <- unlist(QTL.list.file$Trait[rownames(QTL.list.file$Trait)==QTL.table.file$trait[i],])
                                         ###run model
                                         output[i,colnames(output)%in%c("SS_int_mrk1","SS_int_mrk2","SS_int","SS_int_res")] <- anova(lm(y ~ mrk1*mrk2))$"Sum Sq"
                                         output[i,colnames(output)%in%c("SS_add_mrk1","SS_add_mrk2","SS_add_res")] <- anova(lm(y ~ mrk1+mrk2))$"Sum Sq"

                                     }
                                     output <- dplyr::mutate(output,r2_add=(SS_add_mrk1+SS_add_mrk2)/(SS_add_mrk1+SS_add_mrk2+SS_add_res),r2_int=SS_int/(SS_int_mrk1+SS_int_mrk2+SS_int+SS_int_res))
                                 }
                                 return(output)
                                }

