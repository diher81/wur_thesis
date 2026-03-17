################################################################################################################################################################
### uses bioconductor limma package
### agilent.limma.norm(targets,bg.method,nWA.method,nBA.method,bad.arrays, filename) v1.2
###   returns normalized data (matrix), also saves the normalized data files
###   obligatory input
###     targets: a target file, filename of FE output txt + Cy3 and Cy5 sample description
###   optional input
###     bg.method, standard is "none"
###     nWA.method, standard is "loess"
###     nBA.method, standard is "quantile"
###     bad.arrays, standard is empty. Give the sample number that should be removed (Cy3 comes before Cy5)
###     filename, output file name, default is "obj_MAnBWABA.[date].out"
### v1.1 the option of choosing a location for array loading other than the standard wd is added
### v1.2 the standard save location is now the workwd
### v1.21 function returns to original workwd

normalize.agilent.transcriptomics <- function(targets,bg.method,nWA.method,nBA.method,bad.arrays,array_dir,save_dir,filename){
                                if(missing(targets)){                                   stop("Require targets for loading MA data")}
                                if(missing(bg.method)){                                 bg.method <- "none"}
                                if(missing(nWA.method)){                                nWA.method <- "loess"}
                                if(missing(nBA.method)){                                nBA.method <- "quantile"}
                                if(missing(bad.arrays)){                                rem.bad <- FALSE}
                                if(missing(bad.arrays)==FALSE){                         rem.bad <- TRUE}
                                if(missing(array_dir)){                                 array_dir <- getwd()}
                                if(missing(save_dir)){                                  save_dir <- getwd()}
                                if(missing(filename)){                                  filename <- ""}

                                workwd <- getwd()
                                setwd(array_dir)
                              	     RG <<- read.maimages(targets$FileName, source="agilent",columns=list(G = "gMeanSignal",Gb = "gBGMedianSignal", R= "rMeanSignal", Rb = "rBGMedianSignal"))
                                setwd(save_dir)
                                ###Remove bad arrays
                                if(rem.bad){
                                    bad.arrays <- sort(bad.arrays)
                                    ###only remove if both channels of an array are bad: otherwise the normalizeWithin does not work
                                    select <- c(bad.arrays[bad.arrays<=ncol(RG$G)] %in% (bad.arrays[bad.arrays>ncol(RG$R)]-ncol(RG$R)),
                                                (bad.arrays[bad.arrays>ncol(RG$R)]-ncol(RG$R)) %in% bad.arrays[bad.arrays<=ncol(RG$G)])
                                    bad.arrays.single <- bad.arrays[select==FALSE]
                                    bad.arrays.PreNorm <- bad.arrays[select]

                                    for(i in 1:length(bad.arrays.single)){
                                        bad.arrays.single[i] <- bad.arrays.single[i] - sum(bad.arrays.PreNorm < bad.arrays.single[i])
                                    }
                                    ###condition needed, otherwise if 0, everything will be removed
                                    if(length(bad.arrays.PreNorm)>0){
                                        RG$G <- RG$G[,-bad.arrays.PreNorm[bad.arrays.PreNorm<=ncol(RG$G)]]
                                        RG$Gb <- RG$Gb[,-bad.arrays.PreNorm[bad.arrays.PreNorm<=ncol(RG$Gb)]]
                                        RG$R <- RG$R[,-(bad.arrays.PreNorm[bad.arrays.PreNorm>ncol(RG$R)]-ncol(RG$R))]
                                        RG$Rb <- RG$Rb[,-(bad.arrays.PreNorm[bad.arrays.PreNorm>ncol(RG$Rb)]-ncol(RG$Rb))]
                                    }
                                }

                                RGnb <- backgroundCorrect(RG,method=bg.method)
                                MAnbWA <- normalizeWithinArrays(RGnb,method=nWA.method)
                                MAnbWABA <- normalizeBetweenArrays(MAnbWA,method=nBA.method)

                                if(filename==""){save(MAnbWABA,file=paste(save_dir,"/obj_MAnBWABA.",date.now(),".Rdata",sep=""))}
                                if(filename!=""){save(MAnbWABA,file=paste(save_dir,"/obj_MAnBWABA.",filename,".",date.now(),".Rdata",sep=""))}

                                rg.int <- RG.MA(MAnbWABA)
                                rg.norm <- cbind(rg.int$G,rg.int$R)

                                cy3.ids <- as.character(targets$Cy3)
                                cy5.ids <- as.character(targets$Cy5)
                                if(rem.bad){
                                    if(!is.na(bad.arrays.single[1])){
                                        rg.norm <- rg.norm[,-bad.arrays.single]
                                    }
                                    colnames(rg.norm) <- paste(c(1:(nrow(targets)*2)),c(cy3.ids,cy5.ids),sep=":")[-bad.arrays]}
                                if(rem.bad==F){
                                    colnames(rg.norm) <- paste(c(1:(nrow(targets)*2)),c(cy3.ids,cy5.ids),sep=":")}

                                if(filename==""){save(rg.norm,file=paste(save_dir,"/obj_rg.norm.",date.now(),".Rdata",sep=""))}
                                if(filename!=""){save(rg.norm,file=paste(save_dir,"/obj_rg.norm",".",filename,".",date.now(),".Rdata",sep=""))}

                                setwd(workwd)
                                return(rg.norm)
                               }
