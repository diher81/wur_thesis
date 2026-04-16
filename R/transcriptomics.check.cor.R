################################################################################################################################################################
### correlation.check(trans.int,filename) v1.2
###   returns - the mean of the correlations of sample x with all other samples
###           - a plot of the correlations of sample x with all other samples, with an indication which samples are < 0.9 (saved as pdf) and a heatmap of the correlations (saved as pdf)
###   obligatory input
###     trans.int, output from transform.int function. Or any list file which contains a matrix under $log2.int
transcriptomics.check.cor <- function(trans.int,filename,save_dir){
                                     if(missing(trans.int)){                          stop("Require output of transform.int")}
                                     if(missing(filename)){                           filename <- ""}
                                     if(missing(save_dir)){                           save_dir <- getwd()}
      
                                     ###Correlation
                                     correlations <- cor(trans.int$log2_intensities)
                                         write.table(correlations,file=paste(save_dir,date.now(),"_correlations.txt",sep=""),sep="\t",quote=F)
                                     correlsum <- apply(correlations,1,mean)
      
                                     if(filename==""){pdf(file=paste(save_dir,"/",date.now(),"_correlation.check.pdf",sep=""),width=20,height=9)}
                                     if(filename!=""){pdf(file=paste(save_dir,"/",date.now(),"_",filename,"_correlation.check.pdf",sep=""),width=20,height=9)}
                                         ###Log2int sample correlations
                                         axis.nu <- axis.output(1,length(correlsum))
                                         par(fig=c(0.3,0.7,0.1,0.9))
                                         plot(sort(correlsum),pch=19,axes=FALSE,xlab="",ylab="",ylim=c(0.75,1))
                                             axis(1,axis.nu[[2]],axis.nu[[2]])
                                             mtext(side=1,text="Samples",cex=1.5,line=3)
                                             axis(2,seq(0.75,1,by=0.05),seq(0.75,1,by=0.05),las=2)
                                             mtext(side=2,text="Mean rho",cex=1.5,line=3)
                                             abline(h=0.9,lwd=2,col="grey",lty=2)
                                             if(sum(correlsum < 0.9) > 0){
                                                 for(i in 1:sum(correlsum < 0.9)){
                                                     lines(x=c(i,length(correlsum)/4),y=c(sort(correlsum)[sort(correlsum)<0.9][i],seq(0.745,0.895,length=sum(correlsum < 0.9))[i]))
                                                     text(x=length(correlsum)/4,y=seq(0.745,0.895,length=sum(correlsum < 0.9))[i],labels=colnames(trans.int$log2_intensities)[order(correlsum)][i],pos=4,cex=0.7)
                                                 }
                                             }
                                         box()
      
                                         ###ratio mean sample correlations
                                         correlations <- cor(trans.int$log2_ratio_mean)
                                         colscale <- color.scale(colors.use=c("#8E0152","#C51B7D","#DE77AE","#F1B6DA","#FDE0EF","#F7F7F7","#E6F5D0","#B8E186","#7FBC41","#4D9221","#276419"),
                                                                 y.val=as.numeric(correlations))
      
                                         heatmap(correlations,scale="none",col=colscale,margins=c(10,8))
                                         par(fig=c(0,0.3,0.5,1),new=T)
                                         scale.legend(input=as.numeric(correlations),start.zero=F,col.scale=colscale,line.col="black")
                                     dev.off()
      
                                     return(correlsum)
                                    }