###Generates a nice axis scale (transcriptomics plotting auxiliary)

###input
#min.nu
#max.nu

###output
#small and large axis

###Description


axis.output <- function(min.nu,max.nu){
                        output <- NULL; output <- as.list(output)
                        for(i in sort(c(10^c(seq(-2,8,by=1)),2*10^c(seq(-2,8,by=1)),4*10^c(seq(-2,8,by=1)),6*10^c(seq(-2,8,by=1)),8*10^c(seq(-2,8,by=1))),decreasing=T)){
                            ###determine the scale of the differences
                            if(max.nu-min.nu < i*2){
                                small.axis <- seq(floor(min.nu/i)*i,ceiling(max.nu/i)*i*2,by=i/20)
                                large.axis <- seq(floor(min.nu/i)*i,ceiling(max.nu/i)*i*2,by=i/4)[-1]
                            }
                        }
                        output[[1]] <- small.axis
                        output[[2]] <- large.axis
                        return(output)
                       }