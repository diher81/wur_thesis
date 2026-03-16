###Generates a color scale (transcriptomics auxiliary)

###input
#

###output
#date

###Description




###Generates colors for a gradient
color.scale <- function(colors.use,y.val){
                        if(class(y.val) == "matrix"){
                            scale.use <-  c(as.numeric(rownames(y.val))[nrow(y.val)],as.numeric(rownames(y.val))[1],as.numeric(rownames(y.val))[1]-as.numeric(rownames(y.val))[2])
                        }
                        if(class(y.val) == "numeric"){
                            scale.use <- c(min(y.val,na.rm=T),max(y.val,na.rm=T),diff(c(min(y.val,na.rm=T),max(y.val,na.rm=T)))/100)
                        }
                        seq1 <- seq(scale.use[1],scale.use[2],scale.use[3])    ###adjust for optimal effect
                        ramp <- colorRamp(colors.use, bias=1)
                        colors.out <- rgb(ramp(seq(0, 1, length = length(seq1))), max = 255)
                        return(colors.out)
                       }
