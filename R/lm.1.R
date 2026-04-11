###lm function for NA handling

###input
#data.input (matrix)
#variable (vector)

###output
#matrix with -log10 transformed p-values and effects

###Description
# Slow linear model for NA handling
# ignores heterozygous strains in effect determination



    lm.1 <- function(data.input,variable){
                     retext <- function(x){return(c(first(x),last(x)))}
                     
                     out <- NULL
                     for(i in 1:nrow(data.input)){
                         out <- rbind(out,
                                      cbind(LOD=-log10(cor.test(data.input[i,!is.na(variable)],variable[!is.na(variable)],na.action=na.exclude())$p.value),
                                            effect=diff(retext(tapply(data.input[i,!is.na(variable)],variable[!is.na(variable)],mean,na.rm=T))))
                                     )
                     }

                     return(out)
                    }

