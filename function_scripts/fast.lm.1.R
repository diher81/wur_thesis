###fast lm function

###input
#traits (matrix)
#variable (vector)

###output
#matrix with -log10 transformed p-values and effects

###Description
# Fast linear model, 1 variable
# Uses fast.summ, can't handle NA's!

    fast.lm.1 <- function(traits,variable){
                          if(missing(traits)){        stop("specify traits")}
                          if(missing(variable)){      stop("specify variable")}

                          output <- fast.summ(lm(t(traits)~variable))
                          return(output)
                         }

