###calculates R-squared as in a linear model

###input
#x (vector)
#y (vector)


###output
#squared correlation coeffcient

###Description

    R.sq <- function(x,y){return(cor(x,y,use="na.or.complete")^2)}
    
      