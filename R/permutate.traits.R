###permutates a trait matrix

###input
#trait.matrix (matrix strains in columsn, traits in rows)


###output
#by-row permutated matrix

###Description


    permutate.traits <- function(trait.matrix){
                                 perm.list <- rep(c(1:nrow(trait.matrix)),times=ncol(trait.matrix)) + runif(length(trait.matrix),0.25,0.75)
                                 perm.list <- order(perm.list)
                                 perm.data <- matrix(as.numeric(trait.matrix)[perm.list],nrow=nrow(trait.matrix),ncol=ncol(trait.matrix),byrow=TRUE)
                                 return(perm.data)
    }

    