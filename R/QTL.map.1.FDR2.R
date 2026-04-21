###FDR calculation function for single marker mapping
#Snoek & Sterken, 2017; Sterken, 2017

###input
#map1.output (output from the QTL.map.1 function)
#filenames.perm (files containing permutated data)
#FDR_dir (location of the files with permutated data)
#q.value (desired q-value)
#small (logic; set to TRUE for a small set)

###output
#a list with names: LOD, Effect, Trait_perm, Map, and Marker.

###Description
# Based on Benjamini and Yakuteili multiple testing under dependency



QTL.map.1.FDR2 <- function(map1.output,
                           filenames.perm,
                           FDR_dir,
                           q.value,
                           small) {
  if (missing(map1.output)) {
    stop("Missing map1.output ([[LOD]],[[Effect]],[[Trait]],[[Map]],[[Marker]])")
  }
  if (missing(filenames.perm)) {
    stop("Need a vector with filenames for permutation")
  }
  if (missing(FDR_dir)) {
    FDR_dir <- getwd()
  }
  if (missing(q.value)) {
    q.value <- 0.025
  }
  if (missing(small)) {
    small <- FALSE
  }
  
  output <- as.list(NULL)
  
  output[[2]] <- matrix(0, ncol = length(filenames.perm), nrow =
                          5000)
  rownames(output[[2]]) <- seq(0.1, 500, by = 0.1)
  colnames(output[[2]]) <- filenames.perm
  
  ### In case of one trait (small = TRUE)
  if (small) {
    if (length(filenames.perm) > 1) {
      stop("needs one permutation file, with n permutations")
    }
    
    # Load permutation object
    obj_name <- load(file = paste(FDR_dir, filenames.perm, sep = ""))
    obj <- get(obj_name)
    
    # Extract LOD matrix (should be n_perm x 1 or n_perm x n_traits)
    tmp <- obj$LOD
    
    # Handle NAs
    tmp[is.na(tmp)] <- 0.1
    
    # Convert to vector if it's a matrix with one column
    if (is.matrix(tmp) && ncol(tmp) == 1) {
      tmp <- as.vector(tmp)
    }
    
    # CRITICAL FIX: Keep all permutation max LODs, don't take max()
    output[[2]] <- tmp  # Now a vector of length = number of permutations
    
    # Warning if too few permutations for stable FDR
    if ((length(output[[2]]) * q.value) < 10) {
      warning(
        paste(
          "inadvisable to calculate q=",
          q.value,
          "based on ",
          length(output[[2]]),
          " permutations. Increase to at least ",
          ceiling(10 / q.value)
        )
      )
    }
    
    # Calculate threshold: (1 - q.value) quantile of null max LOD distribution
    # Multiply by 10 for consistency with eQTL methodology (per 0.1 LOD unit)
    sorted_nulls <- rev(sort(output[[2]]))
    idx <- ceiling(q.value * length(output[[2]]))
    output[[1]] <- sorted_nulls[idx] * 10
    
    # Average false discovery (for reporting)
    output[[3]] <- mean(output[[2]], na.rm = TRUE)
    
    # Real data max LOD
    real_max <- max(map1.output$LOD, na.rm = TRUE)
    output[[4]] <- real_max
    
    names(output) <- c(
      "Significance_threshold",
      "Null_maxLOD_distribution",
      "Null_maxLOD_mean",
      "Observed_maxLOD"
    )
  }
  
  if (!small) {
    for (i in 1:length(filenames.perm)) {
      tmp <- load(file = paste(FDR_dir, filenames.perm[i], sep = ""))
      ###No NA's
      tmp <- get(tmp)[[1]]
      tmp[is.na(tmp)] <- 0.1
      
      tmp <- table(apply(round(tmp, digits =
                                 1), 1, max, na.rm = T))
      
      output[[2]][(as.numeric(names(tmp)) * 10), i] <- tmp
      output[[2]][, i] <- rev(cumsum(rev(output[[2]][, i])))
    }
    ###Take the average over permutations
    output[[3]] <- apply(output[[2]], 1, mean)
    
    ###Calculate the real discoveries, no NA's
    tmp <- map1.output$LOD
    tmp[is.na(tmp)] <- 0.1
    
    RDS.tmp <- table(round(apply(tmp, 1, max, na.rm =
                                   T), digits = 1))
    RDS <- rep(0, times = 5000)
    RDS[(as.numeric(names(RDS.tmp)) * 10)] <- RDS.tmp
    RDS <- rev(cumsum(rev(RDS)))
    
    output[[4]] <- RDS
    
    output[[1]] <- which(round(((
      dim(map1.output$LOD)[1] - RDS
    ) / dim(
      map1.output$LOD
    )[1]) * q.value * log10(dim(
      map1.output$LOD
    )[1]), digits = 3) - round(output[[3]] / output[[4]], digits = 3) > 0)[1]
    names(output) <- c(
      "Significance_threshold",
      "False_discoveries",
      "False_discoveries_average",
      "Real_discoveries"
    )
  }
  
  ###Return
  return(output)
}
