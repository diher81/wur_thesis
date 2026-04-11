###rewrites the genetic map to a dataframe plottable in ggplot


###input
#strain.map (of the investigated strains: markers in rows, genotypes in columns)
#strain.marker (of the investigated strains: markers in rows, with columns name, chromosome and position)

###output
#dataframe with columns: strain, chromosome, position_start, position_end, and genotype

###Description

prep.ggplot.genetic.map <- function(strain.map, strain.marker){
          if(missing(strain.map)){                        
            stop("Requires genetic map with markers per row and strains per column")}
          if(missing(strain.marker)){                     
            stop("Requires marker file with 3 columns: name, chromosome, position")}
          if(nrow(strain.map) != nrow(strain.marker)){    
            stop("map does not match markers")}
          if(ncol(strain.marker) != 3){                   
            stop("Requires marker file with 3 columns: name, chromosome, position")}

  # strain.map: input matrix
  # 2: operate column-wise (i.e., per strain).
  # diff: computes differences between consecutive rows.
  # TRUE → genotype changes between adjacent markers
  # FALSE → genotype remains the same
          diffmap <- apply(strain.map, 2, diff) != 0
  # Returns TRUE when the chromosome changes between rows
          chrmap <- diff(as.numeric(as.factor(strain.marker[,2]))) != 0

          tmp <- as.data.frame(cbind(name = as.character(unlist(strain.marker[,1])),
                                 chromosome = as.character(unlist(strain.marker[,2])),
                                 position = as.numeric(as.character(unlist(strain.marker[,3]))),
                                 strain = rep(colnames(strain.map), each = nrow(strain.map)),
                                 genotype = as.numeric(strain.map),
                                 # Genotype change indicators
                                 diff_left = as.numeric(rbind(diffmap, FALSE)),
                                 diff_right = as.numeric(rbind(FALSE, diffmap)),
                                 # Chromosome boundary indicators
                                 chr_left = as.numeric(c(chrmap, TRUE)),
                                 chr_right = as.numeric(c(TRUE, chrmap))))
          for(i in c(3,5:9)){
              tmp[,i] <- as.numeric(as.character(unlist(tmp[,i])))
          }

          tmp <- cbind(tmp[,1:5], diffs = apply(tmp[,6:9],1,sum))
          tmp <- tmp[tmp[,6] != 0,]

          tmp <- cbind(tmp[-nrow(tmp),-6],tmp[-1,-6])
          tmp <- tmp[tmp[,2] == tmp[,7] & tmp[,4] == tmp[,9] & tmp[,5] == tmp[,10],]

          output <- as.data.frame(cbind(strain = as.character(unlist(tmp[,4])),
                                        chromosome = as.character(unlist(tmp[,2])),
                                        position_start = as.numeric(as.character(unlist(tmp[,3]))),
                                        position_end = as.numeric(as.character(unlist(tmp[,8]))),
                                        genotype = as.numeric(as.character(unlist(tmp[,5])))))
              output[,1] <- as.character(unlist(output[,1]))
              output[,2] <- as.character(unlist(output[,2]))
              output[,3] <- as.numeric(as.character(unlist(output[,3])))
              output[,4] <- as.numeric(as.character(unlist(output[,4])))

          return(output)
                                   }