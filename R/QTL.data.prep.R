###Data preparation for QTL mapping; start function for most of subsequent functions

###input
#trait.matrix (matrix; traits per row and strains per column)
#strain.trait (vector with strain names per trait. Should match strain.map)
#strain.map (markers in rows, genotypes in columns)
#strain.marker (markers in rows, with columns name, chromosome and position)

###output
#a list witht the entries Trait, Map, and Marker. Where the Traits are aligned with the map

###Description
#Arranges the trait and the map data so they are aligned and runs some basic checks on the input data
#removes the non-present strains it warns about

###See also
#QTL.map.1



###Data preparation
    QTL.data.prep <- function(trait.matrix,strain.trait,strain.map,strain.marker){
                              if(missing(trait.matrix)){                             stop("requires trait data, traits per row and strains per column")}
                              if(missing(strain.trait)){                             stop("requires the strain names matching the trait data")}
                              if(missing(strain.map)){                               stop("requires the genetic map of the population")}
                              if(missing(strain.marker)){                            stop("requires the marker info (name, chromosome, basepair) of the map")}
                              if(ncol(trait.matrix) != length(strain.trait)){        stop("the number of strains does not correspond with the trait matrix")}
                              if(sum(tolower(strain.trait) %in% tolower(colnames(strain.map))) < length(strain.trait)){
                                                                                     warning(paste("The strains: ",paste(strain.trait[!tolower(strain.trait) %in% tolower(colnames(strain.map))],collapse=", ")," are not in the strain.map file"))
                              }
                              if(nrow(strain.map) != nrow(strain.marker)){           stop("the number of markers (strain.map) does not equal the number of marker descriptions (strain.marker)")}
                              if(ncol(strain.marker) != 3){                          stop("make sure the strain.marker file contains a column with the marker name, chromosome, and bp location")}
                              if(nrow(trait.matrix)==1){                             small <- TRUE}
                              if(nrow(trait.matrix)!=1){                             small <- FALSE}
                              
                              strain.map <- data.matrix(strain.map)

                              ###remove strains not present
                              if(small){
                                  rowna <- rownames(trait.matrix)
                                  trait.matrix <- trait.matrix[1,which(tolower(strain.trait) %in% tolower(colnames(strain.map)))]
                                  colna <- names(trait.matrix)
                                  trait.matrix <- t(matrix(trait.matrix))
                                      colnames(trait.matrix) <- colna
                                      rownames(trait.matrix) <- rowna
                              } else {
                                  trait.matrix <- trait.matrix[,which(tolower(strain.trait) %in% tolower(colnames(strain.map)))]
                                  strain.trait <- strain.trait[tolower(strain.trait) %in% tolower(colnames(strain.map))]
                              }
                                  
                              ###Make matrix
                              map.nu <- matrix(NA,nrow=nrow(strain.map),ncol=length(strain.trait))
                              for(j in 1:length(strain.trait)){
                                  map.nu[,j] <- strain.map[,tolower(colnames(strain.map)) == tolower(strain.trait[j])]
                              }
                              colnames(map.nu) <- strain.trait
                              rownames(map.nu) <- rownames(strain.map)

                              ###add rownames to the trait matrix
                              if(length(rownames(trait.matrix)) == 0){
                                  rownames(trait.matrix) <- 1:nrow(trait.matrix)
                              }

                              ###Make output
                              output <- NULL
                              output[[1]] <- trait.matrix
                              output[[2]] <- map.nu
                              output[[3]] <- data.frame(strain.marker)
                              names(output) <- c("Trait","Map","Marker")

                              ###return output
                              return(output)
                             }