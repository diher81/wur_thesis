###Simulate simple additive QTL
#Sterken, 2017

###input
#strain.map (markers in rows, genotypes in columns)
#strain.marker (markers in rows, with columns name, chromosome and position)
#n_per_marker (number of simulated QTL per marker)
#sigma_peak (standard deviation of QTL peak; standard normal)
#sigma_error (standard deviation of noise; standard normal)

###output
#a list with: Simulated traits ([[1]]), Map ([[2]]), and Marker ([[3]]).

###Description
# Noise based on standard normal
# Fractional genotypes are transformed to integer numbers using round for simulating QTL
# Heterozygousity may be simulated using '0' or NA; these are ignored in assigning QTL effects but instead are given
#  an effect based on a uniform distribution.


QTL.sim.data.per.marker <- function(strain.map,strain.marker,n_per_marker,sigma_peak,sigma_error,replicates){
                                    if(missing(replicates)){                               replicates <- 1}
  
                                    if(sum(strain.map%%1>0,na.rm=T)>0){
                                          ###fractional loci are treated as closest genotype
                                          strain.map.use <- round(strain.map,digits = 0)
                                          warning("set fractional genotypes to closest genotype using round() for simulating QTL")
                                      } else {
                                          strain.map.use <- strain.map
                                      }
    
                                      sim.map.file <- list()
                                      
                                      sim.map.file[[2]] <- data.matrix(strain.map.use)
                                          types <- unique(as.numeric(as.character(unlist(sim.map.file[[2]]))))
                                          types <- types[!is.na(types)]
                                          if(length(types) > 2){
                                              types <- sort(types)[c(1,length(types))]
                                              warning(paste("more than 2 genotypes detected, only simulating with the most extreme numbers:",paste(types,collapse=" "),sep=" "))
                                          }
                                          sim.map.file[[2]] <- sim.map.file[[2]][rep(1:nrow(sim.map.file[[2]]),each=n_per_marker),]
                                      
                                      for(i in 1:replicates){
                                          tmp <- sim.map.file[[2]]
                                              tmp[sim.map.file[[2]]==types[1]] <- sigma_peak
                                              tmp[sim.map.file[[2]]==types[2]] <- 0
                                              ###markers of unknown genotype get filled with a random number
                                              tmp[!sim.map.file[[2]]%in%types] <- runif(sum(!sim.map.file[[2]]%in%types),0,sigma_peak)
                                              tmp <- tmp + rnorm(n=length(tmp),mean=0,sd=sigma_error)
                                              rownames(tmp) <- paste(rownames(tmp),rep(1:n_per_marker,times=nrow(sim.map.file[[2]])/n_per_marker),sep="_")
                                      
                                          if(i == 1){
                                              sim.map.file[[1]] <- tmp
                                          } else {
                                              sim.map.file[[1]] <- tmp + sim.map.file[[1]]
                                          }
                                      }
                                      ###average
                                          sim.map.file[[1]] <- sim.map.file[[1]] / replicates
                                          
                                      ###paste original strain map back
                                      sim.map.file[[2]] <- strain.map
                                      
                                      sim.map.file[[3]] <- strain.marker
                                      
                                      return(sim.map.file)
                                     }    
