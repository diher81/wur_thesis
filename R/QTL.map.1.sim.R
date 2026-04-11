###Conduct mapping simulation to investigate power
#Sterken, 2017

###input
#strain.map (markers in rows, genotypes in columns)
#strain.marker (markers in rows, with columns name, chromosome and position)
#n_per_marker (number of simulated QTL per marker)
#sigma_peak (standard deviation of QTL peak; standard normal)
#sigma_error (standard deviation of noise; standard normal)

###output
#a data frame with per rows the outcomes of the simulated peak-sizes. The number of true, false, and undetected
#  QTL; quantiles of the QTL effect-size estimation; quantiles of the QTL location estimation

###Description
# Noise based on standard normal
# Fractional genotypes are transformed to integer numbers using round for simulating QTL
# Heterogeneity may be simulated using '0' or NA; these are ignored in assigning QTL effects but instead are given
#  an effect based on a uniform distribution.

###Updates
# QTL detection was very strict (peak at location of peak, changed that to a peak within a percentage of the genetic map (2.5%))


QTL.map.1.sim <- function(strain.map,strain.marker,n_per_marker,sigma_peak,sigma_error,threshold,peak_within_perc){
                            if(missing(strain.map)){                                        stop("specify strain.map (genetic map)")}
                            if(missing(strain.marker)){                                     strain.marker <- data.frame(cbind(name=1:dim(strain.map)[1],chromosome=NA,bp=1:dim(strain.map)[1]))
                                                                                            rownames(strain.marker) <- 1:dim(strain.map)[1]}
                            if(missing(n_per_marker)){                                      n_per_marker <-10}
                            if(missing(sigma_peak)){                                        sigma_peak <- c(1,1.15,1.3,1.47,1.64,1.82,2,2.2,2.47,2.73,3.05,3.45,4)}
                            if(missing(sigma_error)){                                       sigma_error <- 1}
                            if(missing(threshold)){                                         threshold <- c(3.5)}
                            if(missing(peak_within_perc)){                                  peak_within_perc = 0.0125}
                            ###1 sigma difference on 1 sigma noise is 20% of variance
                            ###2 sigma difference on 1 sigma noise is 50% of variance
                            ###3 sigma difference on 1 sigma noise is 69% of variance
                            ###4 sigma difference on 1 sigma noise is 80% of variance; summary(lm(c(rnorm(10000,0,1),rnorm(10000,4,1))~rep(c(-1,1),each=10000)))
                            ###cbind(var_expl=seq(5,95,by=5),sigma=c(0.46,0.67,0.84,1,1.15,1.3,1.47,1.64,1.82,2,2.2,2.47,2.73,3.05,3.45,4,4.75,6.0,8.75))
    
                            which.max.rev <- function(x){length(x)-which.max(rev(x))+1}
                            
                            out <- NULL
                            
                            pb <- txtProgressBar(min = 0, max = length(sigma_peak), style = 3)
                            counter <- 0
                            
                            for(i in sigma_peak){
                              
                              # progress insight
                              counter <- counter + 1
                              setTxtProgressBar(pb, counter)
                              start_time <- Sys.time()
                              message(sprintf("Start sigma_peak = %s (%d/%d)", 
                                              i, 
                                              which(sigma_peak == i), 
                                              length(sigma_peak)))
                              
                                ###Generate peaked data with random error
                                    sim.map.file <- QTL.sim.data.per.marker(strain.map,strain.marker,n_per_marker,sigma_peak=i,sigma_error)
                                    sim.emp.file <- QTL.sim.data.per.marker(strain.map,strain.marker,n_per_marker,sigma_peak=0,sigma_error)
    
                                ###Map it
                                    sim.map.file <- QTL.map.1(sim.map.file[[1]],sim.map.file[[2]],sim.map.file[[3]])
                                    
                                    ###locate peak
                                    loc_peak1 <- as.numeric(as.character(unlist(apply(sim.map.file$LOD,1,which.max))))
                                    loc_peak2 <- as.numeric(as.character(unlist(apply(sim.map.file$LOD,1,which.max.rev))))                                    
                                    loc_peak <- floor(apply(cbind(loc_peak1,loc_peak2),1,median))
                                    
                                    ###effects and lod
                                    lod_peak <- as.numeric(t(sim.map.file$LOD))[loc_peak+seq(0,((nrow(sim.map.file$Map)^2)*n_per_marker)-nrow(sim.map.file$Map),by=nrow(sim.map.file$Map))]
                                    eff_peak <- as.numeric(t(sim.map.file$Effect))[loc_peak+seq(0,((nrow(sim.map.file$Map)^2)*n_per_marker)-nrow(sim.map.file$Map),by=nrow(sim.map.file$Map))]
                                                                        
                                    ###effects and location of simulated peak
                                    loc_sim <- rep(1:nrow(sim.map.file$Map),each=n_per_marker)

                                    ###peak is true if: above threshold and within % of total genetic map
                                    ###peak is false if: above threshold and outside % of total genetic map
                                    ###peak is not detected if: no significant peak
                                    out.sim <- NULL
                                        out.sim <- c(out.sim,detected_true=sum(lod_peak>=threshold & abs(loc_peak-loc_sim) < floor(nrow(sim.map.file$Map)*peak_within_perc)))
                                        out.sim <- c(out.sim,detected_false=sum(lod_peak>=threshold & abs(loc_peak-loc_sim) >= floor(nrow(sim.map.file$Map)*peak_within_perc)))
                                        out.sim <- c(out.sim,detected_not = sum(lod_peak<threshold))
                                            selc <- lod_peak>=threshold & abs(loc_peak-loc_sim) < floor(nrow(sim.map.file$Map)*peak_within_perc)
                                        out.sim <- c(out.sim,detected_true_eff_est=quantile(abs(eff_peak[selc])/(i/2),probs=c(0.05,0.25,0.5,0.75,0.95)))
                                        out.sim <- c(out.sim,detected_true_loc_est=quantile(abs(as.numeric(as.character(unlist(sim.map.file$Marker[loc_sim,3])))-as.numeric(as.character(unlist(sim.map.file$Marker[loc_peak+ceiling((loc_peak2-loc_peak)/2),3]))))[selc],probs=c(0.05,0.25,0.5,0.75,0.95)))
    
                                ###Map False file
                                    sim.map.file <- QTL.map.1(sim.emp.file[[1]],sim.emp.file[[2]],sim.emp.file[[3]])
    
                                    loc_peak <- as.numeric(as.character(unlist(apply(sim.map.file$LOD,1,which.max))))
                                    lod_peak <- as.numeric(as.character(unlist(apply(sim.map.file$LOD,1,max,na.rm=TRUE))))
                                    eff_peak <- t(sim.map.file$Effect)[(ncol(sim.map.file$Effect)*(0:(nrow(sim.map.file$Effect)-1)))+loc_peak]
    
                                    loc_sim <- rep(1:nrow(sim.map.file$Map),each=n_per_marker)

                                    out.emp <- NULL
                                        out.emp <- c(out.emp,detected_true=0)
                                        out.emp <- c(out.emp,detected_false=sum(lod_peak >= threshold))
                                        out.emp <- c(out.emp,detected_not = sum(lod_peak < threshold))
    
                                        out.emp <- c(out.emp,detected_true_eff_est=quantile(rep(NA,times=100),na.rm=T,probs=c(0.05,0.25,0.5,0.75,0.95)))
                                        out.emp <- c(out.emp,detected_true_loc_est=quantile(rep(NA,times=100),na.rm=T,probs=c(0.05,0.25,0.5,0.75,0.95)))
                                ###Combine and summarize
                                    out.sim <- data.frame(cbind(name=c(paste("peak_sim_",i,sep=""),paste("noise_sim_",sigma_error,sep="")),
                                                                var_explained=c(round(summary(lm(c(rnorm(100000,0,sigma_error),rnorm(100000,i,sigma_error))~rep(c(-1,1),each=100000)))$adj.r.squared,digits=2),round(summary(lm(c(rnorm(100000,0,sigma_error),rnorm(100000,0,sigma_error))~rep(c(-1,1),each=100000)))$adj.r.squared,digits=2)),
                                                                rbind(out.sim,out.emp)))
    
                                    out <- rbind(out,out.sim)
                                    
                                    # progress insight
                                    message(sprintf("Finished sigma_peak = %s | elapsed: %.2f sec",
                                                    i,
                                                    as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
                            }
                            
                            ###Make sure numbers are numeric and characters are characteristic                            
                            for(i in 2:ncol(out)){
                                out[,i] <- as.numeric(as.character(unlist(out[,i])))
                            }
                            out[,1] <- as.character(unlist(out[,1]))


                            ###summarize all noise simulations
                            out[max(grep("noise_sim",out[,1])),-1] <- apply(out[grep("noise_sim",out[,1]),-1],2,sum)
                            out <- out[c(grep("peak_sim",out[,1]),max(grep("noise_sim",out[,1]))),]
                            
                            ##Add percentage detected & percentage false
                            output <- cbind(out,
                                            percentage_detected=round(out$detected_true/(out$detected_true+out$detected_false+out$detected_not),digits=3)*100,
                                            percentage_false=round(out$detected_false/(out$detected_true+out$detected_false+out$detected_not),digits=3)*100,
                                            type=unlist(strsplit(out$name,split="_"))[seq(1,nrow(out)*3,by=3)],
                                            sigma=unlist(strsplit(out$name,split="_"))[seq(3,nrow(out)*3,by=3)])
                                                            
                            output <- output[,c(18,19,2,3:5,16,17,6:15)]
                            colnames(output)[grepl("_est.",colnames(output))] <- paste("quantile_",unlist(strsplit(colnames(output)[grepl("_est.",colnames(output))],split="_est."))[seq(2,2*sum(grepl("_est.",colnames(output))),by=2)],sep="")
                            colnames(output) <- gsub("\\.","",colnames(output))
                            
                            colnames(output) <- paste(c(rep("Simulation",times=3),rep("QTL_detection",times=5),rep("Effect_size_estimation_(quantiles)",times=5),rep("QTL-true_peak_location_distance_(quantiles)",times=5)),colnames(output),sep=".")
                            

                            return(output)
}
