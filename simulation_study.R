#### SIMULATION STUDY #####


source("sim_fun.R")  #for simulation
source("aux_fun.R")  #for plots  
source("DB_est.R")   #depth-based estimation main function

source("est_comp.R")  #Implementation of CM estimation method. NOT PROVIDED: To be downloaded from 
### https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data


###### SET PARAMETERS
m = 101
p_values = c(4,30)
t = seq(0, 1, length.out = m)
n_values= c(50,100)

swarp_values = c(0.5, 1, 5,7.5)
sdist_values = c(0,0.5,1)
serror_values = c(0,1)
out_values=0                   #No outlier contamination -> switch to 1 for outliers (see line 106 below)

n_ite= 1:100    #nb of simulation runs per setting

df= expand.grid(n_ite,out_values,serror_values, sdist_values , swarp_values, n_values,p_values)[,7:1]
names(df) <- c("p","n", "sigma_warp","sigma_dist","sigma_error","outmod","n_ite")
n_runs=nrow(df)

LISE = data.frame(df, "DB_est"=rep(0,n_runs),"CM_est"=rep(0,n_runs))
HMISE = LISE
XMISE = LISE
PMISE = LISE
times = LISE
total_it = nrow(LISE)

k=1
  
    for (p in p_values){
 
        k_psi = ifelse(p==4,1,2)
      
        for (n in n_values){ 

            for (sigma_warp in swarp_values){ 
 
                for (sigma_dist in sdist_values){ 
          
                     for (sigma_error in serror_values){ 

                         for (ite in n_ite){

                              #Data simulation
                              simdata = sim_data(n, sigma_amp=4, sigma_warp, sigma_error, sigma_dist, 
                                                 ameans=rep(100,p), p, t, m, k_psi, out=F, out_prop=0.1)

                         
                              #Depth-based estimation            
                              st <- Sys.time()
                              DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")
                              tDB <- Sys.time()-st
                                   
                              #CM estimation            
                              st <- Sys.time()
                              test = est_components(data = simdata$simdata, smooth = TRUE)
                              tCM <- Sys.time()-st
                
                              #Estimation errors
                              LISE_DB = sum((L(t)-DBest$L_hat)^2)/m
                              LISE_CM = sum((L(t)-test$L.hat)^2)/m

                              HMISE_DB = sum(sum((simdata$hmat-DBest$h_hat)^2)/m)/n
                              HMISE_CM = sum(sum((simdata$hmat-test$H.hat)^2)/m)/n
 
                              XMISE_DB = sum(sum((simdata$simdata[,3:(p+2)]-DBest$x_hat[,3:(p+2)])^2)/m)/(n*p)
                              fitted_CM = do.call("rbind",test$X.hat)
                              XMISE_CM = sum(sum((simdata$simdata[,3:(p+2)]-fitted_CM[,3:(p+2)])^2)/m)/(n*p)
                
                              psi = do.call("cbind",simdata$psi)
                              psi_CM = do.call("cbind",test$Psi.hat)
                              PMISE_DB = sum(sum((psi-DBest$psi_hat)^2)/m)/p
                              PMISE_CM = sum(sum((psi-psi_CM)^2)/m)/p
                              
                              LISE[k,8:9] = c(LISE_DB, LISE_CM)
                              HMISE[k,8:9] = c(HMISE_DB, HMISE_CM)
                              XMISE[k,8:9] = c(XMISE_DB, XMISE_CM)
                              PMISE[k,8:9] = c(PMISE_DB, PMISE_CM)
                              times[k,8:9] = c(tDB,tCM)

                              k = k+1
 
                              cli::cli_alert_info("k:{k} of {total_it} n:{n}, k_psi:{k_psi}, sigma_warp:{sigma_warp}, 
                                                  sigma_dist:{sigma_dist}, sigma_error:{sigma_error}, ite:{ite}")


               }
              save(LISE,HMISE,PMISE,XMISE,times,file="results.RData")

            }
          }
        }
      }
   }



####################### WITH OUTLIER CONTAMINATION ######################
p_values = 4
k_psi=1
n_values= 50

swarp_values = c(0.5, 1, 5,7.5)
sdist_values = c(0,0.5,1)
serror_values = 0
cont_rate = c(0,0.05,0.1)      #Outlier contamination rates

n_ite= 1:100 #nb of iterations

df= expand.grid(n_ite,cont_rate,serror_values, sdist_values , swarp_values, n_values,p_values)[,7:1]
names(df) <- c("p","n", "sigma_warp","sigma_dist","sigma_error","cont_rate","n_ite")
n_runs=nrow(df)

LISE = data.frame(df, "DB_est"=rep(0,n_runs),"CM_est"=rep(0,n_runs))
HMISE = LISE
XMISE = LISE
PMISE = LISE
times = LISE
total_it = nrow(LISE)

k=1

for (cont in cont_rate){
  
      for (sigma_warp in swarp_values){ 
        
        for (sigma_dist in sdist_values){ 
          
            for (ite in n_ite){ #Loop over nite runs
              
              simdata = sim_data(n, sigma_amp=4, sigma_warp, sigma_error=0, sigma_dist, 
                                 ameans=rep(100,p), p, t, m, k_psi, out=T, out_prop=cont)
              
              
              #Depth-based estimation            
              st <- Sys.time()
              DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")
              tDB <- Sys.time()-st
              
              #CM estimation            
              st <- Sys.time()
              test = est_components(data = simdata$simdata, smooth = TRUE)
              tCM <- Sys.time()-st
              
              #Estimation errors
              LISE_DB = sum((L(t)-DBest$L_hat)^2)/m
              LISE_CM = sum((L(t)-test$L.hat)^2)/m
              
              HMISE_DB = sum(sum((simdata$hmat-DBest$h_hat)^2)/m)/n
              HMISE_CM = sum(sum((simdata$hmat-test$H.hat)^2)/m)/n
              
              XMISE_DB = sum(sum((simdata$simdata[,3:(p+2)]-DBest$x_hat[,3:(p+2)])^2)/m)/(n*p)
              fitted_CM = do.call("rbind",test$X.hat)
              XMISE_CM = sum(sum((simdata$simdata[,3:(p+2)]-fitted_CM[,3:(p+2)])^2)/m)/(n*p)
              
              psi = do.call("cbind",simdata$psi)
              psi_CM = do.call("cbind",test$Psi.hat)
              PMISE_DB = sum(sum((psi-DBest$psi_hat)^2)/m)/p
              PMISE_CM = sum(sum((psi-psi_CM)^2)/m)/p
              
              LISE[k,8:9] = c(LISE_DB, LISE_CM)
              HMISE[k,8:9] = c(HMISE_DB, HMISE_CM)
              PMISE[k,8:9] = c(PMISE_DB, PMISE_CM)
              XMISE[k,8:9] = c(XMISE_DB, XMISE_CM)
              times[k,8:9] = c(tDB,tCM)
              
              k = k+1
              
              cli::cli_alert_info("k:{k} of {total_it} cont_rate:{cont}, sigma_warp:{sigma_warp}, 
                                                  sigma_dist:{sigma_dist}, ite:{ite}")
              
              
            }
            save(LISE,HMISE,PMISE,XMISE,times,file="results_outliers.RData")
            
          }
        }
      
    
  }



