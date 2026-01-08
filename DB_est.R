##### Depth-based estimation in the Latent Deformation Model
source("DB_est_functions.R")   #depth-based estimation auxiliar functions

## Depth-based estimation in the latent deformation model
DB_est <- function(data, smooth=TRUE, method = "MA"){ 
          #data is a long-format data frame with a 1st column id, a 2nd column t, and an
          #additional column for each new variable
  
  #obtain strictly monotone curves
  D <- MonTrans(data,smooth=TRUE, method="MA") #smoothing with moving average
  
  #Lambda estimate as the deepest curve in the monotone global sample
  L_med <- median_overall_monot(D$mon.data.nor.df,D$data.nor.df)
  L_hat <- L_med$median_curve
  
  #Component amplitude function estimates as the deepest curves in each monotone sample
  comp <- comp_est(D)
  comp_hat <- comp$comp
  
  #Individual warping function estimates
  h <- warp_est(D, comp$comp_mon)
  h_hat_list <- h$h_list
  h_hat <- h$h_mean
  
  #Component-distortion estimates
  psi_hat <- psi_est(comp$comp_mon, L_hat, t)
  
  #Fitted curves
  x_hat <- fitted_curves(comp_hat, h_hat, t, D$amp)
  
  #Registered curves
  x_reg <- registered_curves(data, h_hat, t)


return(list(D=D,L_med=L_med, L_hat=L_hat, comp = comp, h_hat_list = h_hat_list, h_hat = h_hat, 
             psi_hat = psi_hat, x_hat = x_hat, x_reg = x_reg, a=colMeans(D$amp)))

}
