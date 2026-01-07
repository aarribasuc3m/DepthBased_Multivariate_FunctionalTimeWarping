##### Auxiliar functions for Depth-based estimation in the Latent Deformation Model

# Monotone transformation of an individual curve (as in Dupuy, Loubes and Maza, Statistics and Computing 2011)
monot <- function(x) { 
  #x is a vector representing a single realization of a curve
  m <- length(x)
  x.mon <- x
  if(m>1){
    for (k in 2:m){
      x.mon[k]<-x.mon[k-1]+abs(x[k]-x[k-1])
    }
  }
  return(x.mon) 
}

# Moving average smoothing of a sample of curves with adaptive window length
moving_average <- function(xmat, order,thresh=10^-5){
  cmx = colMeans(xmat)
  xmat_ma = t(apply(xmat,1,forecast::ma,order=order))
  cmx_new = colMeans(xmat_ma)
  
  m=ncol(xmat)
  
  k = round(order/2)
  i1= k+1
  i2= m-k
  msre= mean( ((cmx-cmx_new)/cmx)[i1:i2]^2)
  
  if ( msre > thresh){
    or=order-1
    while(msre > thresh){
      xmat_ma = t(apply(xmat,1,forecast::ma,order=or))
      cmx_new = colMeans(xmat_ma)
      k = round(or/2)
      i1= k+1
      i2= m-k
      msre= mean( ((cmx-cmx_new)/cmx)[i1:i2]^2)

      or=or-1
    }
    or = or +1
  }else{
    or=order+1
    while(msre< thresh & or<m/3){
      xmat_ma = t(apply(xmat,1,forecast::ma,order=or))
      cmx_new = colMeans(xmat_ma)
      k = round(or/2)
      i1= k+1
      i2= m-k
      msre= mean( ((cmx-cmx_new)/cmx)[i1:i2]^2)

      or=or+1
    }
    or = or -1
  }
  
  or = ifelse(or==1,order,or)     #in case no smoothing order preserves the curves well enough, we keep the original one (m/5)
  xmat_ma = t(apply(xmat,1,mov_avg,order=or))
  
  return(xmat_ma)
}
mov_avg <- function(x,order){
  
  x_ma <- forecast::ma(x,order=order)
  k = round(order/2)
  or=k
  for (j in 1:k){ #we fill in values in the extremes
    j1=k-j+1
    j2=m-k+j
    x_ma[j1] = mean(c(x[1:j1],x_ma[(j1+1):(j1+or)]))
    x_ma[j2] = mean(c(x_ma[(j2-or):(j2-1)],x[j2:m])) 
    or=or-1
  }
  x_ma[1]=x[1]
  x_ma[m]=x[m]
  return(x_ma)
}

# Monotone transformation and normalization of the multivariate functions, possibly including smoothing
MonTrans <- function(x,smooth=TRUE, method="MA", bw=0.01) { 
             # x is a long-format data frame with a 1st column id, a 2nd column t, and an
             #additional column for each new variable 
  
  n <- length(unique(x$id)); m <- length(unique(x$t)) # n: number of observations (curves);  m: length of time grid
  p <- ncol(x)-2 #p: dimension of the data
  
  # if (smooth==TRUE & method=="CM"){
  #   x <- smooth_curves(x,p,bw=bw)
  # }
  if (smooth==TRUE & method=="MA"){
    for (j in 1:p){
      xmat = matrix(x[,j+2],nrow=n,ncol=m, byrow=T)
      xm = moving_average(xmat,order=round(m/5))
      x[,j+2] = as.numeric(t(xm))
    }
  }
  
  mon.data.df <- x        #monotone multivariate data set
  mon.data.nor.df <- x    #monotone normalized multivariate data set
  data.nor.df <- x        #normalized multivariate data set
  x <- as.matrix(x) 
  
  amp <- matrix(0, nrow=n, ncol=p)
  for (j in 1:p) {
     y = matrix(x[,j+2], nrow=n, ncol=m, byrow=T)                      
     amp[,j] <- abs(apply(y,1,max))                                       
     amp_j_mat <- matrix(rep(amp[,j],each=m),nrow=n, byrow=T)
     y_nor <- y/amp_j_mat
     data.nor.df[,j+2] <- as.numeric(t(y_nor))
    
    # Monotone transformation of the curves (as in Dupuy, Loubes and Maza, Statistics and Computing 2011)
    x.mon<- matrix(,n,m)
    for (i in 1:n) {
       x.mon[i,] <- monot(y_nor[i,])
    }
    
    min <- apply(x.mon,1,min)             #  We have to first monotonize, then normalize again
    min_j_mat <- matrix(rep(min,each=m),nrow=n, byrow=T)   # Otherwise the final curves might not be normalized. 
    
    amp2 <- apply(x.mon- min_j_mat,1,max)      
    amp_j_mat <- matrix(rep(amp2,each=m),nrow=n, byrow=T)   # Otherwise the final curves might not be normalized. 
    xmon_nor <- (x.mon- min_j_mat)/amp_j_mat +  min_j_mat

    mon.data.nor.df[,j+2] <- as.numeric(t(xmon_nor))
    mon.data.df[,j+2] <- as.numeric(t(x.mon))

  }
  
  return(list(mon.data.df = mon.data.df, mon.data.nor.df = mon.data.nor.df, data.nor.df = data.nor.df, amp=amp, smooth=data.frame(x)))

}

# Deepest curve in the pooled monotone sample (or trimmed mean of deepest curves)
median_overall_monot <- function(xmon,x,trim=0){ 
  # x is a long-format data frame with a 1st column id, a 2nd column t, and an
  #additional column for each new variable 
  # xmon contains the monotone sample
  
  n <- length(unique(x$id)); m <- length(unique(x$t)) # n: number of observations (curves);  m: length of time grid
  p <- ncol(x)-2 #p: dimension of the data
  t <- x$t[1:m]
  
  y <- matrix(as.matrix(x[,3:(p+2)]),nrow=n*p,ncol=m,byrow=T)
  ymon <- matrix(as.matrix(xmon[,3:(p+2)]),nrow=n*p,ncol=m,byrow=T)
  
  y_fd <- roahd::fData(t, ymon)
  mbd <- roahd::MBD(y_fd)
  # med <- which.max(mbd)
  med <- sort(mbd,decreasing = T, index.return=T)$ix[1:(trim*n+1)]
  if (length(med)>1){
    median_curve <- colMeans(y[med,])        #trimmed mean of deepest curves
  }else{
    median_curve <- y[med,]
  }
  
  return(list(median_curve=median_curve,med=med))
}

# Deepest curve in each component of the monotone sample
median_component_wise_monot<- function(x){ 
  # x is a long-format data frame with a 1st column id, a 2nd column t, and an
  #additional column for each new variable 
  
  n <- length(unique(x$id)); m <- length(unique(x$t)) # n: number of observations (curves);  m: length of time grid
  p <- ncol(x)-2 #p: dimension of the data
  t <- x$t[1:m]
  x <- as.matrix(x) # Data matrix
  
  med <- vector("numeric", length=p)
  
  for (j in 1:p) {
    y = matrix(x[,j+2], nrow=n, ncol=m, byrow=T)
    y_fd <- roahd::fData(t, y)
    mbd <- roahd::MBD(y_fd)
    med[j] <- which.max(mbd)
  }
  return(med)
}

# Component amplitude function estimates as the deepest curves in each monotone sample
comp_est <- function(D){ #D is the output of MonTrans
  
  median_comp <- median_component_wise_monot(D$mon.data.nor.df) #Based on marginal component-wise medians
  m <- length(unique(D$mon.data.nor.df$t))
  index <- (median_comp -1)*m +1
  medians_nor <- c()
  mon_medians_nor <- c()
  for (j in 1:length(index)){
    medians_nor = cbind(medians_nor, D$data.nor.df[index[j]:(index[j]+m-1), j+2])
    mon_medians_nor = cbind(mon_medians_nor, D$mon.data.nor.df[index[j]:(index[j]+m-1), j+2])
  }
  
  return(list(comp=medians_nor,comp_mon=mon_medians_nor,index_med=median_comp))
}

function_composition <- function(x,y,t){ #return y^-1(x(t)), y is a strictly monotone function
  #y is a vector of length m
  #x is a vector of length m
  #t is a vector of length m, with the time domain of x and y
  
  m = length(x)
  mi = min(y)
  ma = max(y)
   
  xout = mi+(ma-mi)*(x-min(x))/(max(x)-min(x))  #adjust to ensure xout and y have the same range 
  
  idd = seq(1,m,20)   #low resolution for linear interpolation to avoid problems at "almost plateaux"
  tt = c(min(t),t[idd],max(t) ) 
  
  vec = c(min(xout), xout[idd], max(xout) )  #low resolution of the values in xout
  xx = unique(vec )
  idx=which(!duplicated(vec)) #which indexes in vec have been retained
  tt0 = tt[idx]
  
  zz=approx(y, t , xx, rule=2)$y  #low_res estimate of y^-1(x(t)) at tt0 by linear interpolation
  
  z <- demography::cm.spline(tt0,zz,n=m)$y  #cubic spline monotonic interpolation of  y^-1(x(t)) at t
  z[c(1,m)] <- t[c(1,m)] # We force extremes to be tied to the starting and end points of the interval 
  
  return(z)
}

# Individual warping function estimates
warp_est <- function(D, comp_mon){ #D is the output of MonTrans; comp_mon is a matrix with p-columns and the deepest monotone function in each component
  
  x <- D$mon.data.nor.df   #monotone sample
  
  n <- length(unique(x$id)); m <- length(unique(x$t)) # n: number of observations (curves);  m: length of time grid
  p <- ncol(x)-2 #p: dimension of the data
  t <- x$t[1:m]
  
  h = list() 
  
  for (j in 1:p){ #In each component, estimate h_ij = gamma_j^-1 (x_ij(t)), where gamma_j is comp
    
    y = comp_mon[,j]    #monotone version of the component target j   

    h[[j]] <- matrix(0, nrow=n, ncol=m)
    
    for(i in 1:n){    
      
      xout = x[((i-1)*m+1):(i*m),j+2]  # monotone curve x_ij 
 
      h[[j]][i,] <- function_composition(xout,y,t)
    }
  }
  
  h_mean = apply(simplify2array(h), 1:2, mean) #mean over p components
  
  return(list(h_list=h, h_mean=h_mean))
}

# Component-distortion estimates
psi_est <- function(comp_mon, L_hat, t){
  #comp_mon is a matrix with p-columns containing the deepest monotone function in each component
  #L_hat is a vector containing the amplitude function estimate
  
  m <- length(L_hat)  #m: length of time grid
  p <- ncol(comp_mon) #p: dimension of the data
  
  psi = matrix(0, nrow=m, ncol=p)
  
  L_mon = monot(L_hat)

  for (j in 1:p){  #In each component, estimate psi_j = lambda^-1 (gamma_j(t)), where gamma_j is comp
    xout = comp_mon[,j]

    psi[,j] <-function_composition(xout,L_mon,t)

  }
  
  return(psi)
}

# Fitted curves
fitted_curves <- function(comp_hat, h_hat, t, amat){
  #comp_hat is a matrix with p-columns containing the estimate of each component amplitude function
  #h_hat is an n*m matrix with the individual warping function estimates 
  #t is a vector of length m, with the observation time poins
  #amat is an n*p matrix with the estimates of the a_ij amplitude factors
  
  n <- nrow(h_hat); m <- ncol(h_hat) # n: number of observations (curves);  m: length of time grid
  p <- ncol(comp_hat) #p: dimension of the data
  
  x_fitted = matrix(0, nrow = n*m, ncol = p+2)
  
  x_fitted[,1] = rep(1:n,each = m)
  x_fitted[,2] = rep(t, times = n)
  
  for(j in 1:p){
    complist = list()
    for(i in 1:n){
      fun = approx(t,comp_hat[,j],h_hat[i,],rule=2)$y  #composition of gamma_j with h_ij
      complist[[i]] =  amat[i, j] * fun
    }
    x_fitted[,j+2] = unlist(complist) 
  }
  
  x_fitted = as.data.frame(x_fitted)
  names(x_fitted)[1:2] = c("id", "t")
  
  return(x_fitted)
}

# Registered curves
registered_curves <- function(x, h_hat, t){
  # x is a long-format data frame with a 1st column id, a 2nd column t, and an
  #additional column for each new variable 
  #h_hat is an n*m matrix with the individual warping function estimates 
  #t is a vector of length m, with the observation time poins
  
  n <- length(unique(x$id)); m <- length(unique(x$t)) # n: number of observations (curves);  m: length of time grid
  p <- ncol(x)-2 #p: dimension of the data
  
  x_registered = matrix(0, nrow = n*m, ncol = p+2)
  
  x_registered[,1] = rep(1:n,each = m)
  x_registered[,2] = rep(t, times = n)
  
  for(j in 1:p){
    y = matrix(x[,j+2], nrow=n, ncol=m, byrow=T)
    complist = list()
  
      for(i in 1:n){  

        h_inv = approx(h_hat[i,], t, t,rule=2)$y    #h_ij^-1
        fun = approx(t,y[i,],h_inv,rule=2)$y        #x_ij(h_ij^-1)
        complist[[i]] = fun
        
     }
    x_registered[,j+2] = unlist(complist) 
  }
  
  x_registered = as.data.frame(x_registered)
  names(x_registered) = names(x)
  
  return(x_registered)
}


