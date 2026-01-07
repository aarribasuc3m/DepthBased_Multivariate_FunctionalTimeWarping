### Functions used to generate data in the simulation study
### Whenever a reference to Carroll&Muller,Biometrics, 2023 is made, the original code is available at
### https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data

# normalized amplitude function. As in Carroll&Muller 2023
L = function(t){
  Lshape = 20-(5)*cos(4*pi*t)+(3)*sin(pi*t^2)+15*(t)^2 
  (20-(5)*cos(4*pi*t)+(3)*sin(pi*t^2)+15*(t)^2) / max(Lshape)
}

# normalized amplitude function for outliers
L_out1 = function(t){
  Lshape_out = exp(-25*(t-0.5)^2)+t  
  Lshape_out/ max(Lshape_out) *diff(range(L(t)))+min(L(t))
}

# deterministic component based distortions - For k_psi=1 and p=4. As in Carroll&Muller 2023 
psi_function_CM <- function(p,t){
  
  psi = list() 
  psiinv = list()
  
  if(p==4){
  alpha = c(1,2)
  beta = c(1/2,2)
  }else{
    alpha = c(1,2,1/2,3/4,3/2)
    beta = c(1/2,2,1,3/2,2)
  }
  
  lambda = .5
  
  for(j in 1:(p/2)){
    t = seq(0, 1, length.out = m)
    psi[[j]] = lambda*pbeta(t, alpha[j], beta[j]) + (1-lambda)*t
    psiinv[[j]] = approx(psi[[j]], t, t)$y
    
    psiinv[[j+(p/2)]] = 2*t - psiinv[[j]]
    psi[[j+(p/2)]] = approx(psiinv[[j+(p/2)]], t, t)$y
    
  }
  return(list(psi_true = psi,psi_inv = psiinv))
}

# random warping function generation as proposed in Dupuy, Loubes, Maza, Statistics and Computing, 2009
warps_DLM<- function(t,n, eps=0.005,N=25){  
  
  p = length(t)
  
  #Initial warping function: identity
  h <- matrix(rep(t,n),nrow=n,ncol=p,byrow=T)
  
  #Iterative procedure to generate warping functions
  for (i in 1:N) {
    u<-runif(1)*(1-20*eps)+10*eps
    v<-runif(n)*2*eps+(u-eps)
    pu<-floor(u*(p-1))+1
    
    
    for (j in 1:n){
      w<-v[j]*t[1:pu]/u
      w<-append(w, ((1-v[j])/(1-u))*t[(pu+1):p]+(v[j]-u)/(1-u) )
      
      h[j,]<-approx(t, w, h[j,], method="linear",rule=2)$y   ##composition of w and h in t
      
    }
  }
  
  return(h)  
}

# function to generate curve with individual warping + component-based warping + aditional distorsion. As in Carroll&Muller 2023
makecurve = function(i, j, sigma_dist, hmat, psi_true, t, fun=L){
  if(sigma_dist==0){
    dist=t
  }else{
    d = rnorm(n = 1, 0, sigma_dist^2)
    dinv = (exp(t*d)-1)/(exp(d)-1)
    dist =  approx(dinv, t, t)$y
  }
  hd = approx(t,  hmat[i,], xout = dist)$y
  PoH = approx(t, psi_true[[j]], xout = hd)$y
  fun(PoH)
}


 
#### Main function for data simulation
sim_data = function(n, sigma_amp, sigma_warp, sigma_error = 0, sigma_dist = 0, ameans, p, t, m, k_psi, out=F, out_prop=0.05){
  
  set.seed(Sys.time()) #porque hay alguna función a la que se llama después que parece fijar la semilla y hacía que todas las simulaciones salieran iguales

  #Individual warping functions generation
  if (sigma_warp%in% c(0.5,1)){
     hinvmat = matrix(0, nrow = n, ncol = m)
     hmat = matrix(0, nrow = n, ncol = m)
     z = rnorm(n = n, 0, sigma_warp^2)
    for(i in 1:n){
       if (sigma_warp==0){
         hinvmat[i,] = t
       }else{
          hinvmat[i,] = (exp(t*z[i])-1)/(exp(z[i])-1)
       }
          hmat[i,] =  approx(hinvmat[i,], t, t)$y
    }
  }else{
    hmat <- warps_DLM(t,n,eps=sigma_warp*10^-3, N=2500)   #we get H st E[H]=id
    hinvmat = matrix(0, nrow = n, ncol = m)
    for (i in 1:n){
       hinvmat[i,] =  approx(hmat[i,], t, t,rule = 2)$y
    }
    hmat=hinvmat       #since we want E[H^1]=id
  }
  
  #Random amplitudes generation
  amat = matrix(0, nrow = n , ncol = p)
  for(j in 1:p){
    amat[,j] = pmax(rnorm(n, mean = ameans[j], sd = sigma_amp^2), 5)
  }
  
  #Component-based distortions generation
  if (k_psi==1){
    psi_true = psi_function_CM(p,t)$psi_true   #Deterministic functions
  }else{
     psi_inv <- warps_DLM(t,p,eps=0.005, N=2500)  # Random component based distortions such that  E[psi_inv]=id
     psi = matrix(0, nrow = p, ncol = m)
      for (i in 1:p){
         psi[i,] =  approx(psi_inv[i,], t, t,rule = 2)$y
      }
    psi_true <- lapply(seq_len(nrow(psi)), function(i) psi[i,])  #list from matrix rows
  }
  
  #Component amplitude functions
  makecurvenowarp = function(j){
    PoI = approx(t, psi_true[[j]], xout = t)$y
    L(PoI)
  }

  ###### If outlier contamination  
  n_out=round(out_prop*n)
  
  fL <- function(i){  #function to create shape outliers.
    if (out==T & i> n-n_out){
        return(L_out1)
    }else{
      return(L)
    }
  }
  
  if(out==T){
    outs=cbind(rep((n-n_out+1):n,p), rep(1:p, each=n_out)) #indicator of outlying curve
  }else{
    outs <- c()
    }
  #####
  
  #Multivariate curves generation
  simdata = matrix(0, nrow = n*m, ncol = p+2)
  simdata[,1] = rep(1:n,each = m)
  simdata[,2] = rep(t, times = n)
  
  for(j in 1:p){
    complist = list()
    for(i in 1:n){
      complist[[i]] =  amat[i, j] * makecurve(i, j, sigma_dist, hmat, psi_true, t, fun=fL(i)) + 
        rnorm(n = m, mean = 0,  sd = sigma_error)
      complist[[i]] = pmax(complist[[i]], .1) #force positive; avoid measurement error leading to negatives
    }
    simdata[,j+2] = unlist(complist) 
  }
  
  simdata = as.data.frame(simdata)
  names(simdata)[1:2] = c("id", "t")
  
  return(list(simdata = simdata, hmat = hmat, amat = amat, compfun=makecurvenowarp, psi=psi_true, outs=outs))
}
