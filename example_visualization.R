### Example of depth-based estimation in the latent deformation model

library(gridExtra)
source("sim_fun.R")  #for simulation
source("aux_fun.R")  #for plots  
source("DB_est.R")   #depth-based estimation main function


## 1. Data simulation
m = 101
p=4
t = seq(0, 1, length.out = m)
n=50

ameans=rep(100,p)


simdata = sim_data(n, sigma_amp=4, sigma_warp=1, sigma_error=0, sigma_dist=0, ameans=ameans, p=p, t, m, k_psi=1, 
                   out=F, out_prop=0.05)

p1 <- plot_ggplot_target(t,L, main="Target function",yl=expression(lambda~(t)),lt="solid")
p2 <- psi_plot(t, matrix(unlist(simdata$psi),ncol=4) ,ggcolors, title="Component-based distorsions",sim=F)
p3 <- hplot(t,simdata$hmat, title="Individual warping functions",sim=F)
p4 <- plot_ggplot_funct(simdata$compfun, a=c(1,1,1,1), m, p, t, main="Component target curves",
                        names_comp=c("t", "gamma[1]","gamma[2]","gamma[3]","gamma[4]"))
p5 <- plot_ggplot_fun_df(simdata$simdata, simdata$compfun, a=ameans, m, p, t, 
                         main="Observed curves",names_comp=c("t", "X[1]","X[2]","X[3]","X[4]"))

grid.arrange(
  grobs = list(p1,p2,p3,p4,p5),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 4, 5),
                        c(2, 4, 5),
                        c(3, 4, 5))
)

## 2. Depth-based estimation

DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

est <- data.frame(time=t,value=DBest$L_hat,est=c(rep("depth",m)))

p1<- plot_ggplot_target(t,L, main="Target (--) + estimate (-)",yl=expression(hat(lambda)~(t)),
                        col="black",sz=0.3, lt="dashed", tsz=10) + 
  geom_line(data=est, aes(x=time,y=value,group=est,colour = est, linetype = est), size=1) + 
  scale_linetype_manual(values=1) +
  scale_colour_manual(values="darkgray") 

p2 <- psi_plot(t,DBest$psi_hat,ggcolors, title="Component-based dist.(--) + estimates (-)",yl=bquote(hat(psi)[j]~(t)),
               sim=F, true_psi=simdata$psi)
p3 <- hplot(t,DBest$h_hat, title="Ind. warping functions estimates",yl=bquote(hat(h)[i]~(t)),sim=F)

func_comp_est<-function(j){DBest$comp$comp[,j]}

p4 <- plot_ggplot_funct(func_comp_est, a=c(1,1,1,1), m, p, t, main="Comp. target curves (--) + estimates (-)",
                        names_comp=c("t", "hat(gamma)[1]","hat(gamma)[2]","hat(gamma)[3]","hat(gamma)[4]"),
                        true_fun=simdata$compfun, tsz=10)
p5 <- plot_ggplot_fun_df(DBest$x_reg, func_comp_est, a=DBest$a, m, p, t, 
                         main="Registered curves",names_comp=c("t", "Y[1]","Y[2]","Y[3]","Y[4]"), tsz=10)

grid.arrange(
  grobs = list(p1,p2,p3,p4,p5),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 4, 5),
                        c(2, 4, 5),
                        c(3, 4, 5))
)


