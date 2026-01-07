#### Paper figures illustrating the latent deformation model and the simulation study ###########

library(gridExtra)
source("sim_fun.R")  #for simulation
source("aux_fun.R")  #for plots  
source("DB_est.R")   #depth-based estimation main function


ggcolors = c("aquamarine3", "deepskyblue", "blue3", "purple")


#### FIGURE 1 - MODEL ILLUSTRATION

m = 101
t = seq(0, 1, length.out = m)
p=4
ameans = rep(100,p)

simdata = sim_data(n=50, sigma_amp=4, sigma_warp=1, sigma_error=0, sigma_dist=0, ameans=ameans, p=4, t, m, k_psi=1, 
                   out=F)

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


### FIGURE 2 - WHyRA PLOT
p=2
n=50
# Top plot: sigma_dist=0      # No additional distortion other than individual warping functions
simdata = sim_data(n, sigma_amp=4, sigma_warp=1, sigma_error=0, sigma_dist=0, ameans=rep(100,p), p=2, t, m, k_psi=2, 
                   out=F)

  #Depth-based warping estimation
DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

  #MHI calculation (roahd package)
mhi=roahd::MHI(roahd::fData(t,DBest$h_hat_list[[1]]))

  #Warping estimates plot
id = sort(mhi, index.return=T,decreasing=T)$ix
col2=rep(1,n)
colores = rev(colorspace::heat_hcl(n,h=c(0,-100),l=c(75,40),c=c(40,80),power=1))
col2[id] = colores
names(col2) <- 1:n 

h_df=data.frame( id = rep(simdata$simdata[,1],p), t = rep(simdata$simdata[,2],p), group = rep(1:p,each=m*n),
                 warps = c(as.numeric(t(DBest$h_hat_list[[1]])),as.numeric(t(DBest$h_hat_list[[2]])) ) )

h_df$group <- factor(h_df$group,levels=c(1,2), labels=c(expression(X[1]),expression(X[2]) ))

pl<-ggplot(data=h_df, aes(x=t, y=warps,group=id,color=as.factor(id))) +
  geom_line() + 
  scale_color_manual(values=col2,guide="none")+labs(x="t" , y=expression(hat(h)[ij]), title="")+
  facet_wrap(~group,labeller = label_parsed) +
  theme_light() 


  #Warping MHI's scatter plot
Warps_mfData = roahd::mfData(t,DBest$h_hat_list)
mhi = data.frame(mhi1=roahd::MHI(Warps_mfData$fDList[[1]]), 
                 mhi2=roahd::MHI(Warps_mfData$fDList[[2]]))

pw<- ggplot(mhi, aes(x=mhi1,y=mhi2)) + 
  geom_point(size=3,color=col2) +
  xlim(0,1) + ylim(0,1)+
  geom_abline(intercept = 0, slope = 1) +
  labs(x=bquote(MHI~hat(h)[i1]) , y=bquote(MHI~hat(h)[i2]), title="WHyRA plot") +
  theme_light()

grid.arrange(
  grobs = list(pl,pw),
  widths = c(2, 1.2),
  layout_matrix = rbind(c(1, 2))
)

# Bottom plot: sigma_dist=1    # Additional individual time distortion 
simdata = sim_data(n, sigma_amp=4, sigma_warp=1, sigma_error=0, sigma_dist=1, ameans=rep(100,p), p=2, t, m, k_psi=2, 
                   out=F)

  #Depth-based warping estimation
DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

  #MHI calculation (roahd package)
mhi=roahd::MHI(roahd::fData(t,DBest$h_hat_list[[1]]))

  #Warping estimates plot
id = sort(mhi, index.return=T,decreasing=T)$ix
col2=rep(1,n)
colores = rev(colorspace::heat_hcl(n,h=c(0,-100),l=c(75,40),c=c(40,80),power=1))
col2[id] = colores
names(col2) <- 1:n 

h_df=data.frame( id = rep(simdata$simdata[,1],p), t = rep(simdata$simdata[,2],p), group = rep(1:p,each=m*n),
                 warps = c(as.numeric(t(DBest$h_hat_list[[1]])),as.numeric(t(DBest$h_hat_list[[2]])) ) )

h_df$group <- factor(h_df$group,levels=c(1,2), labels=c(expression(X[1]),expression(X[2]) ))

pl<-ggplot(data=h_df, aes(x=t, y=warps,group=id,color=as.factor(id))) +
  geom_line() + 
  scale_color_manual(values=col2,guide="none")+labs(x="t" , y=expression(hat(h)[ij]), title="")+
  facet_wrap(~group,labeller = label_parsed) +
  theme_light() 


  #Warping MHI's scatter plot
Warps_mfData = roahd::mfData(t,DBest$h_hat_list)
mhi = data.frame(mhi1=roahd::MHI(Warps_mfData$fDList[[1]]), 
                 mhi2=roahd::MHI(Warps_mfData$fDList[[2]]))

pw<- ggplot(mhi, aes(x=mhi1,y=mhi2)) + 
  geom_point(size=3,color=col2) +
  xlim(0,1) + ylim(0,1)+
  geom_abline(intercept = 0, slope = 1) +
  labs(x=bquote(MHI~hat(h)[i1]) , y=bquote(MHI~hat(h)[i2]), title="WHyRA plot") +
  theme_light()

grid.arrange(
  grobs = list(pl,pw),
  widths = c(2, 1.2),
  layout_matrix = rbind(c(1, 2))
)


### FIGURE 3 - SIMULATION SETTINGS
p=4
ameans = rep(100,p)
simdata = sim_data(n, sigma_amp=4, sigma_warp=0.5, sigma_error=0, sigma_dist=0, ameans=ameans, p=4, t, m, k_psi=1, 
                   out=F)

p1 <- plot_ggplot_target(t,L, main="Target function",yl=expression(lambda~(t)),lt="solid")

p2 <- psi_plot(t,matrix(unlist(simdata$psi),ncol=p),ggcolors, title="Component-based distorsions, setting 1", sim=T)

p3 <- hplot(t,simdata$hmat, title=bquote(atop(Warping~funct.~setting~1,sigma[W] == 0.5)), sim=T)

p=30
ameans = rep(100,p)
simdata = sim_data(n, sigma_amp=4, sigma_warp=1, sigma_error=0, sigma_dist=0, ameans=ameans, p=30, t, m, k_psi=2, 
                    out=F)

p4 <- psi_plot(t,matrix(unlist(simdata$psi),ncol=p),ggcolors, title="Component-based distorsions, setting 2", sim=T)
p5 <- hplot(t,simdata$hmat, title=bquote(atop(Warping~funct.~setting~1,sigma[W] == 1)), sim=T)

simdata = sim_data(n, sigma_amp=4, sigma_warp=5, sigma_error=0, sigma_dist=0, ameans=ameans, p=30, t, m, k_psi=2, 
                    out=F)

p6 <- hplot(t,simdata$hmat, title=bquote(atop(Warping~funct.~setting~2, epsilon[W] == 0.005) ), sim=T)

simdata = sim_data(n, sigma_amp=4, sigma_warp=7.5, sigma_error=0, sigma_dist=0, ameans=ameans, p=30, t, m, k_psi=2, 
                    out=F)

p7 <- hplot(t,simdata$hmat, title=bquote(atop(Warping~funct.~setting~2,epsilon[W] == 0.0075)), sim=T)

grid.arrange(
  grobs = list(p1,p2,p3,p4,p5,p6,p7),
  widths = c(1, 1, 1, 1,1,1, 1,1,1, 1,1,1),
  layout_matrix = rbind(c(1, 1,1,1,2,2,2,2,4,4,4,4),
                        c(3,3,3,5,5,5,6,6,6,7,7,7))
)


#### FIGURE 4 - SIMULATION RESULTS
source("est_comp.R")  #Implementation of CM estimation method. NOT PROVIDED: To be downloaded from 
### https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data

### 1st row
p=4
ameans = rep(100,p)

 #Data simulation
simdata = sim_data(n, sigma_amp=4, sigma_warp=0.5, sigma_error=0, sigma_dist=0.5, ameans=ameans, p=4, t, m, k_psi=1, 
                    out=F)

 #Simulated curves plot
p1 <- plot_ggplot_fun_df(simdata$simdata, f=NULL, a=ameans, m, p, t, 
                         main=bquote(Simulated~curves*": "*psi[j]~set.~1~p==4*", "*h[i]~set.~1~sigma[W]==0.5*", "*~sigma[D]==0.5*", "*~sigma[E] ==0),
                         names_comp=c("t", "X[1]","X[2]","X[3]","X[4]"),hor=T)

  #Depth-based warping estimation
DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

  #Estimation method by Carroll & Muller (2023). Code not provided. To be downloaded from
  #https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data
test = est_components(data = simdata$simdata, smooth = TRUE)

 #Target function estimates plot 
est <- data.frame(time=rep(t,2),value=c(DBest$L_hat,test$L.hat),est=c(rep("depth",m),rep("CM",m)))

p1_est <- plot_ggplot_target(t,L, lt="solid",col="darkgray") + 
  geom_line(data=est, aes(x=time,y=value,group=est,colour = est, linetype = est), size=0.3) + 
  scale_linetype_manual(values=c(1,1)) +
  scale_colour_manual(values=c("red","black")) +   theme(strip.text.x = element_text(margin = margin(0.08,0,0.05,0, "cm")), #red es CM (1era por orden alfabético), black es depth
                                                         plot.margin = unit(c(0.1,0.2,0,-0.5), "cm")) #+

### 2nd row

  #Data simulation
simdata = sim_data(n, sigma_amp=4, sigma_warp=1, sigma_error=1, sigma_dist=0, ameans=ameans, p=4, t, m, k_psi=1, 
                    out=F)

  #Simulated curves plot
p2 <- plot_ggplot_fun_df(simdata$simdata, f=NULL, a=ameans, m, p, t, 
                             main=bquote(Simulated~curves*": "*psi[j]~set.~1~p==4*", "*h[i]~set.~1~sigma[W]==1*", "*~sigma[D]==0*", "*~sigma[E] ==1),
                             names_comp=c("t", "X[1]","X[2]","X[3]","X[4]"), hor=T)

  #Depth-based warping estimation
DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

  # Estimation method by Carroll & Muller (2023). Code not provided. To be downloaded from
  # https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data
test = est_components(data = simdata$simdata, smooth = TRUE)

  #Target function estimates plot 
est <- data.frame(time=rep(t,2),value=c(DBest$L_hat,test$L.hat),est=c(rep("depth",m),rep("CM",m)))

p2_est <- plot_ggplot_target(t,L, lt="solid",col="darkgray") + 
  geom_line(data=est, aes(x=time,y=value,group=est,colour = est, linetype = est), size=0.3) + 
  scale_linetype_manual(values=c(1,1)) +
  scale_colour_manual(values=c("red","black")) +   theme(strip.text.x = element_text(margin = margin(0.08,0,0.05,0, "cm")), #red es CM (1era por orden alfabético), black es depth
                                                         plot.margin = unit(c(0.1,0.2,0,-0.5), "cm")) #+

### 3rd row

  #Data simulation
simdata = sim_data(n, sigma_amp=4, sigma_warp=5, sigma_error=0, sigma_dist=1, ameans=ameans, p=4, t, m, k_psi=1, 
                    out=F)

  #Simulated curves plot
p3 <- plot_ggplot_fun_df(simdata$simdata, f=NULL, a=ameans, m, p, t, 
                             main=bquote(Simulated~curves*": "*psi[j]~set.~1~p==4*", "*h[i]~set.~2~epsilon[W]==0.005*", "*~sigma[D]==1*", "*~sigma[E] ==0),
                             names_comp=c("t", "X[1]","X[2]","X[3]","X[4]"), hor=T)


  #Depth-based warping estimation
DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

  # Estimation method by Carroll & Muller (2023). Code not provided. To be downloaded from
  # https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data
test = est_components(data = simdata$simdata, smooth = TRUE)

  #Target function estimates plot 
est <- data.frame(time=rep(t,2),value=c(DBest$L_hat,test$L.hat),est=c(rep("depth",m),rep("CM",m)))

p3_est <- plot_ggplot_target(t,L, lt="solid",col="darkgray") + 
  geom_line(data=est, aes(x=time,y=value,group=est,colour = est, linetype = est), size=0.3) + 
  scale_linetype_manual(values=c(1,1)) +
  scale_colour_manual(values=c("red","black")) +   theme(strip.text.x = element_text(margin = margin(0.08,0,0.05,0, "cm")), #red es CM (1era por orden alfabético), black es depth
                                                         plot.margin = unit(c(0.1,0.2,0,-0.5), "cm")) #+

### 4th row

  #Data simulation
simdata = sim_data(n, sigma_amp=4, sigma_warp=7.5, sigma_error=1, sigma_dist=0, ameans=ameans, p=4, t, m, k_psi=1, 
                    out=F)

  #Simulated curves plot
p4 <- plot_ggplot_fun_df(simdata$simdata, f=NULL, a=ameans, m, p, t, 
                             main=bquote(Simulated~curves*": "*psi[j]~set.~1~p==4*", "*h[i]~set.~2~epsilon[W]==0.0075*", "*~sigma[D]==0*", "*~sigma[E] ==1),
                             names_comp=c("t", "X[1]","X[2]","X[3]","X[4]"),hor=T)

  #Depth-based warping estimation
DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

  # Estimation method by Carroll & Muller (2023). Code not provided. To be downloaded from
  # https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data
test = est_components(data = simdata$simdata, smooth = TRUE)

  #Target function estimates plot 
est <- data.frame(time=rep(t,2),value=c(DBest$L_hat,test$L.hat),est=c(rep("depth",m),rep("CM",m)))

p4_est <- plot_ggplot_target(t,L, lt="solid",col="darkgray") + 
  geom_line(data=est, aes(x=time,y=value,group=est,colour = est, linetype = est), size=0.3) + 
  scale_linetype_manual(values=c(1,1)) +
  scale_colour_manual(values=c("red","black")) +   theme(strip.text.x = element_text(margin = margin(0.08,0,0.05,0, "cm")), #red es CM (1era por orden alfabético), black es depth
                                                         plot.margin = unit(c(0.1,0.2,0,-0.5), "cm")) #+
### 5th row

  #Data simulation
p=30
simdata = sim_data(n, sigma_amp=4, sigma_warp=1, sigma_error=0, sigma_dist=0, ameans=rep(100,30), p=30, t, m, k_psi=2, 
                    out=F)

  #Simulated curves plot
p5 <- plot_ggplot_fun_df(simdata$simdata, f=NULL, a=rep(100,30), m, p=30, t, 
                             main=bquote(Simulated~curves*": "*psi[j]~set.~2~p==30*", "*h[i]~set.~1~sigma[W]==1*", "*~sigma[D]==0*", "*~sigma[E] ==0),
                             names_comp=c("t", "X[1]","X[2]","X[3]","X[4]"), hor=T)

  #Depth-based warping estimation
DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

  #Estimation method by Carroll & Muller (2023). Code not provided. To be downloaded from
  #https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data
test = est_components(data = simdata$simdata, smooth = TRUE)

  #Target function estimates plot 
est <- data.frame(time=rep(t,2),value=c(DBest$L_hat,test$L.hat),est=c(rep("depth",m),rep("CM",m)))

p5_est <- plot_ggplot_target(t,L, lt="solid",col="darkgray") + 
  geom_line(data=est, aes(x=time,y=value,group=est,colour = est, linetype = est), size=0.3) + 
  scale_linetype_manual(values=c(1,1)) +
  scale_colour_manual(values=c("red","black")) +   theme(strip.text.x = element_text(margin = margin(0.08,0,0.05,0, "cm")), #red es CM (1era por orden alfabético), black es depth
                                                         plot.margin = unit(c(0.1,0.2,0,-0.5), "cm")) #+

grid.arrange(
  grobs = list(p1,p1_est,p2,p2_est,p3, p3_est, p4, p4_est, p5, p5_est),
  widths = c(4,1),
  layout_matrix = rbind(c(1,2),c(3,4),c(5,6),c(7,8),c(9,10)))


#### FIGURE 11 - SIMULATION RESULTS WITH OUTLIERS

### 1st row 

  #Data simulation
p=4
ameans = rep(100, p)
outmod=T 

simdata = sim_data(n, sigma_amp=4, sigma_warp=0.5, sigma_error=0, sigma_dist=0.5, ameans=ameans, p=4, t, m, k_psi=1, 
                    out=outmod, out_prop=0.05)

  #Simulated curves plot
p1 <- plot_ggplot_fun_df(simdata$simdata, f=NULL, a=ameans, m, p, t, out=T,outs = simdata$outs, 
                        main=bquote(Simulated~curves*": "*psi[j]~set.~1~p==4*", "*h[i]~set.~1~sigma[W]==0.5*", "*~sigma[D]==0.5*", "*~sigma[E] ==0*", "*~c ==0.05),
                        names_comp=c("t", "X[1]","X[2]","X[3]","X[4]"), hor=T)

  #Depth-based warping estimation
DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

  #Estimation method by Carroll & Muller (2023). Code not provided. To be downloaded from
  #https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data
test = est_components(data = simdata$simdata, smooth = TRUE)

  #Target function estimates plot 
est <- data.frame(time=rep(t,2),value=c(DBest$L_hat,test$L.hat),est=c(rep("depth",m),rep("CM",m)))

p1_est <-  plot_ggplot_target(t,L, lt="solid",col="darkgray") + 
  geom_line(data=est, aes(x=time,y=value,group=est,colour = est, linetype = est), size=0.3) + 
  scale_linetype_manual(values=c(1,1)) +
  scale_colour_manual(values=c("red","black")) +   theme(strip.text.x = element_text(margin = margin(0.08,0,0.05,0, "cm")), #red es CM (1era por orden alfabético), black es depth
                                                         plot.margin = unit(c(0.1,0.2,0,-0.5), "cm")) #+

### 2nd row 

  #Data simulation

simdata = sim_data(n, sigma_amp=4, sigma_warp=0.5, sigma_error=0, sigma_dist=0.5, ameans=ameans, p=4, t, m, k_psi=1, 
                    out=outmod, out_prop=0.1)

  #Simulated curves plot
p2 <- plot_ggplot_fun_df(simdata$simdata, f=NULL, a=ameans, m, p, t, out=T,outs = simdata$outs, 
                             main=bquote(Simulated~curves*": "*psi[j]~set.~1~p==4*", "*h[i]~set.~1~sigma[W]==0.5*", "*~sigma[D]==0.5*", "*~sigma[E] ==0*", "*~c ==0.1),
                             names_comp=c("t", "X[1]","X[2]","X[3]","X[4]"), hor=T)

  #Depth-based warping estimation
DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

  #Estimation method by Carroll & Muller (2023). Code not provided. To be downloaded from
  #https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data
test = est_components(data = simdata$simdata, smooth = TRUE)

  #Target function estimates plot 
est <- data.frame(time=rep(t,2),value=c(DBest$L_hat,test$L.hat),est=c(rep("depth",m),rep("CM",m)))

p2_est <-  plot_ggplot_target(t,L, lt="solid",col="darkgray") + 
  geom_line(data=est, aes(x=time,y=value,group=est,colour = est, linetype = est), size=0.3) + 
  scale_linetype_manual(values=c(1,1)) +
  scale_colour_manual(values=c("red","black")) +   theme(strip.text.x = element_text(margin = margin(0.08,0,0.05,0, "cm")), #red es CM (1era por orden alfabético), black es depth
                                                         plot.margin = unit(c(0.1,0.2,0,-0.5), "cm")) #+

### 3rd row 

  #Data simulation
simdata = sim_data(n, sigma_amp=4, sigma_warp=5, sigma_error=0, sigma_dist=1, ameans=ameans, p=4, t, m, k_psi=1, 
                    out=outmod, out_prop=0.05)

  #Simulated curves plot
p3 <- plot_ggplot_fun_df(simdata$simdata, f=NULL, a=ameans, m, p, t, out=T,outs = simdata$outs,
                             main=bquote(Simulated~curves*": "*psi[j]~set.~1~p==4*", "*h[i]~set.~2~epsilon[W]==0.005*", "*~sigma[D]==1*", "*~sigma[E] ==0*", "*~c ==0.05),
                             names_comp=c("t", "X[1]","X[2]","X[3]","X[4]"), hor=T)

  #Depth-based warping estimation
DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

  #Estimation method by Carroll & Muller (2023). Code not provided. To be downloaded from
  #https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data
test = est_components(data = simdata$simdata, smooth = TRUE)

  #Target function estimates plot 
est <- data.frame(time=rep(t,2),value=c(DBest$L_hat,test$L.hat),est=c(rep("depth",m),rep("CM",m)))

p3_est <-plot_ggplot_target(t,L, lt="solid",col="darkgray") + 
  geom_line(data=est, aes(x=time,y=value,group=est,colour = est, linetype = est), size=0.3) + 
  scale_linetype_manual(values=c(1,1)) +
  scale_colour_manual(values=c("red","black")) +   theme(strip.text.x = element_text(margin = margin(0.08,0,0.05,0, "cm")), #red es CM (1era por orden alfabético), black es depth
                                                         plot.margin = unit(c(0.1,0.2,0,-0.5), "cm")) #+

### 4th row 

  #Data simulation

simdata = sim_data(n, sigma_amp=4, sigma_warp=5, sigma_error=0, sigma_dist=1, ameans=ameans, p=4, t, m, k_psi=1, 
                    out=outmod, out_prop=0.1)

  #Simulated curves plot
p4 <- plot_ggplot_fun_df(simdata$simdata, f=NULL, a=ameans, m, p, t, out=T,outs = simdata$outs,
                             main=bquote(Simulated~curves*": "*psi[j]~set.~1~p==4*", "*h[i]~set.~2~epsilon[W]==0.005*", "*~sigma[D]==1*", "*~sigma[E] ==0*", "*~c ==0.1),
                             names_comp=c("t", "X[1]","X[2]","X[3]","X[4]"), hor=T)
  #Depth-based warping estimation
DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

  #Estimation method by Carroll & Muller (2023). Code not provided. To be downloaded from
  #https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data
test = est_components(data = simdata$simdata, smooth = TRUE)

  #Target function estimates plot 
est <- data.frame(time=rep(t,2),value=c(DBest$L_hat,test$L.hat),est=c(rep("depth",m),rep("CM",m)))

p4_est <- plot_ggplot_target(t,L, lt="solid",col="darkgray") +
  geom_line(data=est, aes(x=time,y=value,group=est,colour = est, linetype = est), size=0.3) + 
  scale_linetype_manual(values=c(1,1)) +
  scale_colour_manual(values=c("red","black")) +   theme(strip.text.x = element_text(margin = margin(0.08,0,0.05,0, "cm")), #red es CM (1era por orden alfabético), black es depth
                                                         plot.margin = unit(c(0.1,0.2,0,-0.5), "cm")) #+

grid.arrange(
  grobs = list(p1,p1_est,p2,p2_est,p3, p3_est, p4, p4_est),
  widths = c(4,1),
  layout_matrix = rbind(c(1,2),c(3,4),c(5,6),c(7,8)))
