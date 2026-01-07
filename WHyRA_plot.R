#### WHyRA plot #####

source("sim_fun.R")  #for simulation
source("aux_fun.R")  #for plots  
source("DB_est.R")   #depth-based estimation main function

## 1. Data simulation
m = 101
t = seq(0, 1, length.out = m)

k_psi=2
p=2
n=50
sigma_error=0
sigma_warp=1

# sigma_dist=0      # No additional distortion other than individual warping functions
sigma_dist=1    # Additional individual time distortion 
simdata = sim_data(n, sigma_amp=4, sigma_warp=1, sigma_error=0, sigma_dist=1, ameans=rep(100,p), p=2, t, m, k_psi=2, 
                   out=F)

## 2. Depth-based warping estimation
DBest <- DB_est(simdata$simdata, smooth=TRUE, method = "MA")

## 3. MHI calculation (roahd package)
mhi=roahd::MHI(roahd::fData(t,DBest$h_hat_list[[1]]))

## 4. Warping estimates plot
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


## 5. Warping MHI's scatter plot
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