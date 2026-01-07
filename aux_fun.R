#### Auxiliar functions for plots 
#### Modified versions of functions from Carroll&Muller,Biometrics, 2023. Original code is available at
### https://academic.oup.com/biometrics/article/79/4/3345/7585581?login=false#supplementary-data

library(reshape2)
library(ggplot2)
ggcolors = c("aquamarine3", "deepskyblue", "blue3", "purple")


plot_ggplot_target<- function(t,L, main="Target function",yl="",xl="t",col="black",sz=1.5, lt="dashed",tsz=12){
  
  df_plot = data.frame(t,L(t))
  names(df_plot) <- c("t","L")
  
  
  ggp = ggplot(df_plot) + geom_line(aes(x=t, y = L), color=col, size=sz, linetype = lt) +
    xlab(xl) +
    ylab(yl)+
    #ylim(c(0.5,1.5)) +
    ggtitle(main) +
    theme(legend.position = "none",
          #axis.text.y = element_blank(),
          axis.ticks = element_blank(),plot.title=element_text(size=tsz))
  ggp
  return(ggp)
  
}

plot_ggplot_funct <- function(f, a=rep(1,p), m, p, t, main=NULL,xl="t",names_comp=c("t", "A","B","C","D"),true_fun=NULL,
                              tsz=12){ ## f is a function  to evaluate on 1:p
  df_plot = matrix(0, nrow = m, ncol = p+1)
  
  df_plot[,1] = t
  
  for(j in 1:p){
    df_plot[,j+1] = a[j]*f(j)
  }
  
  df_plot = as.data.frame(df_plot)
  # names(df_plot) = c("t", "A","B","C","D")
  names(df_plot) = names_comp
  ggdata =reshape2::melt(df_plot , id.vars = "t")
  
  ggp = ggplot(ggdata, aes(x=t, y = value, color = variable))
  
  
  if(!is.null(true_fun)){
    ggdata$true = unlist(lapply(1:p,true_fun))  #true_fun is a function with argument j=1,...,p
    
    ggp  = ggp  + geom_line(data = ggdata, aes(x = t, y = true,
                                                 group = variable), linetype = "dashed",color="black") 
  }  
  
  ggp = ggp + xlab(xl) + ylab("")+
    #ylim(c(0.5,1.5)) +
    geom_line(size=1.15, alpha = .8) +
    facet_wrap(~variable) +
    scale_color_manual(values=ggcolors) +
    facet_wrap(~variable, ncol = 1, labeller=label_parsed) + #label_parsed to get subscript in variable names
    ggtitle(main) +
    theme(legend.position = "none",
          #axis.text.y = element_blank(),
          axis.ticks = element_blank(), plot.title=element_text(size=tsz))
  ggp
  return(ggp)
}

plot_ggplot_fun_df <- function(df, f=NULL, a=rep(1,p), m, p, t, main=NULL,xl="t",names_comp=c("t", "A","B","C","D"),
                               out=F,outs=c(),tsz=12,hor=F){ ## df1 is a data.frame with individuals, f is a function with with target curves
  names(df) <- c("id",names_comp)
  
  if (p<=4){
    x=1:4
  }else{
    x=sort(sample(1:p,4,replace=F))
    p=4
    names_comp=c("t", paste0("X[",x,"]") )
    df = df[,c(1:2,x+2)]
    names(df)[3:6] <- paste0("X[",x,"]")
  }
  
  ggdata=reshape2::melt(df, id.vars = c("id","t"))   ##col 1 y 2 son id y t
  
  if (!is.null(f)){
  df_plot = matrix(0, nrow = m, ncol = p+1)
  
  df_plot[,1] = t
  
  for(j in 1:p){
    df_plot[,j+1] = a[j]*f(j)
  }
  
  df_plot = as.data.frame(df_plot)
  # names(df_plot) = c("t", "A","B","C","D")
  names(df_plot) = names_comp
  ggdata2 = reshape2::melt(df_plot , id.vars = "t")
  }
  

  ggsim = ggplot()+
    geom_line(data = ggdata,
              aes(x=t, y = value, group=id, color = variable),
              size=.6, alpha = 0.2) 
  if (out==T){
    ggdata_sub <- ggdata %>% dplyr::select(id,variable) %>% mutate(indices=paste(id,as.numeric(variable)))  #as.numeric preserves alphabetical order in factor
    outs_p = paste(outs[,1],outs[,2])
    rows=which(!is.na(match(ggdata_sub$indices,outs_p)))
    ggdata_outs <- ggdata[rows,]
    ggsim = ggsim + geom_line(data = ggdata_outs,
                              aes(x=t, y = value, group=id),
                              color = "darkgray",linetype = "dashed") 
  }
  if (!is.null(f)){
  ggsim = ggsim +
    geom_line(data = ggdata2,
              aes(x=t, y = value),
              color = "black") 
  }
  ggsim = ggsim +
    xlab(xl) +
    ylab("")
  
  if(hor==T){ #horizontal
    ggsim = ggsim + facet_grid(~variable,  labeller=label_parsed)  #label_parsed to get subscript in variable names
  }else{ #vertical
    ggsim = ggsim + facet_wrap(~variable, ncol = 1, labeller=label_parsed) #label_parsed to get subscript in variable names
  }
   ggsim = ggsim + scale_color_manual(values=ggcolors) +
    ggtitle(main) +
    theme(legend.position = "none",
          #axis.text.y = element_blank(),
          axis.ticks = element_blank(),plot.title=element_text(size=tsz)) +
    scale_x_continuous(expand = c(0.01,0.01))
  
  ggsim
  return(ggsim)
}



###  plot psi, h  ##############################################################
psi_plot <- function(t,psi_mat,ggcolors, title="Component-based distorsions, setting 1",
                     yl=bquote(psi[j]~(t)),sim=F,true_psi=NULL){
  #psi_mat is an m*p matrix
  
  psi_df = melt(data.frame(psi_mat))
  p = ncol(psi_mat)
  psi_df$t=rep(t,p)
  

  
  ggpsi = ggplot(data = psi_df,
                 aes(x = t, y = value,
                     group = variable, color = variable)) +
    geom_line(size = ifelse(p==4,1.15,0.6), alpha=ifelse(p==4,1,0.5)) 
  
  if(!is.null(true_psi)){
    psi_df$true = unlist(true_psi)  #true_psi is a p-length list with the generated psi's
  
    ggpsi = ggpsi + geom_line(data = psi_df, aes(x = t, y = true,
                                                 group = variable), linetype = "dashed", color = "black") 
  }  
  
  ggpsi = ggpsi + ylab(yl)+
    scale_color_manual(values=colorRampPalette(ggcolors)(p)) +
    theme(legend.position = "none",plot.title=element_text(size=10),plot.margin = margin(0,0.3,0,0, "cm")) +
    ggtitle(title)
  
  if (sim==T){
    ggpsi = ggpsi  +
    geom_vline(xintercept = 0.9,col="violetred", linetype="dashed") +
    annotate("text", x=0.8, y=0.075, label= "t=0.9",col="violetred") 
    }
  
    ggpsi = ggpsi +
    geom_abline(col="black",size=0.6, linetype="longdash") +
    coord_cartesian(xlim = c(min(t),max(t)), ylim = c(min(t),max(t)),
                    expand = FALSE)
  
  return(ggpsi)
}

hplot <- function(t,hmat, title=bquote(Warping~funct.~setting~1*", "*sigma[W] == 0.5), yl=bquote(h[i]~(t)),sim=F){
  #hmat is an n*m matrix
  
  h_df = melt(data.frame(t(hmat)))
  n = nrow(hmat)
  h_df$t=rep(t,n)
  
  ggh = ggplot(data = h_df,
               aes(x = t, y = value,
                   group = variable) )+
    geom_line(size = 0.5, col="lightgray") +
    ylab(yl) +
    theme(legend.position = "none",plot.title=element_text(size=10),plot.margin = margin(0,0.3,0,0, "cm")) +
    ggtitle(title) 
  
  if (sim==T){
    ggh = ggh  +
      geom_vline(xintercept = 0.9,col="violetred", linetype="dashed") +
      annotate("text", x=0.8, y=0.075, label= "t=0.9",col="violetred") 
  }
  
  ggh = ggh  +
    geom_abline(col="black",size=0.6, linetype="longdash") +
    coord_cartesian(xlim = c(min(t),max(t)), ylim = c(min(t),max(t)),
                    expand = FALSE)
  
  return(ggh)
}

