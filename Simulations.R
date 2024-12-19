source("~/Universitet-mat/.Speciale/R/splines_combined.R")
logistic_bivariate<-function(x,y,rho){
  exp(-x-y)/((1+exp(-x))^2*(1+exp(-y))^2*exp((rho*((x-y))^2)/2))
}

#We consider intervals [-5,5]
gen_sym<-function(n,rho=0){
  i<-1
  M<-2/16 #This just needs to be the maximum, since we do not multiply by the constant g
  sample<-data.frame("x"=NA,"y"=NA)
  while(i<=n){
    X<-runif(1,min=-5,max=5); Y<-runif(1,min=-5,max=5); U<-runif(1)
    if(U<=logistic_bivariate(X,Y,rho)/(M)){
      sample[i,]<-c(X,Y)
      i<-i+1
    }
  }
  return(sample)
}


gen_asym<-function(n,rho=0){
  eps1<-rnorm(n)
  eps2<-rnorm(n)
  Theta<-runif(n,0,2*pi)
  sample<-data.frame("x"=rho*cos(Theta)+eps1/4,"y"=rho*sin(Theta)+eps2/4)
  return(sample)
}

library(VGAM)
asym_density<-function(x,y,rho){
  R<-sqrt((x^2+y^2))
  outcome<-drice(R,1/4,rho)/R
  return(outcome)
}



sample_sym<-gen_sym(1000,rho=0)
x<-sample_sym$x
y<-sample_sym$y
knots_x_inner<-c(min(x),-2,seq(-1,1,length.out=4),2,max(x))
knots_y_inner<-c(min(y),-2,seq(-1,1,length.out=4),2,max(y))

cross_validation<-cross_validate2D(x,y,knots_x_inner=knots_x_inner,knots_y_inner=knots_y_inner,k=2,l=2,u=1,v=1)
library(profvis)
profvis(zbSpline2D(x,y,knots_x_inner,knots_y_inner))

test<-zbSpline2D(x,y,knots_x_inner,knots_y_inner)
test$hist_data
plot(test,scale="density",plot_hist=TRUE)
#plot(cross_validation)


biv_sym<-zbSpline2D(x,y,knots_x_inner=knots_x_inner,alfa=0.9,bin_selection=c(10,10),knots_y_inner=knots_y_inner,k=3,l=3,u=1,v=1,res=100)
plot(biv_sym,scale="density")

trapz(biv_sym$seq_x,biv_sym$Z_spline_x)

biv_sym$rsd

plot(biv_sym,scale="density",what="full",type="interactive",plot_hist=TRUE,title="rho=0")

z_values<-kde2d(x,y,n=100)
grid<-expand.grid(biv_sym$midpoints_x,biv_sym$midpoints_y)
norm_constant_sym<-trapz2d(biv_sym$midpoints_x,biv_sym$midpoints_y,biv_sym$hist_data)
plot_ly(x=z_values$x,y=z_values$y,z=z_values$z,type="surface")%>%
  add_markers(x=grid$Var1,y=grid$Var2,z=as.vector(biv_sym$hist_data)/norm_constant,marker = list(size = 3, color = "black"))

#x_vals<-seq(-5,5,length.out=100)->y_vals
par(mfrow=c(2,3),mar=c(0,0,2,0))


x_vals<-seq(-5,5,length.out=100)->y_vals
z_vals<-outer(x_vals,y_vals,FUN=function(x,y){logistic_bivariate(x,y,rho=0)})
persp3D(x_vals,y_vals,z_vals/trapz2d(x_vals,y_vals,z_vals),theta=30,phi=90,col=viridis(50),ticktype="detailed",axes=FALSE,main=expression(rho==0),colkey=FALSE,cex.main=2)
z_vals<-outer(x_vals,y_vals,FUN=function(x,y){logistic_bivariate(x,y,rho=0.5)})
persp3D(x_vals,y_vals,z_vals/trapz2d(x_vals,y_vals,z_vals),theta=30,phi=90,col=viridis(50),ticktype="detailed",axes=FALSE,main=expression(rho==0.5),colkey=FALSE,cex.main=2)
z_vals<-outer(x_vals,y_vals,FUN=function(x,y){logistic_bivariate(x,y,rho=1)})
persp3D(x_vals,y_vals,z_vals/trapz2d(x_vals,y_vals,z_vals),theta=30,phi=90,col=viridis(50),ticktype="detailed",axes=FALSE,main=expression(rho==1),colkey=FALSE,cex.main=2)

x_vals<-seq(-2,2,length.out=100)->y_vals
z_vals<-outer(x_vals,y_vals,FUN=function(x,y){asym_density(x,y,rho=0)})
persp3D(x_vals,y_vals,z_vals/trapz2d(x_vals,y_vals,z_vals),theta=30,phi=60,col=viridis(50),ticktype="detailed",axes=FALSE,main=expression(rho==0),colkey=FALSE,cex.main=2)
z_vals<-outer(x_vals,y_vals,FUN=function(x,y){asym_density(x,y,rho=0.5)})
persp3D(x_vals,y_vals,z_vals/trapz2d(x_vals,y_vals,z_vals),theta=30,phi=60,col=viridis(50),ticktype="detailed",axes=FALSE,main=expression(rho==0.5),colkey=FALSE,cex.main=2)
z_vals<-outer(x_vals,y_vals,FUN=function(x,y){asym_density(x,y,rho=1)})
persp3D(x_vals,y_vals,z_vals/trapz2d(x_vals,y_vals,z_vals),theta=30,phi=60,col=viridis(50),ticktype="detailed",axes=FALSE,main=expression(rho==1),colkey=FALSE,cex.main=2)



x_vals<-seq(-(1+rho),1+rho,length.out=100)->y_vals
z_vals<-outer(x_vals,y_vals,FUN=function(x,y){asym_density(x,y,rho)})
persp3D(x=x_vals,y=y_vals,z=z_vals/trapz2d(x_vals,y_vals,z_vals),phi=70,col=viridis(50),ticktype="detailed",colkey=F,main=expression(rho==0))


sim<-gen_sym(2000,rho=2)
x<-sim$x
y<-sim$y
z_values<-kde2d(x,y,n=100)
plot_ly(x=z_values$x,y=z_values$y,z=z_values$z,type="surface")
grid<-expand.grid(biv_sym$midpoints_x,biv_sym$midpoints_y)
norm_constant<-trapz2d(biv_sym$midpoints_x,biv_sym$midpoints_y,biv_sym$hist_data)
plot_ly(x=z_values$x,y=z_values$y,z=z_values$z,type="surface")#%>%
  #add_markers(x=grid$Var1,y=grid$Var2,z=as.vector(biv_sym$hist_data)/norm_constant,marker = list(size = 3, color = "black"))


kx<-seq(min(x),max(x),length.out=5)
ky<-seq(min(y),max(y),length.out=5)

library(profvis)
profvis(bivariate(x,y,knots_x_inner=kx,knots_y_inner=ky,bin_selection=scott,alfa=1,k=3,l=3,u=1,v=1,res=50))

biv_asym<-zbSpline2D(x,y,knots_x_inner=kx,knots_y_inner=ky,bin_selection=scott,alfa=1,k=3,l=3,u=1,v=1,res=50)
biv_asym$rsd

cv<-cross_validate2D(x,y,knots_x_inner=kx,knots_y_inner=ky,k=3,l=3,u=3,v=1)
plot(cv)

par(mfrow=c(2,3),mar=c(0,0,1,0))
jpeg("clr_full.jpeg")
plot(biv_asym,scale="clr",what="full",type="static",plot_hist=TRUE,title="clr full")
dev.off()
jpeg("density_full.jpeg")
plot(biv_asym,scale="clr",what="interaction",type="static",plot_hist=FALSE,title="clr interaction")
dev.off()
jpeg("clr_interaction.jpeg")
plot(biv_asym,scale="clr",what="independent",type="static",plot_hist=FALSE,title="clr independent")
dev.off()
jpeg("density_interaction.jpeg")
plot(biv_asym,scale="density",what="full",type="static",plot_hist=TRUE,title="density full")
dev.off()
jpeg("clr_independent.jpeg")
plot(biv_asym,scale="density",what="interaction",type="static",plot_hist=FALSE,title="density interaction")
dev.off()
jpeg("density_independent.jpeg")
plot(biv_asym,scale="density",what="independent",type="static",plot_hist=FALSE,title="density independent")
dev.off()

plot(biv_asym,what="geom_X",scale="density",plot_hist=TRUE)

trapz2d(biv_asym$seq_x,biv_asym$seq_y,biv_asym$C_spline)

x_vals<-seq(-2,2,length.out=50)->y_vals
clr_asym<-clr2d(function(x,y){asym_density(x,y,1)},x_vals,y_vals)
z_vals<-outer(x_vals,y_vals,FUN=clr_asym)
persp3D(x_vals,y_vals,z_vals,theta=30,phi=60,col=viridis(50),ticktype="detailed",main=expression(rho==1),colkey=FALSE) 
#Maybe add the plots of clr-transformed things


#We perform the power analysis on the centered data first
#We take 8 values of rho from 0 to 1, and run 100 tests in each case
#These tests are run over simulations of size 2000

library(copula)
library(dHSIC)

cor.test(x,y,method="spearman")$p.value
dhsic.test(x,y,method="gamma")$p.value

#Det kunne være sejt med et konfidensinterval

ind_sim<-indepTestSim(2000,2,N=1000)
indepTest(sample_asym,ind_sim)

rho_seq<-seq(0,1,length.out=10)
p_vals_asym<-data.frame("rho"=sort(rep(rho_seq,100)),"rsd_scott"=NA,"rsd_doane"=NA,"spearman"=NA,"dhsic"=NA,"copula"=NA,"chisq_scott"=NA,"chisq_doane"=NA)

#Der sker godt nok nogle mærkelige ting for at få det til at køre lidt hurtigere.
#Vi imputerer for eksempel med en lidt hurtigere algoritme

#pb <- txtProgressBar(min = 0, max = nrow(p_vals_asym), style = 3)
for(i in 1:nrow(p_vals_asym)){
  #setTxtProgressBar(pb, i)
  print(i)
  simulation<-gen_asym(n=2000,rho=p_vals_asym$rho[i])
  x<-simulation$x
  y<-simulation$y
  kx<-seq(min(x),max(x),length.out=5)
  ky<-seq(min(x),max(x),length.out=5)
  #biv<-bivariate(x,x,alfa=1,bin_selection=scott,knots_x_inner=kx,knots_y_inner=kx,k=2,l=2,u=1,v=1,res=100)
  biv_scott<-perm_test(x,y,kx,ky,bin_selection=scott,k=2,l=2,u=1,v=1,K=1000)
  hist_scott<-biv_scott$hist
  biv_doane<-perm_test(x,y,kx,ky,bin_selection=doane,k=2,l=2,u=1,v=1,K=1000)
  hist_doane<-biv_doane$hist
  
  p_vals_asym$rsd_scott[i]<-biv_scott$p_val
  p_vals_asym$rsd_doane[i]<-biv_doane$p_val
  p_vals_asym$spearman[i]<-cor.test(x,y,method="spearman")$p.value #Simply use spearman test
  p_vals_asym$dhsic[i]<-dhsic.test(x,y,method="gamma")$p.value #Gaussian kernel with median value as bandwidth
  p_vals_asym$copula[i]<-indepTest(simulation,ind_sim)$pvalues #We use global statistic
  p_vals_asym$chisq_scott[i]<-chisq.test(hist_scott)$p.value
  p_vals_asym$chisq_doane[i]<-chisq.test(hist_doane)$p.value
}

saveRDS(p_vals_asym,"asymetric p-values")
p_vals_asym<-readRDS("asymetric p-values")

getwd()

p_vals_summary <- p_vals_asym %>%
  pivot_longer(cols = -rho, names_to = "test", values_to = "p_value") %>%
  group_by(rho, test) %>%
  summarise(
    mean_p = median(p_value, na.rm = FALSE),
    lower_q = quantile(p_value, 0.25, na.rm = FALSE),
    upper_q = quantile(p_value, 0.75, na.rm = FALSE)
  )

p_vals_power_asym<-p_vals_asym %>%
  pivot_longer(cols = -rho, names_to = "test", values_to = "p_value") %>%
  group_by(rho, test) %>%
  summarise(
    power=length(which(p_value<0.05))
  )

# Plotting the means and interquartile range
ggplot(p_vals_summary, aes(x = rho, color = test)) +
  geom_line(aes(y = mean_p), linewidth = 1) +
  geom_ribbon(aes(ymin = lower_q, ymax = upper_q, fill = test), alpha = 0.2) +
  labs(
    title = "Development of p-values as a function of rho",
    subtitle="For perimeter data",
    x = expression(rho),
    y = "p-value/rsd"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")

p2<-ggplot(p_vals_power_asym) +
  geom_line(mapping = aes(
    x = rho,
    y = power,
    color = case_when(
      test %in% c("rsd_scott", "rsd_doane") ~ "rsd",
      test %in% c("chisq_scott", "chisq_doane") ~ "chisq",
      TRUE ~ test
    ),
    linetype = case_when(
      test %in% c("rsd_scott", "chisq_scott") ~ "Scott",
      test %in% c("rsd_doane", "chisq_doane") ~ "Doane",
      TRUE ~ "Scott"
    )
  ), linewidth = 1) +
  ylab("") +
  ggtitle("Perimeter")+
  xlab(expression(rho)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1",name="Test") +
  scale_linetype_manual(values = c("Scott"="solid","Doane"="dashed"),name="Bin selection")

p2

rho_seq<-seq(0,0.2,length.out=10)
p_vals_sym<-data.frame("rho"=sort(rep(rho_seq,100)),"rsd_scott"=NA,"rsd_doane"=NA,"spearman"=NA,"dhsic"=NA,"copula"=NA,"chisq_scott"=NA,"chisq_doane"=NA)
#pb <- txtProgressBar(min = 0, max = nrow(p_vals_sym), style = 3)
print("CENTERED")
for(i in 1:nrow(p_vals_sym)){
#  setTxtProgressBar(pb, i)
  print(i)
  simulation<-gen_sym(n=2000,rho=p_vals_sym$rho[i])
  x<-simulation$x
  y<-simulation$y
  kx<-seq(min(x),max(x),length.out=5)
  ky<-seq(min(x),max(x),length.out=5)
  #biv<-bivariate(x,x,alfa=1,bin_selection=scott,knots_x_inner=kx,knots_y_inner=kx,k=2,l=2,u=1,v=1,res=100)
  biv_scott<-perm_test(x,y,kx,ky,bin_selection=scott,k=2,l=2,u=1,v=1,K=1000)
  hist_scott<-biv_scott$hist
  biv_doane<-perm_test(x,y,kx,ky,bin_selection=doane,k=2,l=2,u=1,v=1,K=1000)
  hist_doane<-biv_doane$hist
  
  p_vals_sym$rsd_scott[i]<-biv_scott$p_val
  p_vals_sym$rsd_doane[i]<-biv_doane$p_val
  p_vals_sym$spearman[i]<-cor.test(x,y,method="spearman")$p.value #Simply use spearman test
  p_vals_sym$dhsic[i]<-dhsic.test(x,y,method="gamma")$p.value #Gaussian kernel with median value as bandwidth
  p_vals_sym$copula[i]<-indepTest(simulation,ind_sim)$pvalues #We use global statistic
  p_vals_sym$chisq_scott[i]<-chisq.test(hist_scott)$p.value
  p_vals_sym$chisq_doane[i]<-chisq.test(hist_doane)$p.value
}

saveRDS(p_vals_sym,"symmetric pvalues")
p_vals_sym<-readRDS("C:/Users/espri/OneDrive/Dokumenter/symmetric pvalues")


p_vals_summary <- (p_vals_sym[1:400,]) %>%
  pivot_longer(cols = -rho, names_to = "test", values_to = "p_value") %>%
  group_by(rho, test) %>%
  summarise(
    mean_p = median(p_value, na.rm = FALSE),
    lower_q = quantile(p_value, 0.25, na.rm = FALSE),
    upper_q = quantile(p_value, 0.75, na.rm = FALSE)
  )

p_vals_power_sym<-p_vals_sym %>%
  pivot_longer(cols = -rho, names_to = "test", values_to = "p_value") %>%
  group_by(rho, test) %>%
  summarise(
    power=length(which(p_value<0.05))
  )

# Plotting the means and interquartile range
ggplot(p_vals_summary, aes(x = rho, color = test)) +
  geom_line(aes(y = mean_p), linewidth = 1) +
  geom_ribbon(aes(ymin = lower_q, ymax = upper_q, fill = test), alpha = 0.2) +
  labs(
    title = expression("Development of p-values as a function of "*rho),
    subtitle="For centered data",
    x = expression(rho),
    y = "p-value"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")

p1<-ggplot(p_vals_power_sym) +
  geom_line(mapping = aes(
    x = rho,
    y = power,
    color = case_when(
      test %in% c("rsd_scott", "rsd_doane") ~ "rsd",
      test %in% c("chisq_scott", "chisq_doane") ~ "chisq",
      TRUE ~ test
    ),
    linetype = case_when(
      test %in% c("rsd_scott", "chisq_scott") ~ "Scott",
      test %in% c("rsd_doane", "chisq_doane") ~ "Doane",
      TRUE ~ "Scott"
    )
  ), linewidth = 1) +
  ylab("Power (%)") +
  ggtitle("Centered")+
  xlab(expression(rho)) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1",name="Test") +
  scale_linetype_manual(values = c("Scott"="solid","Doane"="dashed"),name="Bin selection")
p1

legend <- cowplot::get_legend(p1)
grid.arrange(
  p1+theme(legend.position="none"), 
  arrangeGrob(p2+theme(legend.position="none"), legend, ncol = 2, widths = c(3, 1)),
  ncol = 2,
  widths = c(1, 1.35)  # Adjust these as needed to balance plot widths
)


#We see also the development of RSD as a function of rho
rho_seq<-seq(0,10,length.out=20)
rsd_data<-data.frame("rho"=sort(rep(rho_seq,100)),"rsd"=NA)
pb <- txtProgressBar(min = 0, max = nrow(rsd_data), style = 3)
for(i in 1:nrow(rsd_data)){
  setTxtProgressBar(pb, i)
  simulation<-gen_asym(n=5000,rho=rsd_data$rho[i])
  x<-simulation$x
  y<-simulation$y
  kx<-seq(min(x),max(x),length.out=4+ceiling(rsd_data$rho[i]))
  ky<-seq(min(y),max(y),length.out=4+ceiling(rsd_data$rho[i]))
  biv<-zbSpline2D(x,y,alfa=1,bin_selection=scott,knots_x_inner=kx,knots_y_inner=kx,k=2,l=2,u=1,v=1,res=100)
  rsd_data$rsd[i]<-biv$rsd
}


stats<- rsd_data%>%group_by(rho)%>%summarise("median_rsd"=median(rsd),
                                             "LQR"=quantile(rsd,0.025),
                                             "UQR"=quantile(rsd,0.75))


ggplot(stats,mapping=aes(x=rho))+
  geom_line(mapping=aes(y=median_rsd))+
  geom_ribbon(mapping=aes(ymin=LQR,ymax=UQR),alpha=0.2)+
  labs(title="RSD as a function of dependency parameter for perimeter data",
       x=expression(rho),
       y="median RSD")+
  ylim(c(0,1))+
  theme_minimal()
  
  library(microbenchmark)
  library(dHSIC)
  library(copula)
#And we make a study of the computational time. For this it is necessary to have the histogram generation as part of the chisquare computation
gen_hist<-function(x,y,bin_selection=doane){
  
  m<-bin_selection(x)
  n<-bin_selection(y)
  
  breaks_x<-seq(min(x),max(x),length.out=m+1)
  breaks_y<-seq(min(y),max(y),length.out=n+1)
  levels_x <- 1:m
  levels_y <- 1:n
  
  # Binning
  x_bins <- cut(x, breaks = breaks_x, include.lowest = TRUE, right = FALSE, labels = FALSE)
  y_bins <- cut(y, breaks = breaks_y, include.lowest = TRUE, right = FALSE, labels = FALSE)
  
  midpoints_x <- 0.5 * (breaks_x[-1] + breaks_x[-length(breaks_x)])
  midpoints_y <- 0.5 * (breaks_y[-1] + breaks_y[-length(breaks_y)])
  
  # Factors with levels
  x_bins_factor <- factor(x_bins, levels = levels_x)
  y_bins_factor <- factor(y_bins, levels = levels_y)
  
  # Counts
  counts <- table(x_bins_factor, y_bins_factor) #So each row is an observation from x
  counts<-matrix(counts,nrow=m,ncol=n)
  
  #Impute the counts
  
  hist_data<-impute_zeros(counts)
  return(hist_data)
}
N<-c(50,100,500,1000,2000,3000)
run_times<-data.frame("N"=N,"rsd_scott"=NA,"rsd_doane"=NA,"spearman"=NA,"dhsic"=NA,"copula"=NA,"chisq_scott"=NA,"chisq_doane"=NA)
for(i in 1:length(N)){
  print(N[i])
  samp<-gen_sym(n=N[i],rho=0.1)
  x<-samp$x
  y<-samp$y
  kx<-seq(min(x),max(x),length.out=5)
  ky<-seq(min(y),max(y),length.out=5)
  b<-microbenchmark(perm_test(x,y,kx=kx,ky=ky,bin_selection=scott,K=1000),
                    perm_test(x,y,kx=kx,ky=ky,bin_selection=doane,K=1000),
                    cor.test(x,y,method="spearman"),
                    dhsic.test(x,y),
                    indepTest(samp,indepTestSim(N[i],2,N=1000)),
                    chisq.test(gen_hist(x,y,scott),simulate.p.value=TRUE,B=1000),
                    chisq.test(gen_hist(x,y,doane),simulate.p.value=TRUE,B=1000),times=1)
  run_times[i,2:8]<-as.data.frame(summary(b,unit="s"))$median
} #Do note that the hsic test is used with an asymptotic result, so for lower sample sizes we would have to use a permutation test.
#This would make hsic slower at this stage.
#The constant runtime of rsd is because the most computationally demanding aspect is the imputation procedure, and this is barely impacted by sample size

saveRDS(run_times,"run_times_3000")

run_times<-readRDS("run_times_3000")

run_data <- run_times %>%
  pivot_longer(cols = -N, values_to = "run_time", names_to = "test")


p_run<-ggplot(run_data) +
  geom_line(mapping = aes(
    x = N,
    y = run_time,
    color = case_when(
      test %in% c("rsd_scott", "rsd_doane") ~ "rsd",
      test %in% c("chisq_scott", "chisq_doane") ~ "chisq",
      TRUE ~ test
    ),
    linetype = case_when(
      test %in% c("rsd_scott", "chisq_scott") ~ "Scott",
      test %in% c("rsd_doane", "chisq_doane") ~ "Doane",
      TRUE ~ "Scott"
    )
  ), linewidth = 1) +
  ylab("Run time (seconds, sqrt-scale)") +
  xlab("sample size") +scale_y_sqrt()+
  theme_minimal() +
  scale_color_brewer(palette = "Set1",name="Test") +
  scale_linetype_manual(values = c("Scott"="solid","Doane"="dashed"),name="Bin selection")
p_run



#Since this shows us that the runtime is much lower for rsd, but its power grows slower as a function of rho,
#it is interesting to see how the power grows as a function of sample size
p_vals_samp<-data.frame("N"=sort(rep(c(50,100,500,1000,2000,3000),100)),"rsd_scott"=NA,"rsd_doane"=NA,"spearman"=NA,"dhsic"=NA,"copula"=NA,"chisq_scott"=NA,"chisq_doane"=NA)
unique_Ns <- unique(p_vals_samp$N)
ind_sim_list <- list()
for (N in unique_Ns) {
  ind_sim_list[[as.character(N)]] <- indepTestSim(N, 2, N = 1000)
}

# If you want to save all x and y values, initialize lists to store them
x_list <- vector("list", nrow(p_vals_samp))
y_list <- vector("list", nrow(p_vals_samp))

pb <- txtProgressBar(min = 0, max = nrow(p_vals_samp), style = 3)
# Start the main loop
for (i in 1:nrow(p_vals_samp)) { #Skal lige køres forfra :(
  setTxtProgressBar(pb, i)
  
  success <- FALSE
  attempt <- 0
  max_attempts <- 10  # Set a maximum number of attempts
  
  # Declare x and y before the while loop
  x <- NULL
  y <- NULL
  
  while (!success && attempt < max_attempts) {
    attempt <- attempt + 1
    
    # Generate simulation data
    simulation <- gen_sym(n = p_vals_samp$N[i], rho = 0.2) #We use rho=0.2 and the symmetric
    #To know that all tests will eventually reach full power, also Spearman
    x <- simulation$x
    y <- simulation$y
    
    # Compute knots
    kx <- seq(min(x), max(x), length.out = 4)
    ky <- seq(min(y), max(y), length.out = 4)
    
    # Perform permutation test with tryCatch
    # Create objects for both bin_selection = doane and scott
    biv_doane <- tryCatch(
      {
        perm_test(
          x, y, kx, ky, bin_selection = doane,
          k = 2, l = 2, u = 1, v = 1, res = 100, K = 1000
        )
      },
      error = function(e) {
        if (grepl("system is computationally singular", e$message)) {
          # Error occurred, simulate new variables and retry
          return(NULL)
        } else {
          # Re-throw other errors
          stop(e)
        }
      }
    )
    
    biv_scott <- tryCatch(
      {
        perm_test(
          x, y, kx, ky, bin_selection = scott,
          k = 2, l = 2, u = 1, v = 1, res = 100, K = 1000
        )
      },
      error = function(e) {
        if (grepl("system is computationally singular", e$message)) {
          # Error occurred, simulate new variables and retry
          return(NULL)
        } else {
          # Re-throw other errors
          stop(e)
        }
      }
    )
    
    # Process results for biv_doane
    if (!is.null(biv_doane)) {
      p_vals_samp$rsd_doane[i] <- biv_doane$p_val
      success <- TRUE
    } else if (attempt == max_attempts) {
      warning(paste("Maximum attempts reached for doane at i =", i, "with N =", p_vals_samp$N[i]))
      p_vals_samp$rsd_doane[i] <- NA
      success <- TRUE  # Exit loop after maximum attempts
    }
    
    # Process results for biv_scott
    if (!is.null(biv_scott)) {
      p_vals_samp$rsd_scott[i] <- biv_scott$p_val
      success <- TRUE
    } else if (attempt == max_attempts) {
      warning(paste("Maximum attempts reached for scott at i =", i, "with N =", p_vals_samp$N[i]))
      p_vals_samp$rsd_scott[i] <- NA
      success <- TRUE  # Exit loop after maximum attempts
    }
  }  
  
  # Alternatively, store x and y in lists
  x_list[[i]] <- x
  y_list[[i]] <- y
  
  # Proceed with the rest of the tests using the latest x and y
  # Perform Spearman's test
  p_vals_samp$spearman[i] <- cor.test(x_list[[i]], y_list[[i]], method = "spearman")$p.value
  
  # Perform DHSIC test
  p_vals_samp$dhsic[i] <- dhsic.test(x_list[[i]], y_list[[i]], method = "permutation")$p.value
  
  # Retrieve the precomputed ind_sim
  ind_sim <- ind_sim_list[[as.character(p_vals_samp$N[i])]]
  
  # Perform the copula independence test
  p_vals_samp$copula[i] <- indepTest(data.frame(x_list[[i]],y_list[[i]]), ind_sim)$pvalues
  
  p_vals_samp$chisq_doane[i]<- chisq.test(gen_hist(x_list[[i]],y_list[[i]],doane),simulate.p.value=TRUE,B=1000)$p.value
  
  p_vals_samp$chisq_scott[i]<-chisq.test(gen_hist(x_list[[i]],y_list[[i]],scott),simulate.p.value=TRUE,B=1000)$p.value
}


p_vals_power_N<-p_vals_samp %>%
  pivot_longer(cols = -N, names_to = "test", values_to = "p_value") %>%
  group_by(N, test) %>%
  summarise(
    power=length(which(p_value<0.05))
  )

saveRDS(p_vals_samp,"p_vals_samp")

p_vals_samp<-readRDS("p_vals_samp")

p_pow<-ggplot(p_vals_power_N) +
  geom_line(mapping = aes(
    x = N,
    y = power,
    color = case_when(
      test %in% c("rsd_scott", "rsd_doane") ~ "rsd",
      test %in% c("chisq_scott", "chisq_doane") ~ "chisq",
      TRUE ~ test
    ),
    linetype = case_when(
      test %in% c("rsd_scott", "chisq_scott") ~ "Scott",
      test %in% c("rsd_doane", "chisq_doane") ~ "Doane",
      TRUE ~ "Scott"
    )
  ), linewidth = 1) +
  ylab("Power (%)") +
  xlab("sample size") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1",name="Test") +
  scale_linetype_manual(values = c("Scott"="solid","Doane"="dashed"),name="Bin selection")


#We can combine the two plots to see power as a function of computational time
get_legend(p_run)
ggarrange(p_run,p_pow,ncol=1,common.legend=TRUE,legend="right")
grid.arrange(p_run,p_pow)


run_times
comb_datas<-merge(run_data,p_vals_power_N,by=c("test","N"),all.x=FALSE,all.y=FALSE)
ggplot(comb_datas,aes(x=run_time,color=test))+
  geom_line(aes(y=power))+
  theme_minimal()+
  labs(title="Centered data",
       x="run time (seconds)",
       y="power (%)")+
  scale_color_brewer(palette="Set1")+ylim(c(0,100))+xlim(c(0,10))
