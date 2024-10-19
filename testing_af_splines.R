setwd("~/Universitet-mat/.Speciale/R")
source("splines_source1.R")
library(gridExtra)
library(pracma)

#Tjekker at den virker ved at forsøge og genskabe figur 5 og 6 fra mechalova
set.seed(3)
x_norm<-rnorm(1000)

Z_fig5<-zbSpline(x_norm,degree=2,l=1,alfa=0.9,knots_inner=c(min(x_norm),-2,-1,0,1,2,max(x_norm)),res=1000)

#We add the clr-density of the standard normal
f<-function(x){dnorm(x)}
clr_dens<-clr(f,x_norm)
x_seq<-seq(min(x_norm),max(x_norm),length.out=length(x_norm))
dens_data<-data.frame("x"=x_norm,"clr_density"=clr_dens(x_norm),"density"=dnorm(x_norm))


plot(Z_fig5,what="spline",include_hist=TRUE)+
  geom_hline(mapping=aes(yintercept=0),linetype="dotted")+
  geom_line(dens_data,mapping=aes(x=x,y=clr_density))



plot(Z_fig5,what="basis")

Z_O<-orthogonalise(Z_fig5,dec="svd")

plot(Z_O,what="O-basis")
plot(Z_O,what="O-spline")

Z_C<-cbSpline(Z_O,what="original")

plot(Z_C,what="C-basis")
plot(Z_C,what="C-spline")+geom_line(dens_data,mapping=aes(x=x,y=density))


cross_validation<-cross_validate(x_norm,k=2,l=1,knots_inner=c(min(x_norm),-2,-1,0,1,2,max(x_norm))) #Skal lige fikses
cross_validation

plot(cross_validation)

#With the same breaks we get the same behaviour
breaks<-seq(min(x_norm),max(x_norm),length.out=doane(x_norm)+1)
hist_info<-hist(x_norm,breaks=breaks,plot=FALSE,include.lowest = TRUE, right = FALSE)
hist_data<-data.frame("x"=hist_info$mids,"y"=hist_info$density)
debug(zbSpline)
Z_fig5<-zbSpline(x=hist_data,degree=3,l=2,alfa=0.9,knot=NULL,bin=scott)
undebug(zbSpline)
#And we try with self-selection of knots

Z_fig5<-zbSpline(x_norm,alfa=0.9,l=2,degree=3,knots=NULL,bi=scott)


plot(Z_fig5,what="spline")

#Skriv et afsnit herfra, hvor vi genskaber figur 1 og 6 fra mechalova

x_fig1<-seq(0,20,length.out=1000)
Z_fig1<-zbSpline(x_fig1,degree=3,l=2,knots_inner=c(0,2,5,9,14,20))
Z_fig1<-orthogonalise(Z_fig1,"svd")
plot(Z_fig1,what="O-basis")

#Og så låner vi plot metoden til at ændre lidt i hvad Z_splinen er
Z_fig1$Z_spline<-Z_fig1$Z_basis%*%c(0.5,-1,2,3,-8,9,1)
plot(Z_fig1,what="spline",include_hist=FALSE)+geom_hline(yintercept=0)
grid.arrange(plot(Z_fig1,what="spline",include_hist=FALSE)+geom_hline(yintercept=0),plot(Z_fig1),nrow=1)


#Genskaber
x_fig3<-seq(0,3,length.out=1000)
z_fig3<-zbSpline(x_fig3,degree=2,l=1,knots_inner=c(0,1,2,3))
o_fig3<-orthogonalise(z_fig3)
o_fig3_svd<-orthogonalise(z_fig3,dec="svd")
grid.arrange(plot(z_fig3,what="basis")+ggtitle("Original ZB-spline"),plot(o_fig3,what="O-basis")+ggtitle("Orthogonal cholesky ZB-spline"),plot(o_fig3_svd,what="O-basis")+ggtitle("Orthogonal SVD ZB-spline"),nrow=1)


set.seed(2)
x<-rchisq(1000,df=3)

breaks<-seq(min(x),max(x),length.out=doane(x)+1)
hist_info<-hist(x,breaks=breaks,plot=FALSE,include.lowest = TRUE, right = FALSE)
hist_data<-data.frame("x"=hist_info$mids,"y"=hist_info$counts)
Z<-zbSpline(x,alfa=0.9,degree=3,l=2,knots_inner=seq(min(x),max(x),length.out=20),bin_selection = 70)
Z<-zbSpline(x,alfa=0.9,degree=3,l=2,knots_inner=c(0,0.5,1,2,4,8,16),bin_selection = 50)
#Der er nogle basis-splines der skal fikses lidt
#I dette tilfælde kan det være svært at bruge histogram data
#da vi skal have mange bins for at fange den hurtige vækst af chisquared,
#men så mangler der data mod slutningen

#plot(Z,what="basis")
f<-function(x){dchisq(x,df=3)}
clr_dens<-clr(f,x)
x_seq<-seq(min(x),max(x),length.out=length(x))
dens_data<-data.frame(x,"clr_density"=clr_dens(x))
p<-plot(Z,what="spline")
p+geom_line(dens_data,mapping=aes(x=x,y=clr_density))


