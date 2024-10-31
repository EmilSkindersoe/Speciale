setwd("~/Universitet-mat/.Speciale/R")
source("splines_source1.R")
source("splines_bivariate.R")
library(gridExtra)
library(pracma)

#Tjekker at den virker ved at forsøge og genskabe figur 5 og 6 fra mechalova
set.seed(3)
x_norm<-rnorm(20000)

Z_fig5<-zbSpline(x_norm,degree=2,l=1,alfa=0.9,knots_inner=c(min(x_norm),-2,-1,0,1,2,max(x_norm)),res=1000,bin_selection=33)
#We add the clr-density of the standard normal
Z_fig5<-cbSpline(Z_fig5,what="original")
plot(Z_fig5,what="C-spline",include_hist=TRUE)
f<-function(x){dnorm(x)}
clr_dens<-clr(f,x_norm)
x_seq<-seq(min(x_norm),max(x_norm),length.out=length(x_norm))
dens_data<-data.frame("x"=x_norm,"clr_density"=clr_dens(x_norm),"density"=dnorm(x_norm))


p3<-plot(Z_fig5,what="spline",include_hist=TRUE)+
  geom_hline(mapping=aes(yintercept=0),linetype="dotted")+
  geom_line(dens_data,mapping=aes(x=x,y=clr_density))
Z_C<-cbSpline(Z_fig5,what="original")
p4<-plot(Z_C,what="C-spline")+geom_line(dens_data,mapping=aes(x=x,y=density))

row1_title <- textGrob("N=1000, n=33", gp = gpar(fontsize = 10, fontface = "bold"))
row2_title <- textGrob("N=20000, n=33", gp = gpar(fontsize = 10, fontface = "bold"))

row1 <- arrangeGrob(
  row1_title, 
  arrangeGrob(p1, p2, ncol = 2), 
  nrow = 2,  # 2 rows: one for the title, one for the plots
  heights = c(0.2, 4)  # Adjust heights
)

# Arrange the second row (title + plots)
row2 <- arrangeGrob(
  row2_title, 
  arrangeGrob(p3, p4, ncol = 2), 
  nrow = 2,  # 2 rows: one for the title, one for the plots
  heights = c(0.2, 4)  # Adjust heights
)

# Combine everything into a final layout
grid.arrange(row1, row2, nrow = 2)


plot(Z_fig5,what="basis")

Z_O<-orthogonalise(Z_fig5,dec="svd")

plot(Z_O,what="O-basis")
plot(Z_O,what="O-spline")



plot(Z_C,what="C-basis")



cross_validation<-cross_validate(x_norm,k=2,l=1,knots_inner=c(min(x_norm),-2,-1,0,1,2,max(x_norm))) #Skal lige fikses
cross_validation

plot(cross_validation)

#With the same breaks we get the same behaviour
breaks<-seq(min(x_norm),max(x_norm),length.out=doane(x_norm)+1)
hist_info<-hist(x_norm,breaks=breaks,plot=FALSE,include.lowest = TRUE, right = FALSE)
hist_data<-data.frame("x"=hist_info$mids,"y"=hist_info$density)
debug(zbSpline)
Z_fig5<-zbSpline(x=hist_data,degree=3,l=2,alfa=0.9,knots=NULL,bin=scott)
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

breaks<-seq(min(x),max(x),length.out=12)
hist_info<-hist(x,breaks=breaks,plot=FALSE,include.lowest = TRUE, right = FALSE)
hist_data<-data.frame("x"=hist_info$mids,"y"=hist_info$counts)
Z<-zbSpline(hist_data,alfa=0.9,degree=3,l=2,knots_inner=seq(min(x),max(x),length.out=4))
Z<-zbSpline(x,alfa=0.9,degree=3,l=2,knots_inner=c(0,0.5,1,2,4,8,16),bin_selection = 50)
O<-orthogonalise(Z)
grid.arrange(plot(O,what="basis"),plot(O,what="O-basis"),nrow=1)
plot(C,what="C-spline")+geom_line(dens_data,mapping=aes(x=x,y=density))

C<-cbSpline(Z)
plot(C,what="C-spline",include_hist=TRUE)

#Der er nogle basis-splines der skal fikses lidt
#I dette tilfælde kan det være svært at bruge histogram data
#da vi skal have mange bins for at fange den hurtige vækst af chisquared,
#men så mangler der data mod slutningen

#plot(Z,what="basis")
f<-function(x){dchisq(x,df=3)}
clr_dens<-clr(f,x)
x_seq<-seq(min(x),max(x),length.out=length(x))
dens_data<-data.frame(x,"clr_density"=clr_dens(x),"density"=dchisq(x,df=3))
p<-plot(Z,what="spline")
p+geom_line(dens_data,mapping=aes(x=x,y=clr_density))



#We now try reconstructing figure 5,6 and 7 of skorna
a_0<-3; a_1<-3; a_2<-3
f<-function(x,y){
  normaliser<-gamma(a_0+a_1+a_2)/(gamma(a_0)*gamma(a_1)*gamma(a_2))
  #normaliser<-1
  func_value<-x^(a_1 - 1) * (1 - x)^(a_0 + a_2 - 1) * y^(a_2 - 1) * (1 - y)^(a_0 + a_1 - 1) / (1 - x * y)^(a_0 + a_1 + a_2)
  return(func_value*normaliser)
}

M<-4.1
sample_bivariate<-function(n=3000){
  vals<-data.frame("X"=NA,"Y"=NA)
  i<-1
  while(i<=n){
    X<-runif(1);Y<-runif(1);U<-runif(1)
    if (U<=f(X,Y)/M){
      vals[i,]<-c(X,Y)
      i<-i+1
    }
  }
  return(vals)
}




sample<-sample_bivariate()
x<-sample$X
y<-sample$Y
m<-10
n<-10
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
counts <- table(x_bins_factor, y_bins_factor)
counts<-matrix(counts,nrow=m,ncol=n)
counts<-impute_zeros(counts)

list_data<-list("x"=midpoints_x,"y"=midpoints_y,"data"=counts)
knots_x_inner<-c(0,0.25,0.5,0.75,1)->knots_y_inner

hist3D(midpoints_x,midpoints_y,counts/(sum(counts)*0.1*0.1))

Z2D<-bivariate(list_data,rho=0.001,knots_x_inner=knots_x_inner,knots_y_inner=knots_y_inner,k=2,l=2,u=1,v=1,res=50)
Z2D$rsd
plot(Z2D,what="full",scale="density",plot_hist=TRUE)


simpson2d(f,0.001,0.99,0.001,0.99)
x_vals <- seq(0.01, 0.99, length.out = 50)
y_vals <- seq(0.01, 0.99, length.out = 50)
z_vals <- outer(x_vals, y_vals, function(x, y) f(x, y))
plot_ly(x=x_vals,y=y_vals,z=z_vals,type="surface")


