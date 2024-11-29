library(plot3D)

sim2<-function(n,a=1){
  y<-runif(n,min=-pi,max=pi) #Generate y marginally
  U<-runif(n,min=0,max=1) #Used for quantile method
  x<-numeric(n) #Will fill this out one at a time
  for(i in 1:n){
    f<-function(x){(sin(a*y[i])*(cos(a*x)-cos(a*pi))+a*(x+pi))/(2*a*pi)-U[i]}
    #This is the conditional distribution function (-U)
    x[i]<-uniroot(f,interval=c(-pi,pi))$root
  }
  return(data.frame("X"=x,"Y"=y))
}
sim3<-function(n,a){
  Theta<-runif(n,min=0,max=2*pi)
  epsilon1<-rnorm(n)
  epsilon2<-rnorm(n)
  X<-a*cos(Theta)+epsilon1/4
  Y<-a*sin(Theta)+epsilon2/4
  return(data.frame("X"=X,"Y"=Y))
}
sim4<-function(n,a){
  X<-runif(n,min=-1,max=1)
  Y<-abs(X)ˆa*rnorm(n)
  return(data.frame("X"=X,"Y"=Y))
}

data<-sim3(1000,a=1)

#We start with making the compositional spline for the univariate case.
library(splines)
library(splines2)
library(e1071) #Skewness
library(tidyverse)
library(ggpubr)
library(LaplacesDemon)

k<-3 #degree (order=4)
data<-data$X

knots<-function(data,k=3,df=10,knots=NA){
   interval<-range(data)
  lead_knots<-rep(interval[1],k)
  end_knots<-rep(interval[2],k)
  if(is.na(knots[1])){
    mid<-seq(from=interval[1],to=interval[2],length.out=(df+2))
  }
  else{
    mid<-knots
  }
  return(c(lead_knots,mid,end_knots))
}

sturge<-function(data){
  ceiling(1+log2(length(data)))
}
scott<-function(data){
  ceiling(1/(3.5*sd(data)*length(data)^(-1/3)))
}
doane<-function(data){
  n <- length(data)
  g1 <- skewness(data)
  sigma_g1 <- sqrt(6 * (n - 2) / ((n + 1) * (n + 3)))
  k <- 1 + log2(n) + log2(1 + abs(g1) / sigma_g1)
  return(ceiling(k))
}



#We start with replicating example 1 before I extend to smoothing (choice of z)
x<-seq(from=0,to=20,length.out=1000)

g<-4

k<-3 #degree

knot_points<-gen_knots(x,3,knots=c(0,2,5,9,14,20))

B<-splineDesign(knots,x,ord=(k+1),outer.ok=TRUE,derivs=0)


#We can plot it just cus
B_df <- as.data.frame(B)
B_df$x <- x

B_long <- pivot_longer(B_df, cols = -x, names_to = "basis", values_to = "value")

ggplot(B_long, aes(x = x, y = value, color = basis)) +
  geom_line()+theme(legend.position="none")


l<-numeric(g+k+1)
for(i in 1:length(l)){
  l[i]<-1/(knot_points[i+k+1]-knot_points[i]) #Siden lambda[1]=lambda_-k er lambda[k+1]=lambda_1
}
D<-(k+1)*diag(l)

generate_K <- function(g, k) {
  # Initialize the matrix with zeros
  rows <- g + k + 1
  cols <- g + k
  mat <- matrix(0, nrow = rows, ncol = cols)
  
  # Set the first row (leading 1)
  mat[1, 1] <- 1
  
  # Fill the rest of the rows with the -1, 1 pattern
  for (i in 2:rows) {
    if (i - 1 <= cols) {
      mat[i, i - 1] <- -1
    }
    if (i <= cols) {
      mat[i, i] <- 1
    }
  }
  return(mat)
}

K<-generate_K(g,k)

Z<-B%*%D%*%K
Z_df<-as.data.frame(Z)
Z_df$x<-x

Z_long<-pivot_longer(Z_df,cols=-x,names_to="basis",values_to="value")

ggplot()+
#  geom_histogram(data=Z_long,mapping=aes(x=x,y=after_stat(density)),bins=g)+
  geom_line(Z_long,mapping=aes(x=x,y=value,color=basis),linewidth=1.2)+
  theme(legend.position="none")


#We start smoothing

w<-rep(1,length(x)) #Man kan jo tage andre weights, men jeg kan ikke lige se hvorfor

W<-diag(w)

U<-D%*%K

Dj<-function(j,knot_points){
  d<-numeric(g+k+1-j)
  for(i in seq(from=(-k+j),to=g,by=1)){
    d[i]<-1/(knot_points[(i+k+1-j+k+1)]-knot_points[i+k+1]) #k+1 for at vi starter i lambda_1
  }
  return(diag(d)*(k+1-j))
}
Lj<-function(j){
  L_empty<-matrix(0,nrow=(g+k+1-j),ncol=(g+k+2-j))
  for(i in 1:nrow(L_empty)){
    for(j in 1:ncol(L_empty)){
      if(i==j){
        L_empty[i,j]<- -1
      }
      if(i==(j-1)){
        L_empty[i,j]<- 1
      }
    }
  }
  return(L_empty)
}

Sl<-function(l,knot_points){
  prod<-Dj(1,knot_points)%*%Lj(1)
  for(i in 2:l){
    prod<-Dj(i,knot_points)%*%Lj(i)%*%prod
  }
  return(prod)
}
Sl(2,knot_points)

l<-2
M_func <- function(l) {
  M <- matrix(NA, nrow = (g + k + 1-l), ncol = (g + k + 1 - l))
  a<-min(knot_points)
  b<-max(knot_points)
  f<-function(i,j,a){
    splineDesign(knot_points,a,ord=k+1-l)[i]*splineDesign(knot_points,a,ord=k+1-l)[j]
  }
  for (i in 1:nrow(M)) {
    for (j in 1:ncol(M)) {
      M[i, j] <- ((b-a)/6)*(f(i,j,a)+4*f(i,j,(b-a)/2)+f(i,j,b)) #Approximate the integral by Simpson's rule
    }
  }
  return(M)
}

alpha<-0.5

G<-t(U)%*%((1-alpha)*t(Sl(l,knot_points))%*%M_func(l)%*%Sl(l,knot_points)+alpha*t(B)%*%W%*%B)%*%U

classes<- cut(x,breaks=g)
y<-table(classes)
y<-y[classes]/length(x) #Bør nok også clr-transformeres

g_mat<-alpha*t(K)%*%D%*%t(B)%*%W%*%y

z_star<-solve(G)%*%g_mat

Spline<-B%*%D%*%K%*%z_star
Spline_df<-as.data.frame(Spline)
Spline_df$x<-x

Spline_long<-pivot_longer(Spline_df,cols=-x,names_to="basis",values_to="value")

ggplot(Spline_long)+geom_line(mapping=aes(x=x,y=value,color=basis))+theme(legend.position="nonte")

#Vi skal lige finde ud af noget med at lave de der y'ere og clr transformere dem