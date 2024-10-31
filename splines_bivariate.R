setwd("~/Universitet-mat/.Speciale/R")
source("splines_source1.R")
library(gridExtra)
library(pracma)
#We will use simpson2d() for the integration here
library(fastmatrix)
#Useful for the kronecker.prod() function
library(plot3D)
#Useful for static plots
library(plotly)
library(viridis)
#Useful for interactive plots



#The two first functions are for the imputations
get_neighbors <- function(mat, i, j) {
  n_row <- nrow(mat)
  n_col <- ncol(mat)
  
  # Define the range for neighboring indices
  row_range <- (i-1):(i+1)
  col_range <- (j-1):(j+1)
  
  # Keep indices within matrix bounds
  row_indices <- row_range[row_range >= 1 & row_range <= n_row]
  col_indices <- col_range[col_range >= 1 & col_range <= n_col]
  
  # Create all combinations of neighboring positions
  neighbor_positions <- expand.grid(row = row_indices, col = col_indices)
  
  # Exclude the cell itself
  neighbor_positions <- neighbor_positions[!(neighbor_positions$row == i & neighbor_positions$col == j), ]
  
  # Retrieve the values of the neighboring cells
  neighbor_values <- mapply(function(x, y) mat[x, y], neighbor_positions$row, neighbor_positions$col)
  return(neighbor_values)
}

impute_zeros <- function(counts) {
  counts <- as.matrix(counts)
  
  # Continue the process until there are no zeros left
  while(any(counts == 0)) {
    counts_new <- counts
    zero_positions <- which(counts == 0, arr.ind = TRUE)
    
    # Iterate over each zero entry
    for(k in 1:nrow(zero_positions)) {
      i <- zero_positions[k, 1]
      j <- zero_positions[k, 2]
      neighbor_values <- get_neighbors(counts, i, j)
      non_zero_neighbors <- neighbor_values[neighbor_values > 0]
      
      # If there are non-zero neighbors, compute the geometric mean
      if(length(non_zero_neighbors) > 0) {
        gm <- exp(mean(log(non_zero_neighbors)))*2/3
        counts_new[i, j] <- gm
      }
      # Else, leave it as zero for now
    }
    counts <- counts_new/sum(counts_new)
  }
  return(counts)
}

trapz2d <- function(x, y, f) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("The 'pracma' package is required but not installed.")
  }
  nx <- length(x)
  ny <- length(y)
  if (!is.matrix(f) || dim(f)[1] != nx || dim(f)[2] != ny) {
    stop("Dimensions of 'f' must match lengths of 'x' and 'y'")
  }
  
  #Using Fubini, we first integrate over x to make outcomes depending on the value of y
  #And then integrate over those values of y
  # Integrate f over x for each y
  int_over_x <- apply(f, 2, function(fy) pracma::trapz(x, fy))
  # Integrate the result over y
  integral <- pracma::trapz(y, int_over_x)
  
  return(integral)
}

#Denne funktion mangler mulighed for selv at finde knots og mulighed for at specificere antal bins
bivariate<-function(x,y=NULL,alfa=0.5,rho=NULL,bin_selection=doane,knots_x_inner=NULL,knots_y_inner=NULL,k=3,l=3,u=2,v=2,res=200){
  if(is.null(rho)){
    rho=(1-alfa)/alfa
  }
  if(class(x)=="list"){ #Assumes histogram data is given as list of length 3 with the first two being the bins and the third being a matrix of bin values 
    midpoints_x<-x[[1]]
    midpoints_y<-x[[2]]
    m<-length(midpoints_x)
    n<-length(midpoints_y)
    hist_data<-x[[3]]
    if(any(hist_data==0)){
      warning("Some hist_data were 0, imputing...")
      hist_data<-impute_zeros(hist_data)
    }
    F_mat<-CLR(hist_data)
  } else if (class(x)=="numeric" & class(y)=="numeric"){ #Converts to histogram data
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
    counts <- table(x_bins_factor, y_bins_factor)
    counts<-matrix(counts,nrow=m,ncol=n)/sum(counts)
    
    #Impute the counts
    hist_data<-impute_zeros(counts)
    F_mat<-CLR(hist_data)
  }
  else {stop("Error: Data must be either histogram data x or numerics (x,y).")}
  
  #We make the knot sequence
  #We start with assuming that the inner knots are provided
  
  g<-length(knots_x_inner)-2
  h<-length(knots_y_inner)-2
  knots_x<-c(rep(knots_x_inner[1],k),knots_x_inner,rep(knots_x_inner[length(knots_x_inner)],k))
  knots_y<-c(rep(knots_y_inner[1],l),knots_y_inner,rep(knots_y_inner[length(knots_y_inner)],l))
  
  
  a<-min(knots_x); b<-max(knots_x); c<-min(knots_y); d<-max(knots_y)
  
  seq_x<-seq(a,b,length.out=res)
  seq_y<-seq(c,d,length.out=res)
  
  #We generate the Z_x_bar and Z_y_bar matrices
  lambda_x<-numeric(g+k+1)
  for(i in 1:length(lambda_x)){
    lambda_x[i]<-1/(knots_x[i+k+1]-knots_x[i]) #Siden lambda[1]=lambda_-k er lambda[k+1]=lambda_1
  }
  D_x<-(k+1)*diag(lambda_x)
  
  K_x<-matrix(0,nrow=g+k+1,ncol=g+k)
  K_x[1,1]<-1
  for (i in 2:nrow(K_x)) {
    if (i - 1 <= ncol(K_x)) {
      K_x[i, i - 1] <- -1
    }
    if (i <= ncol(K_x)) {
      K_x[i, i] <- 1
    }
  }
  
  lambda_y<-numeric(h+l+1)
  for(i in 1:length(lambda_y)){
    lambda_y[i]<-1/(knots_y[i+l+1]-knots_y[i]) #Siden lambda[1]=lambda_-k er lambda[k+1]=lambda_1
  }
  D_y<-(l+1)*diag(lambda_y)
  
  K_y<-matrix(0,nrow=h+l+1,ncol=h+l)
  K_y[1,1]<-1
  for (i in 2:nrow(K_y)) {
    if (i - 1 <= ncol(K_y)) {
      K_y[i, i - 1] <- -1
    }
    if (i <= ncol(K_y)) {
      K_y[i, i] <- 1
    }
  }
  
#We create the marginal ZB-splines
#They are evaluated in the midpoints for optimisation and along seq_x,seq_y
  Z_x<-splineDesign(knots_x,midpoints_x,ord=k+1,outer.ok=TRUE)%*%D_x%*%K_x
  Z_x_seq<-splineDesign(knots_x,seq_x,ord=k+1,outer.ok=TRUE)%*%D_x%*%K_x
  Z_y<-splineDesign(knots_y,midpoints_y,ord=l+1,outer.ok=TRUE)%*%D_y%*%K_y
  Z_y_seq<-splineDesign(knots_y,seq_y,ord=l+1,outer.ok=TRUE)%*%D_y%*%K_y

#We create the bar matrices used for expressing the Z-spline
  Z_x_bar<-rbind(t(Z_x),1)
  Z_x_bar_seq<-rbind(t(Z_x_seq),1)
  Z_y_bar<-rbind(t(Z_y),1)
  Z_y_bar_seq<-rbind(t(Z_y_seq),1)
  
#And we create the black board Z matrices which are useful for optimisation
  bbZ<-kronecker(Z_y,Z_x)
  bbZ_x<-kronecker(rep(1,n),Z_x)
  bbZ_y<-kronecker(Z_y,rep(1,m))
  bbZ_bar<-t(cbind(bbZ,bbZ_x,bbZ_y))

#Now we just need the bbN matrix to express G(rho)
#First the matrix S
  Dj_x<-function(j){
    d<-numeric(g+k+1-j)
    for(i in seq_along(d)){
      d[i]<-knots_x[i+k+1]-knots_x[i+j]
    }
    return(diag(1/d)*(k+1-j))
  }
  
  Lj_x <- function(j) {
    M <- g + k + 1 - j
    L_empty <- matrix(0, nrow = M, ncol = M + 1)
    L_empty[cbind(1:M, 1:M)] <- -1
    L_empty[cbind(1:M, 2:(M + 1))] <- 1
    return(L_empty)
  }
  # Generate matrix Su
  if(u==0){ Su<-diag(1,g+k+1)}
  Su <- Dj_x(1) %*% Lj_x(1)
  if (u > 1) for (i in 2:u) Su <- Dj_x(i) %*% Lj_x(i) %*% Su
  
  Dj_y<-function(j){
    d<-numeric(h+l+1-j)
    for(i in seq_along(d)){
      d[i]<-knots_y[i+l+1]-knots_y[i+j]
    }
    return(diag(1/d)*(l+1-j))
  }
  
  Lj_y <- function(j) {
    N <- h + l + 1 - j
    L_empty <- matrix(0, nrow = N, ncol = N + 1)
    L_empty[cbind(1:N, 1:N)] <- -1
    L_empty[cbind(1:N, 2:(N + 1))] <- 1
    return(L_empty)
  }
  # Generate matrix Sv
  if(v==0){ Sv<-diag(1,h+l+1)}
  Sv <- Dj_y(1) %*% Lj_y(1)
  if (v > 1) for (i in 2:v) Sv <- Dj_y(i) %*% Lj_y(i) %*% Sv
  
  
  bbS<-kronecker(Sv,Su)
  bbD<-kronecker(D_y,D_x)
  bbK<-kronecker(t(K_y),t(K_x)) #We use a transposed version of K
  
  #And we make the matrix M
  M_x <- matrix(NA, nrow = (g + k + 1-u), ncol = (g + k + 1 - u))
  knots_M_x<-c(rep(knots_x_inner[1],k-u),knots_x_inner,rep(knots_x_inner[length(knots_x_inner)],k-u))
  f<-function(i,j,x){
    mat<-splineDesign(knots_M_x,x,ord=k+1-u,outer.ok=TRUE)
    return(mat[,i]*mat[,j])
  }
  for (i in 1:nrow(M_x)) {
    for (j in 1:ncol(M_x)) {
      M_x[i,j]<- trapz(seq_x,f(i,j,seq_x))
    }
  }
  
  M_y <- matrix(NA, nrow = (h + l + 1-v), ncol = (h + l + 1 - v))
  knots_M_y<-c(rep(knots_y_inner[1],l-v),knots_y_inner,rep(knots_y_inner[length(knots_y_inner)],l-v))
  f<-function(i,j,x){
    mat<-splineDesign(knots_M_y,x,ord=l+1-v,outer.ok=TRUE)
    return(mat[,i]*mat[,j])
  }
  for (i in 1:nrow(M_y)) {
    for (j in 1:ncol(M_y)) {
      M_y[i,j]<- trapz(seq_y,f(i,j,seq_y))
    }
  }
  
  bbM<-kronecker(M_y,M_x)
  bbN<-bbK%*%bbD%*%t(bbS)%*%bbM%*%bbS%*%bbD%*%t(bbK)
  
  col1<-rbind(rho*bbN+t(bbZ)%*%bbZ,t(bbZ_x)%*%bbZ,t(bbZ_y)%*%bbZ)
  col2<-rbind(t(bbZ)%*%bbZ_x,t(bbZ_x)%*%bbZ_x,t(bbZ_y)%*%bbZ_x)
  col3<-rbind(t(bbZ)%*%bbZ_y,t(bbZ_x)%*%bbZ_y,t(bbZ_y)%*%bbZ_y)
  G<-cbind(col1,col2,col3)
  
  #Hmm G er ikke lige invertibel :(
  coefficients<-solve(G)%*%bbZ_bar%*%as.vector(F_mat)
  Z_opt<-matrix(coefficients[1:((g+k)*(h+l))],nrow=k+g,ncol=l+h)
  v_opt<-coefficients[((g+k)*(h+l)+1):((g+k)*(h+l)+k+g)]
  u_opt<-coefficients[((g+k)*(h+l)+k+g+1):((g+k)*(h+l)+k+g+h+l)]
  # R_opt<-rbind(Z_opt,t(u_opt))
  # R_opt<-cbind(R_opt,c(v_opt,0))
  # spline_opt<-t(Z_x_bar_seq)%*%R_opt%*%Z_y_bar_seq
  spline_int<- Z_x_seq%*%Z_opt%*%t(Z_y_seq)
  spline_x<-Z_x_seq%*%v_opt
  spline_y<-Z_y_seq%*%u_opt
  spline_ind<-outer(spline_x,spline_y,FUN="+")
  spline_ind<-matrix(spline_ind,nrow=res,ncol=res)
  spline_opt<-spline_int+spline_ind
  
  #Calculate relative simplical deviance
  simp_d<-trapz2d(seq_x,seq_y,spline_int^2)
  rel_simp_d<-simp_d/trapz2d(seq_x,seq_y,spline_opt^2)
  
  
  #And we map it back to density
  C_spline<-exp(spline_opt)
  C_spline<-C_spline/trapz2d(seq_x,seq_y,C_spline)
  C_spline_ind<-exp(spline_ind)
  C_spline_ind<-C_spline_ind/trapz2d(seq_x,seq_y,C_spline_ind)
  C_spline_int<-exp(spline_int)
  C_spline_int<-C_spline_int/trapz2d(seq_x,seq_y,C_spline_int)
  
  structure(list("Z_spline"=spline_opt,
                 "Z_spline_int"=spline_int,
                 "Z_spline_ind"=spline_ind,
                 "Z_spline_x"=spline_x,
                 "Z_spline_y"=spline_y,
                 "C_spline"=C_spline,
                 "C_spline_ind"=C_spline_ind,
                 "C_spline_int"=C_spline_int,
                 "seq_x"=seq_x,"seq_y"=seq_y,
                 "midpoints_x"=midpoints_x,"midpoints_y"=midpoints_y,
                 "F_mat"=F_mat,"hist_data"=hist_data,
                 "rho"=rho,"alfa"=alfa,
                 "sd"=simp_d,"rsd"=rel_simp_d),class="bivariate_zbSpline")
}





plot.bivariate_zbSpline<-function(Z,type="static",what="full",scale="clr",plot_hist=FALSE){
  #type can be static or interactive
  #What can be either full, independent, interaction or geometric marginals
  #plot_hist is whether we plot the underlying histogram
  #scale chooses between the clr transform or the density
  if(type=="static"){ #Here we use the persp3D function
    #We can update all of these to finding a response and range based on the if statements
    if(what=="full"){
      if(scale=="clr"){
        outcome<-Z$Z_spline
      }
      if(scale=="density"){
        outcome<-Z$C_spline
      }
    }
    if(what=="interaction"){
      if(scale=="clr"){
        outcome<-Z$Z_spline_int
      }
      if(scale=="density"){
        outcome<-Z$C_spline_int
      }
    }
    if(what=="independent"){
      if(scale=="clr"){
        outcome<-Z$Z_spline_ind
      }
      if(scale=="density"){
        outcome<-Z$C_spline_ind
      }
    }
    range<-c(min(outcome)-0.5,max(outcome)+0.5)
    
      persp3D(
        x = Z$seq_x, 
        y = Z$seq_y, 
        z = outcome, 
        col = viridis(50),
        theta = 315,          
        phi = 30,            
        ticktype = "detailed", 
        nticks = 3,
        zlim = range,
        xlab = "",          
        ylab = "",          
        zlab = "",          
        bty = "b2"          
      )
    
    if(plot_hist){
      if(scale=="clr"){
        hist_response<-Z$F_mat
      }
      if(scale=="density"){
        hist_response<-Z$hist_data/trapz2d(Z$midpoints_x,Z$midpoints_y,Z$hist_data)
      }
      points3D(
        x = rep(Z$midpoints_x, each = length(Z$midpoints_y)),   # Repeat x for each y value
        y = rep(Z$midpoints_y, times = length(Z$midpoints_x)),   # Repeat y for each x value
        z = as.vector(hist_response),                # Flatten z matrix to match x and y
        add = TRUE,
        col = "black",
        pch = 19,
        cex = 0.5
      )
    }
  }
}

#We still need a cross_validate function and a bivariate_cbSpline


# plot(Z2D,what="full",plot_hist=TRUE)
# 
# Z<-bivariate(x,y,knots_x_inner=knots_x_inner,knots_y_inner=knots_y_inner)
# 
# all.equal(spline_2d$spline,spline_2d$spline_int+spline_2d$spline_ind)
# 
# 
# 
# persp3D(seq_x,seq_y,spline_ind+spline_int,col="green")
# persp3D(seq_x,seq_y,spline_int,add=TRUE,col="red")
# persp3D(seq_x,seq_y,spline_ind,add=TRUE,col="blue")
# 
# p <- plot_ly()
# 
# # Add the combined surface with a solid color
# p <- p %>% add_surface(x = ~seq_x, y = ~seq_y, z = ~spline_opt, 
#                        colorscale = list(c(0, 1), c("blue", "blue")), 
#                        opacity = 0.6, name = "Combined Surface",showscale=FALSE)
# 
# # Add the spline_int surface with a solid color
# p <- p %>% add_surface(x = ~seq_x, y = ~seq_y, z = ~spline_int, 
#                        colorscale = list(c(0, 1), c("green", "green")), 
#                        opacity = 0.6, showscale = FALSE, name = "Spline Int")
# 
# # Add the spline_ind surface with a solid color
# p <- p %>% add_surface(x = ~seq_x, y = ~seq_y, z = ~spline_ind, 
#                        colorscale = list(c(0, 1), c("red", "red")), 
#                        opacity = 0.6, showscale = FALSE, name = "Spline Ind")
# 
# p
# 
# hist3D(x=midpoints_x,y=midpoints_y,z=(CLR(counts_impute)),add=TRUE)


