library(splines)
library(tidyr)
library(zCompositions)
library(ggplot2)
library(pracma)
library(gridExtra)
library(plot3D)
library(plotly)
library(viridis)

#Is used for the scott rule
skewness<-function(x){
  return(mean((x-mean(x))^3/sd(x)^3))
}

doane<-function(data){
  n <- length(data)
  g1 <- skewness(data)
  sigma_g1 <- sqrt(6 * (n - 2) / ((n + 1) * (n + 3)))
  k <- 1 + log2(n) + log2(1 + abs(g1) / sigma_g1)
  return(ceiling(k))
}
sturge<-function(data){
  ceiling(1+log2(length(data)))
}
scott <- function(data) {
  n <- length(data)
  sigma <- sd(data)
  h <- 3.5*sigma * n^(-1/3)
  k<-(max(data)-min(data))/h
  return(ceiling(k))
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
  int_over_x <- apply(f, 2, function(fy) pracma::trapz(x, fy)) #Integrate over x first
  # Integrate the result over y, by Fubini
  integral <- pracma::trapz(y, int_over_x)
  
  return(integral)
}

#The CLR will work in any dimension due to vectorisation
CLR<-function(x){
  return(log(x)-mean(log(x)))
}

clr <- function(f,x) { #x is used to get the range over which f is evaluated and the resolution
  a<-min(x)
  b<-max(x)
  n<-length(x)
  t_values <- seq(a, b, length.out = n)
  f_values <- f(t_values)
  
  if (any(f_values <= 0)) {
    stop("The function must be positive over the entire interval [a, b].")
  }
  # Calculate the integral of log(f) over [a, b] using the trapezoidal rule
  log_f_values <- log(f_values)
  integral_log_f <- trapz(t_values,log_f_values)
  # Return the clr-transformed function
  clr_function <- function(t) {
    log_f_t <- log(f(t))
    clr_value <- log_f_t - (integral_log_f / (b - a))
    return(clr_value)
  }
  return(clr_function)
}


czm<-function(x){ #Impute 0-entries of 1-D histogram x
  if(any(x<0)){
    stop("Can't contain negative ratios")
  }
  if(any(x==0)){
    n<-sum(x)
    x<-x/n
    x[which(x==0)]<- 2/3*1/n
  }
  return(x/sum(x))
} #Impute 1-D histogram x


#We interpret neightbours as the directly adjacent cells because it is much faster computationally
get_neighbors <- function(mat, i, j) {
  offsets <- c(-1, 0, 1)
  neighbors <- numeric(0)
  n_row <- nrow(mat)
  n_col <- ncol(mat)
  
  for (dx in offsets) {
    for (dy in offsets) {
      if (dx == 0 && dy == 0) next  # skip the cell itself
      x <- i + dx
      y <- j + dy
      if (x >= 1 && x <= n_row && y >= 1 && y <= n_col) {
        neighbors <- c(neighbors, mat[x, y])
      }
    }
  }
  neighbors
}

impute_zeros <- function(counts) {
  counts <- as.matrix(counts)
  repeat {
    zero_positions <- which(counts == 0, arr.ind = TRUE)
    if (nrow(zero_positions) == 0) break
    
    changed <- FALSE
    for (k in seq_len(nrow(zero_positions))) {
      i <- zero_positions[k, 1]
      j <- zero_positions[k, 2]
      # Double-check it's still zero (may have been imputed this iteration)
      if (counts[i, j] == 0) {
        neighbor_values <- get_neighbors(counts, i, j)
        non_zero_neighbors <- neighbor_values[neighbor_values > 0]
        if (length(non_zero_neighbors) > 0) {
          gm <- exp(mean(log(non_zero_neighbors))) * 2/3
          counts[i, j] <- gm
          changed <- TRUE
        }
      }
    }
    # If no zeros changed, we're done
    if (!changed) break
  }
  counts
}



zbSpline1D<-function(x,alfa=0.5,l=2,k=3,bin_selection=doane,knots_inner=NULL,res=1000){ #Kan udvides til at have støj med i form af weights
  #x is a vector of observations or a dataframe. If it is a dataframe, the bin_points are in the first column and bin_values are in the second
  #alfa is smoothing parameter, and may be found by cross-validation (other function)
  #l is the derivative in smoothing term
  #k is the degree of the polynomial
  #bin_selection is either numeric or function for choosing bin_points
  #knots_inner provides either the inner knots, or they will be chosen in the bin_points
  options(warnPartialMatchArgs = TRUE)
  
  if(res<10){warning("Resolutions lower than 10 are likely to give numerical inaccuracies")}
  if (alfa<=0 | alfa>1) stop ("parameter alpha must be from interval (0,1]")
  
  if(class(x)=="numeric"){ #Creates the histogram using the bin_selection described
    if(class(bin_selection)=="function"){
      n_bins<-bin_selection(x)
    } else if(class(bin_selection)=="numeric"){
      n_bins<-bin_selection
    } else{stop("Error: bin_selection must be either a function or numeric")}
    breaks <- seq(min(x), max(x), length.out = n_bins+1)
    hist_info <- hist(x, breaks=breaks, plot = FALSE, include.lowest = TRUE, right = FALSE)
    bin_points<-hist_info$mids
    
    bins<-as.numeric(t(as.matrix(hist_info$counts)))
    bins_new<-czm(bins)
  } else if (class(x)=="data.frame"){ #Otherwise takes the histogram from the 
    n_bins<-nrow(x)
    bin_points<-x[,1]
    if(any(x[,2]==0)){ #Impute if necessary
      print("Some bins had 0 values, imputing by czm. The imputed bins are:")  # Ensure the warning is displayed before the print
      print(x[which(x[,2]==0),1])
      bins_new<-czm(x[,2])
    } else{bins_new<-x[,2]}
  } else stop("Error: x must be numeric or data frame")
  
  
  if(is.null(knots_inner)){
    knots_inner<-bin_points
  }
  g <- length(knots_inner) - 2
  knots <- c(rep(knots_inner[1], k), knots_inner, rep(knots_inner[length(knots_inner)],k))
  
  a<-min(knots); b<-max(knots)
  x_seq<-seq(a,b,length.out=res) #Need this for evaluating M, where higher res gives higher accuracy
  
  # Function factories
  Dj<-function(j){
    d<-numeric(g+k+1-j)
    for(i in seq_along(d)){
      d[i]<-knots[i+k+1]-knots[i+j]
    }
    return(diag(1/d)*(k+1-j))
  }
  
  Lj <- function(j) {
    n <- g + k + 1 - j
    L_empty <- matrix(0, nrow = n, ncol = n + 1)
    L_empty[cbind(1:n, 1:n)] <- -1
    L_empty[cbind(1:n, 2:(n + 1))] <- 1
    return(L_empty)
  }
  
  # Generate matrix Sl
  if(l==0){ Sl<-diag(1,g+k+1)}
  Sl <- Dj(1) %*% Lj(1)
  if (l > 1) for (i in 2:l) Sl <- Dj(i) %*% Lj(i) %*% Sl
  
  
  # Generate matrix M
  M <- matrix(NA, nrow = (g + k + 1-l), ncol = (g + k + 1 - l))
  knots_M<-c(rep(knots_inner[1],k-l),knots_inner,rep(knots_inner[length(knots_inner)],k-l))
  f<-function(i,j,x){
    mat<-splineDesign(knots_M,x,ord=k+1-l,outer.ok=TRUE)
    return(mat[,i]*mat[,j])
  }
  for (i in 1:nrow(M)) {
    for (j in 1:ncol(M)) {
      M[i,j]<- trapz(x_seq,f(i,j,x_seq))
    }
  }
  
  #We create the Z-basis matrix
  B<-splineDesign(knots,bin_points,ord=(k+1),outer.ok=TRUE,derivs=0)
  
  lambda<-numeric(g+k+1)
  for(i in 1:length(lambda)){
    lambda[i]<-1/(knots[i+k+1]-knots[i]) #Siden lambda[1]=lambda_-k er lambda[k+1]=lambda_1
  }
  D<-(k+1)*diag(lambda)
  
  K<-matrix(0,nrow=g+k+1,ncol=g+k)
  K[1,1]<-1
  for (i in 2:nrow(K)) {
    if (i - 1 <= ncol(K)) {
      K[i, i - 1] <- -1
    }
    if (i <= ncol(K)) {
      K[i, i] <- 1
    }
  }
  
  U<-D%*%K
  
  #For evaluating across the whole x_axis
  B_x<-splineDesign(knots,x_seq,ord=(k+1),outer.ok=TRUE)
  
  #We make the basis functions evaluated across the axis and in the histograms
  Z<-B_x%*%U
  Z_hist<-B%*%U
  
  #Find optimal splines
  G<-t(U)%*%((1-alfa)*t(Sl)%*%M%*%Sl+alfa*t(B)%*%B)%*%U
  g_mat<-alfa*t(U)%*%t(B)%*%CLR(bins_new)
  G_inv<-solve(G)
  z_coef<-G_inv%*%g_mat
  
  #Used for cross-validation
  H<- alfa*B%*%U%*%G_inv%*%t(U)%*%t(B)
  
  #And we convert all of it back to probability scale
  
  
  zeta<-exp(Z) #We transform the basis
  #Normalising basis splines may give numerical instability
  #Having the basis on a small domain is also problematic..
  #Perhaps the clr-space is not properly scaled, but it does not seem so
  #Yet this zeta matrix can become numerically unstable
 
  
  terms=matrix(NA,nrow=nrow(zeta),ncol=ncol(zeta))
  for(i in 1:ncol(zeta)){ #Create all of the terms for the spline composition
    terms[,i]<-zeta[,i]^z_coef[i] #Perform powering on the terms
  }
  xi<-terms[,1]
  for (i in 2:ncol(zeta)){ #And combine them with perturbation
    xi<-xi*terms[,i]
  }
  xi<-xi/trapz(x_seq,xi)
  
  structure(list("Z_basis"=Z,"z_coef"=z_coef,"Z_spline"=Z%*%z_coef,"Z_spline_hist"=Z_hist%*%z_coef,
                 "C_basis"=zeta,"C_spline"=xi,
                 "U"=U,"H"=H,"M"=M,"Sl"=Sl,"K"=K,"knots"=knots,"x"=x,"x_seq"=x_seq,
                 "bin_points"=bin_points,"bin_values"=bins_new,"degree"=k,"g"=g,"alfa"=alfa),class="zbSpline1D")
}
zbSpline2D<-function(x,y=NULL,knots_x_inner,knots_y_inner,alfa=0.5,rho=NULL,bin_selection=doane,k=3,l=3,u=2,v=2,res=200,range_x=NULL,range_y=NULL,calc_dens=TRUE){
  #x,y are either vectors of observation (x,y) or x is a list of length 3 with the first two being the midpoints and the third being a matrix of bin_values
  #alfa is the smoothing parameter used in Machalova21 and Hron23, and is converted to the smoothing parameter of skorna unless this is specified directly
  #bin_selection is either a function or numeric. In the first case it will calculate the number of bins by the specified rule-of-thumb. In the second it will use the provided number of bins.
  #bin_selection is not used in x and y are numerics.
  #knots_x_inner and knots_y_inner are the inner knots and boundary knots where the basis splines should be created. Does not include anchor knots, since these are calculated in the function when necessary
  #k,l are the degrees of the basis splines
  #u,v are the orders of derivatives
  #res is what resolution the densities are evaluated in
  #calc_dens is true or false. If false, can save some computational time
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
    
    if(class(bin_selection)=="function"){
      m<-bin_selection(x)
      n<-bin_selection(y)
    }
    if(class(bin_selection)=="numeric"){
      m<-bin_selection[1]
      n<-bin_selection[2]
    }
    
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
    
    F_mat<-CLR(hist_data)
  }else {stop("Error: Data must be either histogram data x or numerics (x,y).")}
  
  #We make the knot sequence
  #We start with assuming that the inner knots are provided
  #Actually there is too much numeric instability if they are not suitable provided up front
  
  g<-length(knots_x_inner)-2
  h<-length(knots_y_inner)-2
  knots_x<-c(rep(knots_x_inner[1],k),knots_x_inner,rep(knots_x_inner[length(knots_x_inner)],k))
  knots_y<-c(rep(knots_y_inner[1],l),knots_y_inner,rep(knots_y_inner[length(knots_y_inner)],l))
  
  

  
  if(is.null(range_x)){
    seq_x<-seq(min(knots_x),max(knots_x),length.out=res)
  } else{seq_x<-seq(range_x[1],range_x[2],length.out=res)
    }
  
  if(is.null(range_y)){
    seq_y<-seq(min(knots_y),max(knots_y),length.out=res)
    } else{seq_y<-seq(range_y[1],range_y[2],length.out=res)
  }
  
  a<-min(seq_x); b<-max(seq_x); c<-min(seq_y); d<-max(seq_y)
  
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
  if(v==0){ Sv<-diag(1,h+l+1)} #This seems to not be the way to go
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
  z_opt<-matrix(coefficients[1:((g+k)*(h+l))],nrow=k+g,ncol=l+h) #The weird indexing is to extract the correct information for creating R_opt
  v_opt<-coefficients[((g+k)*(h+l)+1):((g+k)*(h+l)+k+g)]
  u_opt<-coefficients[((g+k)*(h+l)+k+g+1):((g+k)*(h+l)+k+g+h+l)]
  #This is unnecessary for all but creating
  R_opt<-rbind(z_opt,t(u_opt))
  R_opt<-cbind(R_opt,c(v_opt,0))
  spline_hist<-t(Z_x_bar)%*%R_opt%*%Z_y_bar
  
  spline_int<- Z_x_seq%*%z_opt%*%t(Z_y_seq)
  spline_int<-spline_int-trapz2d(seq_x,seq_y,spline_int)/((b-a)*(d-c))
  spline_x<-Z_x_seq%*%v_opt
  spline_x<-spline_x-trapz(seq_x,spline_x)/(b-a)
  spline_y<-Z_y_seq%*%u_opt
  spline_y<-spline_y-trapz(seq_y,spline_y)/(d-c)
  spline_ind<-outer(spline_x,spline_y,FUN="+")
  spline_ind<-matrix(spline_ind,nrow=res,ncol=res)
  spline_opt<-spline_int+spline_ind
  
  #And we evaluate the residuals
  spline_hist_clr<-t(Z_x_bar)%*%R_opt%*%Z_y_bar #The spline evaluated in the histogram bins
  res_clr<-F_mat-spline_hist_clr #The residuals
  spline_hist_bhs<-exp(spline_hist_clr)
  spline_hist_bhs<-spline_hist_bhs/trapz2d(midpoints_x,midpoints_y,spline_hist_bhs)
  res_bhs<-hist_data-spline_hist_bhs
  
  
  #Calculate relative simplical deviance
  simp_d<-trapz2d(seq_x,seq_y,spline_int^2)
  rel_simp_d<-simp_d/trapz2d(seq_x,seq_y,spline_opt^2)
  
  
  #And we map it back to density
  #This should be done via the B^2 basis functions, but I won't focus too much on it
  if(calc_dens){
    C_spline<-exp(spline_opt)
    C_spline<-C_spline/trapz2d(seq_x,seq_y,C_spline)
    C_spline_ind<-exp(spline_ind)
    C_spline_ind<-C_spline_ind/trapz2d(seq_x,seq_y,C_spline_ind)
    C_spline_int<-exp(spline_int)
    C_spline_int<-C_spline_int/trapz2d(seq_x,seq_y,C_spline_int)
    C_spline_x<-exp(spline_x)
    C_spline_x<-C_spline_x/trapz(seq_x,C_spline_x)
    C_spline_y<-exp(spline_y)
    C_spline_y<-C_spline_y/trapz(seq_y,C_spline_y)
  }

  
  if (calc_dens) {
    C_spline <- exp(spline_opt)
    C_spline <- C_spline / trapz2d(seq_x, seq_y, C_spline)
    
    C_spline_ind <- exp(spline_ind)
    C_spline_ind <- C_spline_ind / trapz2d(seq_x, seq_y, C_spline_ind)
    
    C_spline_int <- exp(spline_int)
    C_spline_int <- C_spline_int / trapz2d(seq_x, seq_y, C_spline_int)
    
    C_spline_x <- exp(spline_x)
    C_spline_x <- C_spline_x / trapz(seq_x, C_spline_x)
    
    C_spline_y <- exp(spline_y)
    C_spline_y <- C_spline_y / trapz(seq_y, C_spline_y)
  }
  
  res_list <- list(Z_spline=spline_opt,Z_spline_int=spline_int,Z_spline_ind=spline_ind,
                   Z_spline_x=spline_x,Z_spline_y=spline_y,
                   seq_x=seq_x,seq_y=seq_y,midpoints_x=midpoints_x,midpoints_y=midpoints_y,
                   F_mat=F_mat,hist_data=hist_data,residuals_clr=res_clr,residuals_bhs=res_bhs,
                   R_opt=R_opt,Z_opt=z_opt,v_opt=v_opt,u_opt=u_opt,Z_x_bar=Z_x_bar,Z_y_bar=Z_y_bar,
                   bbZ_bar=bbZ_bar,G=G,spline_hist=spline_hist,rho=rho,alfa=alfa,sd=simp_d,rsd=rel_simp_d)
  
  if(calc_dens){
    res_list <- c(res_list,list(C_spline=C_spline,C_spline_ind=C_spline_ind,C_spline_int=C_spline_int,
                                C_spline_x=C_spline_x,C_spline_y=C_spline_y,
                                spline_max=which(C_spline==max(C_spline),arr.ind=TRUE)))
  }
  
  structure(res_list,class="zbSpline2D")
}

summary.zbSpline2D <- function(Z) {
  cat("Spline evaluated over the intervals x in (", min(Z$seq_x), ",", max(Z$seq_x), 
      ") and y in (", min(Z$seq_y), ",", max(Z$seq_y), ")\n")
  cat("RSD:", Z$rsd, "\n")
  cat("Arg max:", Z$seq_x[Z$spline_max[1]],Z$seq_y[Z$spline_max[2]], "\n")
}

cross_validate1D <- function(x, alfa_seq = seq(0.01, 1, length.out = 30), k = 3, l = 2, knots_inner = NULL) {
  cv <- numeric(length = length(alfa_seq))
  for(i in seq_along(cv)) {
    Z <- zbSpline1D(x, alfa = alfa_seq[i], k = k, l = l, knots_inner = knots_inner)
    y <- CLR(Z$bin_values)
    n <- length(y)
    B <- splineDesign(Z$knots, Z$bin_points, ord = Z$degree + 1, outer.ok = TRUE) # Using U and spline design to evaluate the basis
    spline <- B %*% Z$U%*%Z$z_coef
    cv[i] <- mean((y - spline)^2) / ((1 - sum(diag(Z$H)) / n)^2)
  }
  structure(list("alfa" = alfa_seq, "cv" = cv, "optimum" = alfa_seq[which.min(cv)]), class = "zbCV")
}
cross_validate2D<-function(x,y=NULL,knots_x_inner,knots_y_inner,alfa_seq=seq(0.01,0.99,length.out=20),res=100,k=2,l=2,u=1,v=1,bin_selection=doane,range_x=NULL,range_y=NULL){
  cv<-numeric(length(alfa_seq)) #Will capture the cross-validations scores
  
  #We make this code just to have the dimensions to make the permutation matrix from.
  #It is very inefficient, but it works so long as no coefficients are repeated
  Z_temp<-try(zbSpline2D(x,y,alfa=0.5,bin_selection=bin_selection,
                         knots_x_inner=knots_x_inner,knots_y_inner=knots_y_inner,
                         k=k,l=l,u=u,v=v,res=200,range_x=range_x,range_y=range_y,calc_dens=FALSE),
              silent=TRUE)
  if (inherits(Z_temp, "try-error") && grepl("Lapack routine", Z_temp)) {
    stop("Error with optimisation. Ensure knots are inside the histogram range.")
  }
  
  dat<-list(Z_temp$midpoints_x,Z_temp$midpoints_y,Z_temp$hist_data)
  
  csR<-as.vector(Z_temp$R_opt)
  csZ<-c(as.vector(Z_temp$Z_opt),Z_temp$v_opt,Z_temp$u_opt)
  P<-matrix(0,nrow=length(csR),ncol=length(csZ))#We make the permutation matrix
  for(i in 1:nrow(P)){
    sigma_i<-which(csZ==csR[i])
    P[i,sigma_i]<-1
  }
  #Compute the cross-validation scores at each alfa
  #Can reuse the same histogram as earlier
  for(i in seq_along(alfa_seq)){
    Z<-zbSpline2D(x=dat,alfa=alfa_seq[i],knots_x_inner=knots_x_inner,
                  knots_y_inner=knots_y_inner,bin_selection=bin_selection,
                  k=k,l=l,u=u,v=v,res=res,range_x=range_x,range_y,range_y,calc_dens=FALSE)
    H<-kronecker(t(Z$Z_y_bar),t(Z$Z_x_bar))%*%P%*%solve(Z$G)%*%Z$bbZ_bar #Definition of H
    m<-length(Z$midpoints_x)
    n<-length(Z$midpoints_y)
    cv[i]<-mean((Z$F_mat-Z$spline_hist)^2)/((1-sum(diag(H))/(m*n))^2)
  }
  structure(list("alfa" = alfa_seq, "cv" = cv, "optimum" = alfa_seq[which.min(cv)]),class="zbCV")
}




plot.zbSpline1D<-function(Z,what="Z-basis",include_knots=TRUE,include_hist=TRUE){
  if(what=="Z-basis"){
    Z_df<-as.data.frame(Z$Z_basis)
    Z_df[Z_df == 0] <- NA #To make sure it behaves as undefined outside the domain
  }
  else if(what=="C-basis"){
    Z_df<-as.data.frame(Z$C_basis)
  }
  else if(what=="Z-spline"){
    Z_df<-as.data.frame(Z$Z_spline)
  }
  else if(what=="C-spline"){
    Z_df<-as.data.frame(Z$C_spline)
  }
  x_seq<-Z$x_seq
  Z_df$x<-x_seq
  Z_long<-pivot_longer(Z_df,cols = -x, names_to = "basis", values_to = "value")
  rep_points<-data.frame("x"=Z$bin_points,"y"=Z$bin_values/(trapz(Z$bin_points,Z$bin_values))) #y is divided with the integral to have it at probability scale
  #alternatively we could use bin_width when normalising
  p<-ggplot() +
    geom_line(Z_long,mapping=aes(x=x, y=value, color=basis)) +
    theme_minimal() +
    theme(legend.position = "none")
  if(include_knots==TRUE){
    p<-p+geom_vline(xintercept=Z$knots,color="gray",linetype="dashed")
  }
  if(include_hist==TRUE & what=="Z-spline"){ #Dont make histogram data when irrelevant
    p<-p+geom_point(data=rep_points,mapping=aes(x=x,y=CLR(y)),shape=8,color="blue")
  }
  if(include_hist==TRUE & what=="C-spline"){
    p<-p+geom_point(data=rep_points,mapping=aes(x=x,y=y),shape=8,color="blue")
  }
  if(what=="C-basis"){
    p<-p+ylim(0,max(Z_long$value))
  }
  return(p)
}
plot.zbSpline2D<-function(Z,type="static",what="full",scale="clr",title="",plot_hist=FALSE,plot=TRUE,xlab="x",ylab="y",xlim=NA,ylim=NA,theta=325,phi=30,point_col="black"){
  #type can be static or interactive
  #what can be either full, independent, interaction, geom_X or geom_Y
  #plot_hist is whether we plot the underlying histogram
  #scale chooses between the clr transform or the density
  if (!type %in% c("static", "interactive")) {
    stop("Error: 'type' must be either 'static' or 'interactive'")
  }
  
  valid_what <- c("full", "independent", "interaction", "geom_X", "geom_Y")
  if (!what %in% valid_what) {
    stop("Error: 'What' must be one of 'full', 'independent', 'interaction', 'geom_X', or 'geom_Y'")
  }
  if (!is.logical(plot_hist)) {
    stop("Error: 'plot_hist' must be a logical value (TRUE or FALSE)")
  }
  if (!scale %in% c("clr", "density")) {
    stop("Error: 'scale' must be either 'clr' or 'density'")
  }
  
  #We can update all of these to finding a response and range based on the if statements
  if (what == "full") {
    outcome <- if (scale == "clr") Z$Z_spline else Z$C_spline
  } else if (what == "interaction") {
    outcome <- if (scale == "clr") Z$Z_spline_int else Z$C_spline_int
  } else if (what == "independent") {
    outcome <- if (scale == "clr") Z$Z_spline_ind else Z$C_spline_ind
  }
  
  if(is.na(xlim[1])){
    xlim<-c(min(Z$seq_x),max(Z$seq_x))
  }
  if(is.na(ylim[1])){
    ylim<-c(min(Z$seq_y),max(Z$seq_y))
  }
  
  
  if (plot_hist) {
    hist_response <- if (scale == "clr") {
      Z$F_mat
    } else {
      Z$hist_data / trapz2d(Z$midpoints_x, Z$midpoints_y, Z$hist_data)
    }
    grid <- expand.grid("x" = Z$midpoints_x, "y" = Z$midpoints_y)
  }
  
  # Plotting based on 'what' and 'type'
  if (what %in% c("full", "independent", "interaction")) {
    if (type == "static") {
      # Static plotting using persp3D
      persp3D(
        x = Z$seq_x, 
        y = Z$seq_y, 
        z = outcome, 
        col = viridis(50),
        theta = theta,          
        phi = phi,            
        ticktype = "detailed", 
        nticks = 3,
        xlab = xlab,          
        ylab = ylab,
        zlab = "",          
        bty = "b2",
        colkey = FALSE,
        main = title,
        cex.axis = 0.5,
        xlim=xlim,
        ylin=ylim
      )
      if (plot_hist) {
        points3D(
          x = grid$x,
          y = grid$y,
          z = as.vector(hist_response),
          add = TRUE,
          col = point_col,
          pch = 19,
          cex = 0.5
        )
      }
    } else {
      # Interactive plotting using plotly
      p <- plot_ly() %>%
        add_surface(x = Z$seq_x, y = Z$seq_y, z = t(outcome))%>%layout(title=title,
                                                                       scene=list(xaxis=list(title=xlab),
                                                                                  yaxis=list(title=ylab)))
      if (plot_hist) {
        p <- p %>%
          add_markers(
            x = grid$x,
            y = grid$y,
            z = as.vector(hist_response),
            marker = list(size = 3, color = "black")
          )
      }
      p
    }
  } else {
    # 2D plotting for 'geom_X' and 'geom_Y'
    axis_var <- if (what == "geom_X") "x" else "y"
    seq_var <- Z[[paste0("seq_", axis_var)]]
    if(scale=="clr"){
      spline_var <- Z[[paste0("Z_spline_", axis_var)]]
    } else if (scale=="density"){
      spline_var <- Z[[paste0("C_spline_", axis_var)]]
    }
    
    plot_data <- data.frame("x" = seq_var, "z" = spline_var)
    p <- ggplot(plot_data) +
      geom_line(aes(x = x, y = z)) +
      xlab(what)
    if (plot_hist) {
      midpoints_var <- Z[[paste0("midpoints_", axis_var)]]
      if(scale=="clr"){
        bins <- if (axis_var == "x") rowMeans(Z$F_mat) else colMeans(Z$F_mat)
      } else if(scale=="density") {
        bins <- if (axis_var == "x") rowMeans(Z$hist_data) else colMeans(Z$hist_data)
        bins<-bins/trapz(midpoints_var,bins)
      }
      plot_hist_data <- data.frame("x" = midpoints_var, "z" = bins)
      p <- p +
        geom_point(aes(x = x, y = z), data = plot_hist_data, shape = 8, color = "blue")
    }
    p
  }
}


plot.zbCV<-function(cv){
  cv_df<-data.frame("x"=cv$alfa,"GCV"=cv$cv)
  ggplot(cv_df)+geom_point(mapping=aes(x=x,y=GCV),color="pink")+geom_line(mapping=aes(x=x,y=GCV))+
    xlab(expression("smoothing parameter " * alpha))+geom_vline(xintercept=cv$optim,color="gray",linetype="dashed")+
    theme_minimal()
}


perm_test<-function(x,y,kx,ky,alfa_seq=seq(0.01,1,length.out=10),bin_selection=doane,k=2,l=2,u=1,v=1,res=100,K=500){
  cv_true<-cross_validate2D(x,y,knots_x_inner=kx,knots_y_inner=ky,alfa_seq=alfa_seq,k=k,l=l,u=u,v=v)
  alfa_true<-cv_true$optimum
  true_spline<-zbSpline2D(x,y,alfa=alfa_true,bin_selection=bin_selection,
                       knots_x_inner=kx,knots_y_inner=ky,
                       k=k,l=l,u=u,v=v,res=res,calc_dens=FALSE)
  rsd_true<-true_spline$rsd
  
  x_new <- sample(x, replace = FALSE)
  y_new <- sample(y, replace = FALSE)
  kx_new <- seq(min(x_new), max(x_new), length.out = length(kx))
  ky_new <- seq(min(y_new), max(y_new), length.out = length(ky))
  
  cv_perm<-cross_validate2D(x_new,y_new,knots_x_inner=kx_new,
                            knots_y_inner=ky_new,alfa_seq=alfa_seq,
                            k=k,l=l,u=u,v=v)
  alfa_perm<-cv_perm$optimum
  
  rsd_vec<-numeric(K)
  for (i in 1:K) {
    repeat {  # Use repeat to retry in case of an error
      tryCatch({
        # permute data and make knots
        x_new <- sample(x, replace = FALSE)
        y_new <- sample(y, replace = FALSE)
        kx_new <- seq(min(x_new), max(x_new), length.out = 4)
        ky_new <- seq(min(y_new), max(y_new), length.out = 4)
        
        # Run the bivariate function and store result
        sim_rsd <- zbSpline2D(x_new, y_new, alfa = alfa_perm, bin_selection = bin_selection, 
                             knots_x_inner = kx_new, knots_y_inner = ky_new, 
                             k = k, l = l, u = u, v = v, res = res,calc_dens=FALSE)
        rsd_vec[i] <- sim_rsd$rsd
        
        break  # Exit repeat if no error occurs
        
      }, error = function(e) {
        # Check if it's the specific error, if so, retry
        if (grepl("system is computationally singular", e$message)) {
          message("Singular matrix encountered, retrying iteration ", i)
        } else {
          stop(e)  # If it's another error, stop execution
        }
      })
    }
  }
  p<-length(which(rsd_vec>=rsd_true))/K
  return(list("p_val"=p,"alfa_true"=alfa_true,"alfa_perm"=alfa_perm,"rsd_true"=rsd_true,"rsd_perms"=rsd_vec,"hist"=true_spline$hist_data))
}
