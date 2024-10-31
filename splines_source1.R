library(splines)
library(tidyr)
library(zCompositions)
library(ggplot2)
library(pracma)

truncate<-function(x,a,b){
  return(x[which(x>a&x<b)])
}


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

czm<-function(x){ #Impute 0-entries of histogram
  if(any(x<0)){
    stop("Can't contain negative ratios")
  }
  if(any(x==0)){
    n<-sum(x)
    x<-x/n
    x[which(x==0)]<- 2/3*1/n
  }
  return(x/sum(x))
}

zbSpline<-function(x,alfa=0.5,l=2,degree=3,bin_selection=doane,knots_inner=NULL,knots_full=NULL,res=1000){ #Kan udvides til at have støj med i form af weights
  #Kan også udvides til at x kan være en dataframe med histogram data
  #Skal skrives lidt om så histogrammernes mid-points of knots ikke nødvendigvis er ens
  #knots is without anchor knots
  #If x is a data.frame/histogram data, it must have column one of bin points and column y of bin values
  #Kan også udvides til at angive z-koefficienterne eller udregne de optimale.
  #knots_full kan angives hvis vi ikke vil have ens anchor knots
  options(warnPartialMatchArgs = TRUE)
  
  if(res<10){warning("Resolutions lower than 10 are likely to give numerical inaccuracies")}
  
  if(class(x)=="numeric"){
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
  } else if (class(x)=="data.frame"){
    n_bins<-nrow(x)
    bin_points<-x[,1]
    if(any(x[,2]==0)){
      print("Some bins had 0 values, imputing by czm. The imputed bins are:")  # Ensure the warning is displayed before the print
      print(x[which(x[,2]==0),1])
      bins_new<-czm(x[,2])
    }
    bins_new<-x[,2]
  } else stop("Error: x must be numeric or data frame")
  
  if(is.null(knots_full)){
    if (class(knots_inner) == "numeric") { 
      # If knots are supplied and are numeric
      g <- length(knots_inner) - 2
      knots <- c(rep(knots_inner[1], degree), knots_inner, rep(knots_inner[length(knots_inner)],degree))
    } else if(is.null(knots_inner)){ #Then we have no better option than to make knots in the histogram bins
      g<-n_bins
      knots<- c(rep(bin_points[1],degree),seq(bin_points[1],bin_points[n_bins],length.out=g+2),rep(bin_points[n_bins],degree))
    } else stop("Error: If knots are supplied, they must be numeric")
  } else if(class(knots_full)=="numeric"){
    knots<-knots_full
    g<-length(knots)-2-2*degree 
  }

  
  a<-min(knots); b<-max(knots)
  x_seq<-seq(a,b,length.out=res) #Need this for evaluating M, where higher res gives higher accuracy
  
  
  # Function factories
  Dj<-function(j){
    d<-numeric(g+degree+1-j)
    for(i in seq_along(d)){
      d[i]<-knots[i+degree+1]-knots[i+j]
    }
    return(diag(1/d)*(degree+1-j))
  }
  
  Lj <- function(j) {
    n <- g + degree + 1 - j
    L_empty <- matrix(0, nrow = n, ncol = n + 1)
    L_empty[cbind(1:n, 1:n)] <- -1
    L_empty[cbind(1:n, 2:(n + 1))] <- 1
    return(L_empty)
  }
  
  # Generate matrix Sl
  if(l==0){ Sl<-diag(1,g+degree+1)}
  Sl <- Dj(1) %*% Lj(1)
  if (l > 1) for (i in 2:l) Sl <- Dj(i) %*% Lj(i) %*% Sl

  
  # Generate matrix M
  # f <- function(i, j, a) splineDesign(knots, a, ord = degree + 1 - l, outer.ok = TRUE)[, i] * splineDesign(knots, a, ord = degree + 1 - l, outer.ok = TRUE)[, j]
  # M <- outer(1:(g + degree + 1 - l), 1:(g + degree + 1 - l), Vectorize(function(i, j) trapz(x_seq, f(i, j, x_seq))))
  
  M <- matrix(NA, nrow = (g + degree + 1-l), ncol = (g + degree + 1 - l))
  knots_M<-c(rep(knots_inner[1],degree-l),knots_inner,rep(knots_inner[length(knots_inner)],degree-l))
  f<-function(i,j,x){
    #mat<-splineDesign(knots[(degree-l+1):(length(knots)-l-1)],x,ord=degree+1-l,outer.ok=TRUE)
    mat<-splineDesign(knots_M,x,ord=degree+1-l,outer.ok=TRUE)
    #I use the lower 
    return(mat[,i]*mat[,j])
  }
  for (i in 1:nrow(M)) {
    for (j in 1:ncol(M)) {
      M[i,j]<- trapz(x_seq,f(i,j,x_seq))
      #M[i, j] <- ((b-a)/6)*(f(i,j,a)+4*f(i,j,(b-a)/2)+f(i,j,b)) #Approximate the integral by Simpson's rule
    }
  }
  
  
  #We create the Z-basis matrix
  B<-splineDesign(knots,bin_points,ord=(degree+1),outer.ok=TRUE,derivs=0)
  
  lambda<-numeric(g+degree+1)
  for(i in 1:length(lambda)){
    lambda[i]<-1/(knots[i+degree+1]-knots[i]) #Siden lambda[1]=lambda_-k er lambda[k+1]=lambda_1
  }
  D<-(degree+1)*diag(lambda)
  
  K<-matrix(0,nrow=g+degree+1,ncol=g+degree)
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
  B_x<-splineDesign(knots,x_seq,ord=(degree+1),outer.ok=TRUE)
  Z<-B_x%*%U
  
  G<-t(U)%*%((1-alfa)*t(Sl)%*%M%*%Sl+alfa*t(B)%*%B)%*%U
  g_mat<-alfa*t(U)%*%t(B)%*%CLR(bins_new)
  
  G_inv<-solve(G)
  
  z_coef<-G_inv%*%g_mat
  
  H<- alfa*B%*%U%*%G_inv%*%t(U)%*%t(B)
  
  structure(list("Z_basis"=Z,"z_coef"=z_coef,"Z_spline"=Z%*%z_coef,
                 "U"=U,"H"=H,"M"=M,"Sl"=Sl,"K"=K,"knots"=knots,"x"=x,"x_seq"=x_seq,
                 "bin_points"=bin_points,"bin_values"=bins_new,"degree"=degree,"g"=g,"alfa"=alfa),class="zbSpline")
}

plot.zbSpline<-function(Z,what="basis",include_knots=TRUE,include_hist=TRUE){
  if(what=="Z-basis"){
    Z_df<-as.data.frame(Z$Z_basis)
    Z_df[Z_df == 0] <- NA
  }
  else if(what=="O-basis"){
    Z_df<-as.data.frame(Z$O)
    Z_df[Z_df == 0] <- NA
  }
  else if(what=="C-basis"){
    Z_df<-as.data.frame(Z$C_basis)
  }
  else if(what=="Z-spline"){
    Z_df<-as.data.frame(Z$Z_spline)
  }
  else if(what=="O-spline"){
    Z_df<-as.data.frame(Z$O%*%Z$z_coef)
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
  if(include_hist==TRUE & what%in% c("spline","O-spline")){ #Dont make histogram data when irrelevant
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


cross_validate <- function(x, alfa_grid = seq(0.01, 1, length.out = 30), k = 3, l = 2, knots_inner = NULL) {
  cv <- numeric(length = length(alfa_grid))
  for(i in seq_along(cv)) {
    Z <- zbSpline(x, alfa = alfa_grid[i], degree = k, l = l, knots_inner = knots_inner)
    y <- CLR(Z$bin_values)
    n <- length(y)
    B <- splineDesign(Z$knots, Z$bin_points, ord = Z$degree + 1, outer.ok = TRUE) # Using U and spline design to evaluate the basis
    spline <- B %*% Z$U%*%Z$z_coef
    cv[i] <- mean((y - spline)^2) / ((1 - sum(diag(Z$H)) / n)^2)
  }
  structure(list("alfa" = alfa_grid, "cv" = cv, "optimum" = alfa_grid[which.min(cv)]), class = "zbCV")
}
#For det mere generelle skal vi sikre at k og l matcher med Z. Det kan gøres ved at arve Z i stedet for x

plot.zbCV<-function(cv){
  cv_df<-data.frame("x"=cv$alfa,"GCV"=cv$cv)
  ggplot(cv_df)+geom_point(mapping=aes(x=x,y=GCV),color="pink")+geom_line(mapping=aes(x=x,y=GCV))+
   xlab(expression("smoothing parameter " * alpha))+geom_vline(xintercept=cv$optim,color="gray",linetype="dashed")+
   theme_minimal()
}

orthogonalise<-function(Z,dec="Cholesky"){
  M <- matrix(NA, nrow = (Z$g + Z$degree + 1), ncol = (Z$g + Z$degree + 1))
  f<-function(i,j,x){
    mat<-splineDesign(Z$knots,x,ord=Z$degree+1,outer.ok=TRUE)
    return(mat[,i]*mat[,j])
  }
  for (i in 1:nrow(M)) {
    for (j in 1:ncol(M)) {
      M[i,j]<- trapz(Z$x_seq,f(i,j,Z$x_seq))
    }
  }
  Sigma<-t(Z$U)%*%M%*%Z$U
  Sigma_inv<-solve(Sigma)
  if(dec=="Cholesky"){
    Phi<-chol(Sigma_inv)
  }
  if(dec=="svd"){
    decomp<-svd(Sigma_inv)
    Phi<-decomp$v%*%diag(sqrt(decomp$d))
  }
  O<-Z$Z_basis%*%t(Phi) #This is the reverse way of what is indicated by the paper
  Z_list<-as.list(Z)
  Z_list$O<-O
  Z_list$Phi<-Phi
  structure(Z_list,class="zbSpline")
} #I reuse the optimal coefficients from another basis, but can't really figure out
#how to optimise when the basis is orthogonal.

cbSpline<-function(Z,what="original"){
 if(what=="original"){
   zeta<-exp(Z$Z_basis) #We transform the basis
   for(i in 1:ncol(zeta)){
     zeta[,i]<-zeta[,i]/trapz(Z$x_seq,zeta[,i]) #Normalise
   }
 } else if(what=="orthogonal"){
   zeta<-exp(Z$O)
   for(i in 1:ncol(zeta)){
     zeta[,i]<-zeta[,i]/trapz(Z$x_seq,zeta[,i])
   }
 } else(stop("Error: ZB-splines must be either original or orthogonal"))
  z<-Z$z_coef #Extract the coefficients
  terms=matrix(NA,nrow=nrow(zeta),ncol=ncol(zeta))
  for(i in 1:ncol(zeta)){ #Create all of the terms for the spline composition
    terms[,i]<-zeta[,i]^z[i] #Perform powering on the terms
  }
  xi<-terms[,1]
  for (i in 2:ncol(zeta)){ #And combine them with perturbation
    xi<-xi*terms[,i]
  }
  xi<-xi/trapz(Z$x_seq,xi)
  Z_list<-as.list(Z)
  Z_list$C_basis<-zeta
  Z_list$C_spline<-xi
  structure(Z_list,class="zbSpline")
}
