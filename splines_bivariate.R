setwd("~/Universitet-mat/.Speciale/R")
source("splines_source1.R")
library(gridExtra)
library(splines)
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
    counts <- counts_new
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

clr2d<-function(f,x,y){
  f_vals<-log(outer(x,y,FUN=f))
  int<-trapz2d(x,y,f_vals)
  area <- (max(x) - min(x)) * (max(y) - min(y))
  clr2_function<-function(x,y){
    log(f(x,y))-int/area
  }
  return(clr2_function)
}

#Denne funktion mangler mulighed for selv at finde knots og mulighed for at specificere antal bins
bivariate<-function(x,y=NULL,alfa=0.5,rho=NULL,bin_selection=doane,knots_x_inner,knots_y_inner,k=3,l=3,u=2,v=2,res=200){
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
  }
  else {stop("Error: Data must be either histogram data x or numerics (x,y).")}
  
  #We make the knot sequence
  #We start with assuming that the inner knots are provided
  #Actually there is too much numeric instability if they are not suitable provided up front
  
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
  Z_opt<-matrix(coefficients[1:((g+k)*(h+l))],nrow=k+g,ncol=l+h) #The weird indexing is to extract the correct information for creating R_opt
  v_opt<-coefficients[((g+k)*(h+l)+1):((g+k)*(h+l)+k+g)]
  u_opt<-coefficients[((g+k)*(h+l)+k+g+1):((g+k)*(h+l)+k+g+h+l)]
  #This is unnecessary for all but creating
  R_opt<-rbind(Z_opt,t(u_opt))
  R_opt<-cbind(R_opt,c(v_opt,0))
  spline_hist<-t(Z_x_bar)%*%R_opt%*%Z_y_bar
  
  
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
  #This should be done via the B^2 basis functions, but I won't focus too much on it
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
  
  #spline_x is the geometric marginal
  structure(list("Z_spline"=spline_opt, "Z_spline_int"=spline_int, "Z_spline_ind"=spline_ind, "Z_spline_x"=spline_x,
                 "Z_spline_y"=spline_y,
                 "C_spline"=C_spline,
                 "C_spline_ind"=C_spline_ind,
                 "C_spline_int"=C_spline_int,
                 "C_spline_x"=C_spline_x,
                 "C_spline_y"=C_spline_y,
                 "seq_x"=seq_x,"seq_y"=seq_y,
                 "midpoints_x"=midpoints_x,"midpoints_y"=midpoints_y,
                 "F_mat"=F_mat,"hist_data"=hist_data,
                 "R_opt"=R_opt,"Z_opt"=Z_opt,"v_opt"=v_opt,"u_opt"=u_opt,
                 "Z_x_bar"=Z_x_bar,"Z_y_bar"=Z_y_bar,"bbZ_bar"=bbZ_bar,"G"=G,"spline_hist"=spline_hist,
                 "rho"=rho,"alfa"=alfa,
                 "sd"=simp_d,"rsd"=rel_simp_d),class="bivariate_zbSpline")
}

plot.bivariate_zbSpline<-function(Z,type="static",what="full",scale="clr",title="",plot_hist=FALSE,plot=TRUE,xlab="x",ylab="y",xlim=NA,ylim=NA){
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
        theta = 325,          
        phi = 30,            
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
          col = "black",
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


#We still need a cross_validate function

cross_validate2d<-function(x,y=NULL,alfa_seq=seq(0.01,0.99,length.out=20),rho=NULL,res=100,k=2,l=2,u=1,v=1,knots_x_inner,knots_y_inner,bin_selection=doane){
  if(!is.null(rho)){
    alfa_seq<-1/(rho+1)
  }
  cv<-numeric(length(alfa_seq)) #Will capture the cross-validations scores
  
  #We make this code just to have the dimensions to make the permutation matrix from.
  #It is very inefficient, but it works so long as no coefficients are repeated
  Z_temp<-bivariate(x,y,alfa=0.5,bin_selection=bin_selection,knots_x_inner=knots_x_inner,knots_y_inner=knots_y_inner,k=k,l=l,u=u,v=v,res=5)
  csR<-as.vector(Z_temp$R_opt)
  csZ<-c(as.vector(Z_temp$Z_opt),Z_temp$v_opt,Z_temp$u_opt)
  P<-matrix(0,nrow=length(csR),ncol=length(csZ))#We make the permutation matrix
  for(i in 1:nrow(P)){
    sigma_i<-which(csZ==csR[i])
    P[i,sigma_i]<-1
  }
  
  #Compute the cross-validation scores at each alfa
  for(i in seq_along(alfa_seq)){
  Z<-bivariate(x,y,alfa_seq[i],knots_x_inner=knots_x_inner,knots_y_inner=knots_y_inner,bin_selection=bin_selection,k=k,l=l,u=u,v=v,res=res)
  H<-kronecker(t(Z$Z_y_bar),t(Z$Z_x_bar))%*%P%*%solve(Z$G)%*%Z$bbZ_bar #Definition of H
  m<-length(Z$midpoints_x)
  n<-length(Z$midpoints_y)
  cv[i]<-mean((Z$F_mat-Z$spline_hist)^2)/((1-sum(diag(H))/(m*n))^2)
  }
  structure(list("alfa" = alfa_seq, "cv" = cv, "optimum" = alfa_seq[which.min(cv)]),class="zbCV")
}

perm_test<-function(x,y,kx,ky,alfa=1,bin_selection=doane,k=2,l=2,u=1,v=1,res=100,K=200){
  rsd_true<-bivariate(x,y,alfa=alfa,bin_selection=bin_selection,knots_x_inner=kx,knots_y_inner=ky,k=k,l=l,u=u,v=v,res=res)$rsd
  rsd_vec<-numeric(K)
  for (i in 1:K) {
    repeat {  # Use repeat to retry in case of an error
      # Try-catch block
      tryCatch({
        # Sample data
        x_new <- sample(x, replace = FALSE)
        y_new <- sample(y, replace = FALSE)
        
        # Define knot sequences
        kx_new <- seq(min(x_new), max(x_new), length.out = 4)
        ky_new <- seq(min(y_new), max(y_new), length.out = 4)
        
        # Run the bivariate function and store result
        sim_rsd <- bivariate(x_new, y_new, alfa = alfa, bin_selection = bin_selection, 
                             knots_x_inner = kx_new, knots_y_inner = ky_new, 
                             k = k, l = l, u = u, v = v, res = res)
        rsd_vec[i] <- sim_rsd$rsd  # Save result
        
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
  return(p)
}


dist_test<-function(x,y,f,samp_func,scale="clr",k=3,l=3,u=2,v=2,alfa=0.5,g=2,h=2,bin_selection=doane){
  if (!scale %in% c("clr", "density")) {
    stop("Error: 'scale' must be either 'clr' or 'density'")
  }
  if(scale=="density"){ #Ensure the function is a clr density
    f<-clr2d(f,x,y)
  }
  #compute the empirical clr-density
  kx<-seq(min(x),max(x),length.out=g+2)
  ky<-seq(min(y),max(y),length.out=h+2)
  fit<-bivariate(x,y,alfa=alfa,knots_x_inner=kx,knots_y_inner=ky,k=k,l=l,u=u,v=v)
  #compute the hypothetical at same domain
  x_vals<-fit$seq_x
  y_vals<-fit$seq_y
  hyp_val<-outer(x_vals,y_vals,FUN=f)
  
  hyp_dist<-trapz2d(x_vals,y_vals,(fit$Z_spline-hyp_val)^2) #Compute the integral over the squared difference
}
