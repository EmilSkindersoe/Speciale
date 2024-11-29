source("splines_combined.R")
library(openxlsx)
#We read all the data
#Per construction of the functions, a row is an observation from x and a column is an observation from y.
#All integrals are divided by the domain which they are taken over, in order to create a probability base measure


fin20<-read.xlsx("frbny-sce-public-microdata-latest.xlsx",colNames=TRUE,rows=2:60084,cols=c(1,2,14,16,19,117,119,122))
fin17<-read.xlsx("FRBNY-SCE-Public-Microdata-Complete-17-19.xlsx",colNames=TRUE,rows=2:48763,cols=c(1,2,12,14,17,108,110,113))
fin13<-read.xlsx("FRBNY-SCE-Public-Microdata-Complete-13-16.xlsx",colNames=TRUE,rows=2:56446,cols=c(1,2,12,14,17,108,110,113))

fin_total<-rbind(fin13,fin17,fin20)
fin_total<-na.omit(fin_total)
rm(fin13,fin17,fin20) #They take up too much memory
gc()

x<-fin_total$Q9_cent50
y<-fin_total$C1_cent50

kx<-c(-36,seq(-25,25,by=10),36)->ky

c<-cross_validate2D(x,y,alfa_seq=seq(0.0001,0.0008,length.out=10),bin_selection=doane,knots_x_inner=kx,knots_y_inner=ky,k=3,l=3,u=2,v=2)
plot(c)
c$optimum

test<-zbSpline2D(x,y,alfa=0.0006,knots_x_inner=kx,knots_y_inner=ky,k=3,l=3,u=2,v=2)
par(mfrow=c(1,1),mar=c(0,0,1,0))
plot(test,what="geom_Y",scale="density",type="static",theta=125,xlab="inflation",ylab="housing price",title="C-spline density")
test$midpoints_y
#Prøv også med kde2d
dens<-kde2d(x,y,n=100)
plot_ly(x=dens$x,y=dens$y,z=dens$z,type="surface")
persp3D(x = dens$x, 
        y = dens$y, 
        z = dens$z, 
        col = viridis(50),
        theta = 125,          
        phi = 30,            
        ticktype = "detailed", 
        nticks = 3,
        xlab="inflation",
        ylab="housing price",          
        zlab = "",          
        bty = "b2",
        colkey = FALSE,
        main = "KDE",
        cex.axis = 0.5)
hist_data<-test$hist_data/trapz2d(test$midpoints_x,test$midpoints_y,test$hist_data)
grid<-expand.grid("x"=test$midpoints_x,"y"=test$midpoints_y)
points3D(
  x = grid$x,
  y = grid$y,
  z = as.vector(hist_data),
  add = TRUE,
  col = "black",
  pch = 19,
  cex = 0.5
)


par(mfrow=c(3,2),mar=c(0,0,1,0))

plot(test,what="full",scale="clr",plot_hist=TRUE,type="static",title="clr full",xlab="inflation",ylab="housing p")
plot(test,what="full",scale="density",plot_hist=TRUE,type="static",title="density full",xlab="inflation",ylab="housing p")
plot(test,what="independent",scale="clr",type="static","clr independent",xlab="inflation",ylab="housing p")
plot(test,what="independent",scale="density",type="static","density independent",xlab="inflation",ylab="housing p")
plot(test,what="interaction",scale="clr",type="static","clr interaction",xlab="inflation",ylab="housing p")
plot(test,what="interaction",scale="density",type="static","density interaction",xlab="inflation",ylab="housing p")

perm_test(fin_total$Q9_cent50,fin_total$C1_cent50,kx,ky,alfa=0.00002)

#We plot the marginals of the entire data as well as their distances
seq_x<-test$seq_x


geom_x<-test$C_spline_x
arithmetic_x<-colSums(test$C_spline)/trapz(seq_x,colSums(test$C_spline))
clr_geom_x<-log(geom_x)-trapz(seq_x,log(geom_x))/(max(seq_x)-min(seq_x))
clr_arithmetic_x<-log(arithmetic_x)-trapz(seq_x,log(arithmetic_x))/(max(seq_x)-min(seq_x))
norm1<-trapz(seq_x,geom_x/arithmetic_x)
norm2<-trapz(seq_x,geom_x-arithmetic_x)
trapz(seq_x,(geom_x-arithmetic_x)/norm2)

plot_data<-data.frame("x"=seq_x,"geom_X"=geom_x,"arith_X"=arithmetic_x,"perturbed_diff"=(geom_x/arithmetic_x)/norm1)
pivot_data<-pivot_longer(plot_data,cols=!x)
ggplot(pivot_data)+geom_line(mapping=aes(x=x,y=value,color=name))

#Make a list of dataframes 

total_months<-unique(fin_total$date)
list_months<-vector("list",length=length(total_months))
range_data<-data.frame("min_x"=rep(NA,length(total_months)),"max_x"=NA,"min_y"=NA,"max_y"=NA)
for(i in 1:length(list_months)){
  new_data<-subset(fin_total,date==total_months[i])
  list_months[[i]]<-new_data[,c(1,2,4,7)] #Take date, userid, median inflation and median housing price
  range_data$min_x[i]<-min(new_data$Q9_cent50)
  range_data$max_x[i]<-max(new_data$Q9_cent50)
  range_data$min_y[i]<-min(new_data$C1_cent50)
  range_data$max_y[i]<-max(new_data$C1_cent50)
}
x_range<-c(max(range_data$min_x),min(range_data$max_x))
y_range<-c(max(range_data$min_y),min(range_data$max_y))


#Make a list of splines
list_of_splines<-vector("list",length=length(total_months))

pb <- txtProgressBar(min = 0, max = length(list_of_splines), style = 3)
for(i in 1:length(list_of_splines)){
  setTxtProgressBar(pb, i)
  data<-list_months[[i]]
  x<-data$Q9_cent50
  y<-data$C1_cent50
  kx<-seq(min(x),max(x),length.out=10)
  ky<-seq(min(y),max(y),length.out=10)
  c<-cross_validate2D(x,y,kx,ky,alfa_seq=seq(0.00001,0.1,length.out=20),k=3,l=3,u=2,v=2) #Setting a range here is bugged for some reason
  alfa<-c$optimum
  list_of_splines[[i]]<-zbSpline2D(x,y,kx,ky,alfa=alfa,range_x=x_range,range_y=y_range)
}

saveRDS(list_of_splines,"splines_list")
list_of_splines<-readRDS("splines_list")

#We plot RSD as a function of time
rsd_data<-data.frame("month_along"=seq(0,length(list_of_splines)-1),
                     "month"=NA,
                     "RSD"=NA)

convert_to_date <- function(date_char) {
  Sys.setlocale("LC_TIME", "C")
  date_obj <- as.Date(paste0(date_char, "01"), format = "%Y%m%d")
  formatted_date <- format(date_obj, "%B %Y")
  return(tolower(formatted_date))
}
rsd_data$month=convert_to_date(total_months)

for(i in 1:nrow(rsd_data)){
  rsd_data$RSD[i]<-list_of_splines[[i]]$rsd
}

lin_mod<-lm(RSD~month_along,data=rsd_data)
summary(lin_mod)
ggplot(rsd_data, aes(x = month_along, y = RSD)) +
  geom_point(color="#619CFF") +       # Points on the line
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +  # Regression line
  scale_x_continuous(
    breaks = seq(0, max(rsd_data$month_along), by = 10), # Show every 10th month_along
    labels = rsd_data$month[seq(1, nrow(rsd_data), by = 10)] # Custom labels
  ) +
  geom_vline(xintercept=105)+
  labs(x = "Months", y = "RSD") + # Labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels



integrate_over_x <- function(fx, x_values) {
  trapz(x_values, fx)
}

spline_test<-list_of_splines[[1]]
plot(colSums(spline_test$C_spline),main="arithmetic marginal") #colsums takes the sums over all columns for a row
plot(spline_test$C_spline_x,main="geometric marginal")

#This function help with integrating out some axis
integrate_over_axis <- function(M, x_values, margin = 1, method = "trapz") {
  # M: Matrix of function values
  # x_values: Vector of variable values corresponding to the margin
  # margin: 1 for rows (integrate over x), 2 for columns (integrate over y)
  # method: Numerical integration method ("trapz", "simpson")
  
  # Load necessary package
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("Package 'pracma' is required. Please install it using install.packages('pracma').")
  }

  if (length(x_values) != dim(M)[margin]) {
    stop("Length of x_values must match the dimension of M over which to integrate.")
  }
  
  # Select integration function based on the specified method
  integration_function <- switch(method,
                                 trapz = pracma::trapz,
                                 simpson = pracma::simps,
                                 stop("Unsupported integration method. Use 'trapz' or 'simpson'."))
  
  # Define the integration function to apply
  integrate_func <- function(f_values) {
    integration_function(x_values, f_values)
  }
  
  # Determine the MARGIN parameter for apply (opposite of the integration axis)
  apply_margin <- ifelse(margin == 1, 2, 1)
  
  # Apply the integration function over the specified axis
  result <- apply(M, MARGIN = apply_margin, integrate_func)
  
  return(result)
}

new_clr_marg<-integrate_over_axis(spline_test$Z_spline_int,x_seq,margin=2)

plot(x_seq,new_clr_marg/max(new_clr_marg))
plot(x_seq,spline_test$Z_spline_x)
plot(x_seq,log(colSums(spline_test$C_spline)))
#Conclusion: Integrating over the independent part gives the same clr-marginal
#as the constructed spline. It is not the same as the defined clr-marginal.
#We choose to recalculate everything from the new Z_splines which come from the permutation (for now).


#We can make a function which takes as input a list of splines and then returns the test statistic
#Make another function which takes as input a list of splines and returns the permuted splines.

#Since we evaluate on the same domain for all functions, these are saved globally.
x_seq<-list_splines[[1]]$seq_x
y_seq<-list_splines[[1]]$seq_y
norm_cons_x<-(x_range[2]-x_range[1])
norm_cons_y<-(y_range[2]-y_range[1])

list_splines<-lapply(list_of_splines,function(df){df$Z_spline})
calc_test<-function(list_splines){
  #Takes a list of clr-spline matrices and computes a test-statistic
  t_i<-vector("numeric",length=length(list_splines))
  for(i in 1:length(t_i)){
    sim_Z<-list_splines[[i]]
    sim_dens<-exp(sim_Z)/trapz2d(x_seq,y_seq,exp(sim_Z))
    a_marginal_x<-rowSums(sim_dens)
    a_marginal_y<-colSums(sim_dens)
    a_marginal_x_clr<-log(a_marginal_x)-(1/norm_cons_x)*trapz(x_seq,log(a_marginal_x))
    a_marginal_y_clr<-log(a_marginal_y)-(1/norm_cons_y)*trapz(y_seq,log(a_marginal_y))
    
    #This is not the same as the clr-marginal according to the spline because we integrate over the entire spline
    #and not just the independent part
    g_marginal_x_clr<-integrate_over_axis(sim_Z,x_seq,margin=2)/(norm_cons_x) 
    g_marginal_y_clr<-integrate_over_axis(sim_Z,y_seq,margin=1)/norm_cons_y
    g_marginal_x<-exp(g_marginal_x_clr)/trapz(x_seq,exp(g_marginal_x_clr))
    D_X<-sqrt(trapz(x_seq,(g_marginal_x_clr-a_marginal_x_clr)^2))
    D_Y<-sqrt(trapz(x_seq,(g_marginal_y_clr-a_marginal_y_clr)^2))
    t_i[i]<-D_X+D_Y
  }
  return(sum(t_i))
}
T_obs<-calc_test(list_splines)  

#list_of_splines are of spline objects
perm_spline<-function(list_of_splines){
  perm<-sample(1:length(list_of_splines),replace=FALSE)
  list_perm<-vector("list",length=length(perm))
  for(i in 1:length(perm)){
    list_perm[[i]]<-list_of_splines[[i]]$Z_spline_ind+list_of_splines[[perm[i]]]$Z_spline_int
  }
  return(list_perm)
}

K<-500
T_new<-rep(NA,length=K)
pb <- txtProgressBar(min = 0, max = K, style = 3)
for(i in 1:K){
  setTxtProgressBar(pb, i)
  perm_new<-perm_spline(list_of_splines)
  T_new[i]<-calc_test(perm_new)
}
p<-length(which(T_new>=T_obs))/K
ggplot()+
  geom_histogram(data=data.frame("T"=T_new),mapping=aes(x=T),bins=doane(T_new),fill="gray",color="black")+
  theme_minimal()

hist_T<-hist(T_new,plot=FALSE)

plot_T<-zbSpline1D(data.frame(hist_T$mids,hist_T$counts))

hist(T_new,xlim=c(min(hist_T$breaks),T_obs+1),probability=T,main="Histogram of simulated T statistics",xlab="T")
abline(v=T_obs,lty=2)
lines(plot_T$x_seq,plot_T$C_spline)
