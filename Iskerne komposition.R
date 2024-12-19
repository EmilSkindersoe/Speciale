source("~/Universitet-mat/.Speciale/R/splines_combined.R")
library(readxl)
ice_ages<-read_excel("Rasmussen_et_al_2014_QSR_Table_2.xlsx",range = "A22:E139")[-(1:8),]
ice_ages<-ice_ages[-(97:nrow(ice_ages)),]
#GS is cold and GI is milder



ice_cores<-read_excel("GICC05modelext_GRIP_and_GISP2_and_resampled_data_series_Seierstad_et_al._2014_version_10Dec2014-2.xlsx", 
                sheet = "3) d18O and Ca 20 yrs mean", 
                range = "A51:M12500")[-1,]

ice_cores<-read_excel("GICC05modelext_GRIP_and_GISP2_and_resampled_data_series_Seierstad_et_al._2014_version_10Dec2014-2.xlsx", 
                                 sheet = "4) d18O and Ca 50 yrs mean", 
                                 range = "A5:M4894")[-1,]

ngrip2<-na.omit(as.data.frame(lapply(ice_cores[,c(1,5:6)],as.numeric),col.names = c("age","d18O","Ca2")))
grip<-na.omit(as.data.frame(lapply(ice_cores[,c(1,8:9)],as.numeric),col.names = c("age","d18O","Ca2")))
gisp<-na.omit(as.data.frame(lapply(ice_cores[,c(1,11:12)],as.numeric),col.names = c("age","d18O","Ca2")))
ngrip2$Ca2<-log(ngrip2$Ca2)
grip$Ca2<-log(grip$Ca2)
gisp$Ca2<-log(gisp$Ca2)


#We will be working only in the pre-holoscene period
process_ice_age_data <- function(data) {
  combined<-data %>%
    mutate(
      # Generalize period labels by capturing only the main part (e.g., GI-1, GS-2)
      Event = sub("^(Start of [A-Z]+-\\d+).*", "\\1", Event)
    ) %>%
    # Select only the Age column and the simplified Event column
    select(age = `Age (a b2k)`, period = Event)
  aggr<-aggregate(age ~ period, data = combined, FUN = max)
  return(aggr[order(aggr$age),])
}
ice_age_proc<-process_ice_age_data(ice_ages)


split_ice<-function(ice_ages,data){
  ice_times<-ice_ages$age
  periods<-ice_ages$period
  subsets_list<-list()
  j<-1
  for(i in length(ice_times):2){
    interval<-c(ice_times[i],ice_times[i-1])
    sub_data<-subset(data,age<interval[1]&age>interval[2])
    if(nrow(sub_data)==0){stop(paste("No data in interval",interval))}
      subsets_list[[j]]<-sub_data
      names(subsets_list)[j]<-periods[i]
      j<-j+1
  }
  subsets_list[[j]]<-subset(data,age<ice_times[i]) #To handle the last period because it is just all data before
  names(subsets_list)[j]<-periods[i-1]
  return(subsets_list)
}
ice_ages<-ice_age_proc
data<-ngrip2
ngrip2_list<-split_ice(ice_age_proc,ngrip2)
grip_list<-split_ice(ice_age_proc,grip)
gisp_list<-split_ice(ice_age_proc,gisp)

#We cross-validate over all observations from all cores
#x<-c(ngrip2$d18O,grip$d18O,gisp$d18O)
#y<-c(ngrip2$Ca2,grip$Ca2,gisp$Ca2)


#Cold periods
nrow(do.call(rbind,ngrip2_list[seq(2,length(ngrip2_list),by=2)]))
nrow(do.call(rbind,grip_list[seq(2,length(ngrip2_list),by=2)]))
nrow(do.call(rbind,gisp_list[seq(2,length(ngrip2_list),by=2)]))

#Mild periods
nrow(do.call(rbind,ngrip2_list[seq(1,length(ngrip2_list),by=2)]))
nrow(do.call(rbind,grip_list[seq(1,length(ngrip2_list),by=2)]))
nrow(do.call(rbind,gisp_list[seq(1,length(ngrip2_list),by=2)]))


ngrip_cold<-do.call(rbind,ngrip2_list[seq(2,length(ngrip2_list),by=2)])
ngrip_cold<-do.call(rbind,ngrip2_list[seq(2,length(ngrip2_list),by=2)])
ngrip_cold<-do.call(rbind,ngrip2_list[seq(2,length(ngrip2_list),by=2)])

#Mild periods

x<-c(ngrip2$d18O)
y<-c(ngrip2$Ca2)

kx<-seq(min(x),max(x),length.out=4)
ky<-seq(min(y),max(y),length.out=4)

#cross<-cross_validate2d(x,y,knots_x_inner=kx,knots_y_inner=ky,bin_selection=doane,k=3,l=3,u=1,v=1)
#plot(cross)
#cross$optimum

biv_full<-bivariate(x,y,alfa=0.9,knots_x_inner=kx,knots_y_inner=ky,k=3,l=3,u=1,v=1)
plot(biv_full,scale="density",plot_hist=TRUE,type="interactive",title="density function of all cores and all observations",xlab="δ18O",ylab="log[Ca2+]")

par(mfrow=c(2,3),mar=c(0,0,1,0))

plot(biv_full,what="full",scale="clr",plot_hist=TRUE,type="static",title="clr full",xlab="δ18O",ylab="log[Ca2+]")
plot(biv_full,what="interaction",scale="clr",type="static","clr interaction",xlab="δ18O",ylab="log[Ca2+]")
plot(biv_full,what="independent",scale="clr",type="static","clr independent",xlab="δ18O",ylab="log[Ca2+]")
plot(biv_full,what="full",scale="density",plot_hist=TRUE,type="static",title="density full",xlab="δ18O",ylab="log[Ca2+]")
plot(biv_full,what="interaction",scale="density",type="static","density interaction",xlab="δ18O",ylab="log[Ca2+]")
plot(biv_full,what="independent",scale="density",type="static","density independent",xlab="δ18O",ylab="log[Ca2+]")

biv_full$rsd




x<-ngrip2_list[[44]]$d18O #Odd numbers are mild periods, evens are ice ages
y<-ngrip2_list[[44]]$Ca2
kx<-seq(min(x),max(x),length.out=4)
ky<-seq(min(y),max(y),length.out=4)
#cross<-cross_validate2d(x,y,knots_x_inner=kx,knots_y_inner=ky,bin_selection=doane)
#plot(cross)
#Det kunne godt være at vi skal medtage en periode kun hvis dens histogramdata faktisk kan fittes
#Og vi skal måske finde et optimalt alpha til istider og varmeperioder

par(mfrow=c(1,2),mar=c(0,1,1,0))
biv_cold<-bivariate(x,y,knots_x_inner = kx,knots_y_inner=ky,bin_selection=doane,alfa=0.9,k=3,l=3,u=1,v=1)
x<-ngrip2_list[[45]]$d18O #Odd numbers are mild periods, evens are ice ages
y<-ngrip2_list[[45]]$Ca2
kx<-seq(min(x),max(x),length.out=4)
ky<-seq(min(y),max(y),length.out=4)
biv_mild<-bivariate(x,y,knots_x_inner = kx,knots_y_inner=ky,bin_selection=doane,alfa=0.9,k=3,l=3,u=1,v=1)
plot(biv_cold,scale="density",what="full",type="static",plot_hist=FALSE,title="(GS-2, NGRIP2)",xlab="δ18O",ylab="log[Ca2+]")
plot(biv_mild,scale="density",what="full",type="static",plot_hist=FALSE,title="(GI-1, NGRIP2)",xlab="δ18O",ylab="log[Ca2+]")


ngrip_ice<-do.call(rbind,ngrip2_list[seq(1,length(ngrip2_list),by=2)])
x<-ngrip2_list[[44]]$d18O #Even numbers are mild periods, odds are ice ages
y<-ngrip2_list[[44]]$Ca2
kx<-seq(min(x),max(x),length.out=4)
ky<-seq(min(y),max(y),length.out=4)
#cross<-cross_validate2d(x,y,knots_x_inner=kx,knots_y_inner=ky,bin_selection=doane)
#plot(cross)

biv_ice<-bivariate(x,y,knots_x_inner = kx,knots_y_inner=ky,bin_selection=doane,alfa=1,k=3,l=3,u=1,v=1)
plot(biv_ice,scale="density",what="full",type="static",plot_hist=FALSE,title="density function during cold period (GS-2, NGRIP2)",xlab="δ18O",ylab="log[Ca2+]")
#The distributions during different ages are rather dissimilar, so it would be unfair to show samples

#ngrip2_rsd<-data.frame("rsd"=rep(NA,length(ngrip2_list)),"period"=rep(c("cold","mild"),length.out=length(ngrip2_list)))
#grip_rsd<-data.frame("rsd"=rep(NA,length(grip_list)),"period"=rep(c("cold","mild"),length.out=length(grip_list)))
#gisp_rsd<-data.frame("rsd"=rep(NA,length(gisp_list)),"period"=rep(c("cold","mild"),length.out=length(gisp_list)))

View(ngrip2_list)
fill_rsd<-function(list){
  ice_rsd<-data.frame("rsd"=NA,"period"=rep(c("mild","cold"),length.out=length(list))) #We start in a mild period
  for(i in 1:nrow(ice_rsd)){
    x<-list[[i]]$d18O
    y<-list[[i]]$Ca2
    kx<-seq(min(x),max(x),length.out=4)
    ky<-seq(min(y),max(y),length.out=4)
    cross<-tryCatch({
      cross_validate2D(x,y,knots_x_inner=kx,knots_y_inner=ky,bin_selection=doane,k=3,l=3,u=1,v=1)},
       error = function(e) {
       message("Skipping iteration ", i, " due to error: ", conditionMessage(e))
       NULL # Return NULL if there's an error
      })
    if(!is.null(cross)){
      dens <- tryCatch({
        zbSpline2D(x, y, alfa = cross$optimum, bin_selection = doane, knots_x_inner = kx, 
                  knots_y_inner = ky, k = 3, l = 3, u = 1, v = 1)
      })
      
      # Only assign if dens was successfully computed
      if (!is.null(dens)) {
        ice_rsd$rsd[i] <- dens$rsd
      }
    }
  }
  return(ice_rsd)
}
ngrip2_box<-fill_rsd(ngrip2_list)
ngrip2_box$source="NGRIP2"
grip_box<-fill_rsd(grip_list)
grip_box$source="GRIP"
gisp_box<-fill_rsd(gisp_list)
gisp_box$source="GISP"

combined_box<-rbind(ngrip2_box,grip_box,gisp_box)
combined_box$source=factor(combined_box$source,levels=c("NGRIP2","GRIP","GISP"))
ggplot(combined_box, aes(x = period, y = rsd, fill = period)) +
  geom_boxplot() +
  facet_wrap(~ source) +
  labs(
    title = "50-Year Average",
    x = "Period",
    y = "RSD"
  ) +
  theme_minimal() +
  ylim(0, 1) +
  scale_fill_manual(
    values = c("mild" = "#f8766d", "cold" = "#00bfc4"),
    name = "Period"
  )

t.test(subset(combined_box,period=="cold")$rsd,subset(combined_box,period=="mild")$rsd)

#We aggregate all cold periods and all mild periods and get p-values for H0 in each

ngrip2_ice<-do.call(rbind,ngrip2_list[seq(1,length(ngrip2_list),by=2)])
grip_ice<-do.call(rbind,grip_list[seq(1,length(ngrip2_list),by=2)])
gisp_ice<-do.call(rbind,gisp_list[seq(1,length(ngrip2_list),by=2)])
ice_sim<-rbind(ngrip2_ice,grip_ice,gisp_ice)
x<-ice_sim$d18O
y<-ice_sim$Ca2
kx<-seq(min(x),max(x),length.out=4)
ky<-seq(min(y),max(y),length.out=4)
sim_fit<-(zbSpline2D(x,y,knots_x_inner=kx,alfa=0.1,knots_y_inner=ky,k=3,l=3,u=1,v=1))
summary(sim_fit)
p_ice<-perm_test(x,y,kx=kx,ky=ky,alfa=0.1,k=3,l=3,u=1,v=1,K=1000)

dhsic.test(x,y)

#Make a histogram from it
hist_rsd<-hist(p_ice$rsd_perms,plot=FALSE)

cross_validate1D(data.frame(hist_rsd$mids,hist_rsd$counts),knots_inner=c(min(hist_rsd$mids),0.01,0.02,0.03))
plot_rsd<-zbSpline1D(data.frame(hist_rsd$mids,hist_rsd$counts),alfa=1,knots_inner=c(min(hist_rsd$mids),0.01,0.02,max(hist_rsd$mids)))
plot(plot_rsd,what="Z-spline")


hist(p_ice$rsd_perms,probability=TRUE,main="Histogram of simulated RSD statistics",xlab="RSD")
abline(v=p_ice$rsd_true,lty=2)
lines(plot_rsd$x_seq,exp(plot_rsd$Z_spline)/trapz(plot_rsd$x_seq,exp(plot_rsd$Z_spline)))



ngrip2_mild<-do.call(rbind,ngrip2_list[seq(2,length(ngrip2_list),by=2)])
grip_mild<-do.call(rbind,grip_list[seq(2,length(ngrip2_list),by=2)])
gisp_mild<-do.call(rbind,gisp_list[seq(2,length(ngrip2_list),by=2)])
mild_sim<-rbind(ngrip2_mild,grip_mild,gisp_mild)
x<-mild_sim$d18O
y<-mild_sim$Ca2
kx<-seq(min(x),max(x),length.out=4)
ky<-seq(min(y),max(y),length.out=4)
p_mild<-perm_test(x,y,kx=kx,ky=ky,alfa=0.1,k=3,l=3,u=1,v=1,K=1000)
p_mild

dhsic.test(x,y,method="gamma")


par(mfrow=c(1,1))

#To compare different bivariate estimators
x<-c(subset(ngrip2,age>11700 & age<104000)$d18O,subset(grip,age>11700 & age<104000)$d18O,subset(gisp,age>11700 & age<104000)$d18O)
y<-c(subset(ngrip2,age>11700 & age<104000)$Ca2,subset(grip,age>11700 & age<104000)$Ca2,subset(gisp,age>11700 & age<104000)$Ca2)
kx<-seq(min(x),max(x),length.out=10)
ky<-seq(min(y),max(y),length.out=10)
biv_fit<-bivariate(x,y,alfa=0.9,bin_selection=scott,knots_x_inner=kx,knots_y_inner=ky,k=3,l=3,u=1,v=1)
plot(biv_fit,scale="density",title="C-spline density",xlab="δ18O",ylab="log[Ca2+]",plot_hist=TRUE) #Vi kan faktisk genkende de kolde og milde perioder her
#Samtidig kan vi se at den kolde bin faktisk nok fittes lidt unaturligt højt
dens<-kde2d(x,y,n=200)
persp3D(x = dens$x, 
                y = dens$y, 
                z = dens$z, 
                col = viridis(50),
                theta = 325,          
                phi = 30,            
                ticktype = "detailed", 
                nticks = 3,
                xlab="δ18O",
                ylab="log[Ca2+]",          
                zlab = "",          
                bty = "b2",
                colkey = FALSE,
                main = "KDE",
                cex.axis = 0.5,
                zlim=c(0,0.21))


x_centers<-biv_fit$midpoints_x
y_centers<-biv_fit$midpoints_y
grid<-expand.grid("X"=x_centers,"Y"=y_centers)



grid$Z<-c(biv_fit$hist_data/trapz2d(x_centers,y_centers,biv_fit$hist_data))
fit<-gam(Z~te(X,Y,bs="cr",k=c(10,10),m=c(3,2)),data=grid)
x_pred <- seq(min(x_centers), max(x_centers), length.out = 100)
y_pred <- seq(min(y_centers), max(y_centers), length.out = 100)
pred_grid <- expand.grid(X = x_pred, Y = y_pred)
pred_grid$Z <- predict(fit, newdata = pred_grid)

# Reshape the prediction data for persp3D
z_matrix <- matrix(pred_grid$Z, nrow = length(x_pred), ncol = length(y_pred))
persp3D(
  x = x_pred, 
  y = y_pred, 
  z = z_matrix, 
  col = viridis(50),        # Use viridis for color gradient
  theta = 325,              # Rotation angle
  phi = 30,                 # Viewing angle
  ticktype = "detailed", 
  nticks = 3,               # Number of ticks
  xlab = "δ18O",            # X-axis label
  ylab = "log[Ca2+]",       # Y-axis label
  zlab = "",                # Z-axis label
  bty = "b2",               # Box type
  colkey = FALSE,           # Disable color legend
  main = "B-spline density",# Title
  cex.axis = 0.5,
  zlim=c(0,0.21)
)
trapz2d(x_pred,y_pred,z_matrix)



#We plot the data
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("rnaturalearth", quietly = TRUE)) install.packages("rnaturalearth")
if (!requireNamespace("rnaturalearthdata", quietly = TRUE)) install.packages("rnaturalearthdata")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)

# Define coordinates and labels with manual nudges
sites <- data.frame(
  name = c("NGRIP2", "GISP", "GRIP"),
  lat = c(75.10, 72.58, 72.58),
  lon = c(-42.32, -38.48, -37.64),
  nudge_x = c(0, -5, 5),  # Adjust horizontal position
  nudge_y = c(0.7, 0.7, 0.7) # Adjust vertical position
)

# Load Greenland map and filter specifically for Greenland
world <- ne_countries(scale = "medium", returnclass = "sf")
greenland <- world %>% filter(admin == "Greenland")

# Create the map
ggplot(data = greenland) +
  geom_sf(fill = "lightblue", color = "black") +
  geom_point(data = sites, aes(x = lon, y = lat), color = "purple", size = 1) +
  geom_text(data = sites, aes(x = lon + nudge_x, y = lat + nudge_y, label = name), 
            hjust = 0.5, vjust = 0.5) +
  coord_sf(xlim = c(-75, -10), ylim = c(60, 85), expand = FALSE) +
  theme_minimal() +
  labs(
    title = "Greenland Map with Selected Sites",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(panel.grid.major = element_line(color = "gray80", linetype = "dotted"))
