library(easyCODA)
library(DirichletReg)
library(vcd)
library(Ternary)
library(ggplot2)
library(ggtern)
library(robCompositions)
library(ggtext)
library(tidyverse)
library(viridis)
library(ggtext)
library(gridExtra)


ArcticLake<-ArcticLake
ArcticLake[24,3]<-0.405
ArcticLake[4,3]<-0.069

plot(AL, cex = 0.5, a2d = list(colored = FALSE, c.grid = FALSE))

par(mfrow=c(1,2))

BAR(ArcticLake[,1:3]*100,
    cols=c("#DAA520","#FF7F50","#008080"),
    eps=0,
    row.names=signif(ArcticLake$depth,2),
    main="Soil composition (depth in meters)",
    ylab="")

TernaryPlot(grid.lines=4,
            atip=colnames(ArcticLake)[1],
            btip=colnames(ArcticLake)[2],
            ctip=colnames(ArcticLake)[3],
            tip.col=c("#DAA520","#FF7F50","#008080"),
            main="Ternary plot of soil data")
TernaryPoints(coordinates=ArcticLake[,1:3],pch=16,cex=0.7)
#TernaryPoints(coordinates=ArcticLake[,1:3]+c(1/3,1/3,1/3),pch=16,cex=0.7,col="darkred")
#TernaryPoints(coordinates=ArcticLake[,1:3]*10,pch=16,cex=0.7,col="lightblue")

library(proxy)

M<-as.matrix(ArcticLake[,1:3])

dist_matrix<-dist(M[,1:3],method="euclidean")
dist_matrix <- as.matrix(dist_matrix)


min_diff <- Inf
best_pair1 <- c()
best_pair2 <- c()

for (i in 1:(nrow(dist_matrix) - 1)) {
  for (j in (i + 1):nrow(dist_matrix)) {
    d1 <- dist_matrix[i, j]
    diff1_x <- abs(M[i, 1] - M[j, 1])
    
    for (k in 1:(nrow(dist_matrix) - 1)) {
      for (l in (k + 1):nrow(dist_matrix)) {
        if (i != k && j != l && i != l && j != k) {
          d2 <- dist_matrix[k, l]
          diff2_x <- abs(M[k, 1] - M[l, 1])
          
          # Find pairs with similar distances and larger first coordinate difference in the first pair
          if (abs(d1 - d2) < 0.1 && diff1_x > diff2_x * 1.5) {
            diff <- abs(d1 - d2)
            
            if (diff < min_diff) {
              min_diff <- diff
              best_pair1 <- c(i, j)
              best_pair2 <- c(k, l)
            }
          }
        }
      }
    }
  }
}

# Output the results
if (length(best_pair1) > 0 && length(best_pair2) > 0) {
  cat("First pair of rows:", best_pair1, "with coordinates", M[best_pair1[1],], "and", M[best_pair1[2],], "\n")
  cat("Second pair of rows:", best_pair2, "with coordinates", M[best_pair2[1],], "and", M[best_pair2[2],], "\n")
  cat("Their respective distances are:", dist_matrix[best_pair1[1], best_pair1[2]], "and", dist_matrix[best_pair2[1], best_pair2[2]], "\n")
  cat("The difference between the distances is:", min_diff, "\n")
  cat("Difference in the first coordinate for pair 1:", abs(M[best_pair1[1], 1] - M[best_pair1[2], 1]), "\n")
  cat("Difference in the first coordinate for pair 2:", abs(M[best_pair2[1], 1] - M[best_pair2[2], 1]), "\n")
} else {
  cat("No suitable pairs found.\n")
}

clr<-function(x){
  log(x)-1/length(x)*sum(log(x))
}

ArcticLake$Adist<-NA
ArcticLake$Edist<-NA
#Tager afstande til sidste punkt
for(i in 1:(nrow(ArcticLake))){
  ArcticLake$Adist[i]<-aDist(ArcticLake[i,1:3],ArcticLake[1,1:3])
  ArcticLake$Edist[i]<-sqrt(sum((ArcticLake[i,1:3]-ArcticLake[1,1:3])^2))
}


ArcticLake_Adist <- ArcticLake %>%
  arrange(Adist) %>%
  mutate(label_A = ifelse(row_number() <= 6 & row_number()>1, as.character(row_number()-1), NA))

ArcticLake_Edist <- ArcticLake %>%
  arrange(Edist) %>%
  mutate(label_E = ifelse(row_number() <= 6 & row_number()>1, as.character(row_number()-1), NA))

# Create the ternary plot for Adist with labels for 1 to 5
ternary_a <- ggtern(data = ArcticLake_Adist, mapping = aes(x = clay, y = sand, z = silt)) +
  theme(
    tern.axis.title.T = element_text(color = "#DAA520"),  # Golden color for the T axis label
    tern.axis.title.L = element_text(color = "#008080"),  # Teal color for the L axis label
    tern.axis.title.R = element_text(color = "#FF7F50"),
      plot.margin = unit(c(1, 1, 1, 1), "lines"),  # Adjust margins: top, right, bottom, left
      legend.position = "bottom",  # Move legend to bottom if it's taking up side space
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12)
  ) +
  geom_point(data=ArcticLake_Adist[-c(1:6),],mapping = aes(color = Adist)) +
  scale_color_viridis(name = "Aitchison") +  # Color scale for Adist
  labs(title = "", color = "Aitchison") +
  geom_point(data = ArcticLake_Adist[1,1:3], col = "red") +
  geom_text(aes(label = label_A), na.rm = TRUE,size=2.5)  # Add labels for the lowest Adist points

# Create the ternary plot for Edist with labels for 1 to 5
ternary_e <- ggtern(data = ArcticLake_Edist, mapping = aes(x = clay, y = sand, z = silt)) +
  theme(
    tern.axis.title.T = element_text(color = "#DAA520"),  # Golden color for the T axis label
    tern.axis.title.L = element_text(color = "#008080"),  # Teal color for the L axis label
    tern.axis.title.R = element_text(color = "#FF7F50"),
      plot.margin = unit(c(1, 1, 1, 1), "lines"),  # Adjust margins: top, right, bottom, left
      legend.position = "bottom",  # Move legend to bottom if it's taking up side space
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12)
  ) +
  geom_point(data=ArcticLake_Edist[-c(1:6),],mapping = aes(color = Edist)) +
  scale_color_viridis(name = "Euclidean",breaks=c(0.4,0.6,0.8)) +  # Color scale for Edist
  labs(title = "", color = "Euclidean") +
  geom_point(data = ArcticLake_Edist[1,1:3], col = "red") +
  geom_text(aes(label = label_E), na.rm = TRUE, size=2.5)  # Add labels for the lowest Edist points

# Combine the two plots into one
combined_plot <- grid.arrange(
  ternary_a, ternary_e, ncol = 2,
  top = textGrob("Plot of Distances", gp = gpar(fontface = "bold", fontsize = 14), hjust = 0.5, vjust = 8)
)
x1_cov<-matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)
x2_cov<-matrix(c(1,0.8,0,0.8,1,0,0,0,1),nrow=3)
x3_cov<-matrix(c(5,0,0,0,5,0,0,0,5),nrow=3)

x1<-data.frame(rmvnorm(10000,sigma=x1_cov))
x1_comp<-data.frame(clrInv(x1))

x2<-data.frame(rmvnorm(10000,sigma=x2_cov))
x2_comp<-data.frame(clrInv(x2))

x3<-data.frame(rmvnorm(10000,sigma=x3_cov))
x3_comp<-data.frame(clrInv(x3))

#Simulations sammenligning
heat_scatter<-function(df){
  ggplot(df, aes(x = X1, y = X2)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_gradientn(colors = c("black", "red", "yellow", "white")) + 
    theme(
      plot.background = element_rect(fill = "white", color = "white"),  # Set the background to black
      panel.background = element_rect(fill = "white", color = "white"),
      legend.position = "none"
    ) +
    labs(title = "",
         x = NULL,
         y = NULL)+
    xlim(c(-5,5))+
    ylim(c(-5,5))
}


heat_ternary<-function(df){
  ggtern(data = df, aes(x = X1, y = X2, z = X3)) +
    stat_density_tern(geom = 'polygon', aes(fill = ..level..), color = NA,bdl=0.05) +
    scale_fill_gradientn(colors = c("black", "red", "yellow", "white")) +
    theme(
      tern.panel.background = element_rect(fill = "black", color = NA),  # Black ternary background
      plot.background = element_rect(fill = "white", color = NA),  # White outer background
      panel.background = element_rect(fill = "white", color = NA),  # White panel background
      tern.axis.line = element_blank(),  # Remove axis lines
      tern.axis.text = element_blank(),  # Remove axis text
      tern.axis.title = element_blank(),  # Remove axis titles
      tern.panel.grid.major = element_blank(),  # Remove major grid lines
      tern.panel.grid.minor = element_blank(),  # Remove minor grid lines
      tern.plot.background = element_rect(fill = "white", color = NA), # Ensure only ternary plot background is black
      legend.position = "none"  # Hide the legend
    ) +
    theme_nolabels()+
    theme_hidegrid_major() +  # Hide major grid lines if present
    theme_hidegrid_minor()# Hide axis labels and titles
}

p1<-heat_scatter(x1)
p2<-heat_ternary(x1_comp)
p3<-heat_scatter(x2)
p4<-heat_ternary(x2_comp)
p5<-heat_scatter(x3)
p6<-heat_ternary(x3_comp)

grid.arrange(p1,p3,p5,p2,p4,p6,nrow=2)


ArcticLake_p<-acomp(ArcticLake,parts=1:3)
perturbation<-as.data.frame(perturbe(ArcticLake_p[,1:3],acomp(c(0.2,0.2,0.6))))
powering<-as.data.frame(power.acomp(ArcticLake_p,0.5))



ternary_p<-ggtern(mapping = aes(x = clay, y = sand, z = silt)) +
  theme(
    tern.axis.title.T = element_text(color = "#DAA520"),  # Golden color for the T axis label
    tern.axis.title.L = element_text(color = "#008080"),  # Teal color for the L axis label
    tern.axis.title.R = element_text(color = "#FF7F50")  # Centers the title horizontally# Coral color for the R axis label
  ) +
  geom_point(data=perturbation,col="cyan")+
  geom_point(data=ArcticLake,col="black")+
  geom_point(data=powering,col="pink")+
  ggtitle("Visualisation of <span style='color:cyan;'>perturbation</span> and <span style='color:pink;'>powering</span>")+
  theme(plot.title = element_markdown(hjust=0.5))


ternary_p

perturb_cont<-function(f,g){
  nom<-function(x){
    f(x)*g(x)
  }
  denom<-integrate(nom,lower=-6,upper=6)$value
  return(function(x){nom(x)/denom})
}

power_cont<-function(f,a=1){
  nom<-function(x){
    f(x)^a
  }
  denom<-integrate(nom,lower=-6,upper=6)$value
  return(function(x){nom(x)/denom})
}

axis_only_theme <- theme(
       panel.background = element_blank(),
       panel.grid = element_blank(),
       axis.line.x = element_line(),
       axis.ticks.x = element_line(), 
       axis.text.x = element_text(),
       axis.title.x = element_text(),
       axis.line.y = element_line(),
       axis.ticks.y = element_line(),
       axis.text.y = element_text(),
       axis.title.y = element_text())
x_axis_only_theme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line.x = element_line(),
  axis.ticks.x = element_line(), 
  axis.text.x = element_text(),
  axis.title.x = element_text(),
  axis.line.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y = element_blank(),
  axis.title.y = element_blank())
library(extrafont)
font_import()

loadfonts(device = "win")



perturb_data<-data.frame("x"=seq(-6,6,length.out=100),"logis"=NA,"gauss"=NA,"sum"=NA,"pertutbation"=NA)
perturb_data$logis<-dlogis(perturb_data$x)
perturb_data$gauss<-dnorm(perturb_data$x,1)
perturb_data$sum<-perturb_data$gauss+perturb_data$logis
perturb_data$pertutbation<-perturb_cont(function(x){dlogis(x)},function(x){dnorm(x,1)})(seq(-6,6,length.out=100))

perturb_plot<-ggplot(perturb_data) +
  geom_line(mapping=aes(x=x, y=logis), color="lightgray", linewidth=1) +
  geom_line(mapping=aes(x=x, y=gauss), color="gray", linewidth=1) +
  geom_line(mapping=aes(x=x, y=sum), linetype="dashed") +
  geom_line(mapping=aes(x=x, y=pertutbation)) +
  annotate("text",x=1,y=0.15,label="g",family="DejaVu Sans")+
  annotate("text",x=0.8,y=0.35,label="f",family="DejaVu Sans")+
  annotate("text",x=0.6,y=0.5,label="f \u2295 g",family="DejaVu Sans")+
  annotate("text",x=-0.6,y=0.55,label="f+g",family="DejaVu Sans")+
  axis_only_theme +
  ylab("Density") +
  ggtitle("Perturbation") +ylim(c(0,0.62))+
  theme(plot.title = element_markdown(hjust = 0.5))

perturb_plot

power_data<-data.frame("x"=seq(-6,6,length.out=100),"logis"=NA,"gauss"=NA,"prod"=NA,"powering"=NA)
power_data$logis<-dlogis(perturb_data$x)
power_data$prod<-perturb_data$logis*2
power_data$powering<-power_cont(function(x){dlogis(x)},a=2)(seq(-6,6,length.out=100))


power_plot<-ggplot(power_data) +
  geom_line(mapping=aes(x=x, y=logis), color="lightgray", linewidth=1) +
  geom_line(mapping=aes(x=x, y=prod), linetype="dashed") +
  geom_line(mapping=aes(x=x, y=powering)) +
  annotate("text",x=0,y=0.23,label="g",family="DejaVu Sans")+
  annotate("text",x=0,y=0.32,label="2 \u2299 g",family="DejaVu Sans")+
  annotate("text",x=0,y=0.45,label="2 \u22C5 g",family="DejaVu Sans")+
  x_axis_only_theme+
  ylab("Density") +
  ggtitle("Powering") + ylim(c(0,0.62))+
  theme(plot.title = element_markdown(hjust = 0.5))
power_plot  

grid.arrange(perturb_plot,power_plot,nrow=1)




f_P<-function(x){
  dnorm(x,mean=1)/dnorm(x,mean=0)
}
g_P<-function(x){
  exp(-x+x^2/2)/(1+exp(-x))^2
}

perturb_cont(f_P,g_P)
P_perturb_data<-data.frame("x"=seq(-6,6,length.out=100),"f_P"=NA,"g_P"=NA,"sum"=NA,"pertutbation"=NA)
P_perturb_data$f_P<-f_P(P_perturb_data$x)
P_perturb_data$g_P<-g_P(P_perturb_data$x)
P_perturb_data$sum<-P_perturb_data$f_P+P_perturb_data$g_P
P_perturb_data$pertutbation<-perturb_cont(f_P,g_P)(seq(-6,6,length.out=100))

P_perturb_plot<-ggplot(P_perturb_data) +
  geom_line(mapping=aes(x=x, y=f_P), color="lightgray", linewidth=1) +
  geom_line(mapping=aes(x=x, y=g_P), color="gray", linewidth=1) +
  geom_line(mapping=aes(x=x, y=sum), linetype="dashed") +
  geom_line(mapping=aes(x=x, y=pertutbation)) +
  #annotate("text",x=1,y=0.15,label="g",family="DejaVu Sans")+
  #annotate("text",x=0.8,y=0.35,label="f",family="DejaVu Sans")+
  #annotate("text",x=0.6,y=0.5,label="f \u2295 g",family="DejaVu Sans")+
  #annotate("text",x=-0.6,y=0.55,label="f+g",family="DejaVu Sans")+
  axis_only_theme +
  ylab("Density") +
  ggtitle("Perturbation") +#ylim(c(0,0.62))+
  theme(plot.title = element_markdown(hjust = 0.5))

P_perturb_plot

power_data<-data.frame("x"=seq(-6,6,length.out=100),"logis"=NA,"gauss"=NA,"prod"=NA,"powering"=NA)
power_data$logis<-dlogis(perturb_data$x)
power_data$prod<-perturb_data$logis*2
power_data$powering<-power_cont(function(x){dlogis(x)},a=2)(seq(-6,6,length.out=100))


power_plot<-ggplot(power_data) +
  geom_line(mapping=aes(x=x, y=logis), color="lightgray", linewidth=1) +
  geom_line(mapping=aes(x=x, y=prod), linetype="dashed") +
  geom_line(mapping=aes(x=x, y=powering)) +
  annotate("text",x=0,y=0.23,label="g",family="DejaVu Sans")+
  annotate("text",x=0,y=0.32,label="2 \u2299 g",family="DejaVu Sans")+
  annotate("text",x=0,y=0.45,label="2 \u22C5 g",family="DejaVu Sans")+
  x_axis_only_theme+
  ylab("Density") +
  ggtitle("Powering") + ylim(c(0,0.62))+
  theme(plot.title = element_markdown(hjust = 0.5))
power_plot  

grid.arrange(perturb_plot,power_plot,nrow=1)



#

# Install and load libraries
# install.packages("ggplot2")
# install.packages("reshape2")
library(ggplot2)
library(reshape2)

# Define the range of x-values
x_values <- seq(-2.5, 3, length.out = 20000)

# Calculate y-values for each function
f1 <- exp(-x_values)
f1<-f1/trapz(x_values,f1)
f2 <- exp(-x_values + x_values^2 / 2) / (1 + exp(-x_values))^2
f2<-f2/trapz(x_values,f2)
f3 <- exp(-2 * x_values + x_values^2 / 2) / (1 + exp(-x_values))^2
f3<-f3/trapz(x_values,f3)
f2_sq<-f2^2
f2_sq<-f2_sq/trapz(x_values,f2_sq)
f_0<-rep(1,length.out=length(x_values))
f_0<-f_0/trapz(x_values,f_0)

# Create a data frame
data <- data.frame(
  x = x_values,
  f1 = f1,
  f2 = f2,
  f3 = f3,
  f2_sq=f2_sq,
  f_0=f_0
)

# Melt the data frame
data_melted <- melt(data, id.vars = 'x', variable.name = 'Function', value.name = 'y')

# Assign colors
data_melted <- data_melted %>%
  mutate(
    Color = case_when(
      Function =="f_0" ~ 'blue',
      Function %in% c('f1', 'f2') ~ 'gray',
      TRUE ~ 'black'
    ),
    Linetype = case_when(
      Function == 'f_0' ~ 'dotted',
      Function %in% c('f2', 'f2_sq') ~ 'dashed',
      TRUE ~ 'solid'
    )
  )

# Plot using ggplot2
pplot1<-ggplot(data_melted, aes(x = x, y = y, color = Color, group = Function,linetype=Linetype)) +
  geom_line(linewidth=0.5) +
  scale_color_identity() +
  labs(title = 'Operations in Gaussian Bayes space',
       x = 'x',
       y = 'y') +
   annotate("text",x=2.5,y=2,label="f\u209A \u2295 g\u209A",family="DejaVu Sans")+
   annotate("text",x=-2.3,y=0.7,label="g\u209A",family="DejaVu Sans")+
   annotate("text",x=2.86,y=0.6,label=" f\u209A",family="DejaVu Sans")+
  annotate("text",x=-1.9,y=1.2,label="2 \u2299 g\u209A",family="DejaVu Sans")+
  axis_only_theme+theme(plot.title=element_markdown(hjust=0.5),legend.position="none")
pplot1

library(pracma)

f2_sq<-f2^2
f2_sq<-f2_sq/trapz(x_values,f2_sq)

data2 <- data.frame(
  x = x_values,
  f2=f2,
  f2_sq=f2_sq
)

# Melt the data frame
data_melted2 <- melt(data2, id.vars = 'x', variable.name = 'Function', value.name = 'y')

# Assign colors
data_melted2$Color <- with(data_melted2, ifelse(
  Function == 'f2', 'gray80',
  ifelse(Function == 'f2_sq', 'black',"black")
))

pplot2<-ggplot(data_melted2, aes(x = x, y = y, color = Color, group = Function)) +
  geom_line(linewidth=0.5) +
  scale_color_identity() +
  labs(title = 'Powering in Gaussian Bayes space',
       x = 'x',
       y = 'y') +
  annotate("text",x=-1.5,y=0.25,label="2 \u2299 g_P",family="DejaVu Sans")+
  annotate("text",x=-1.5,y=0.25,label="g_P",family="DejaVu Sans")+
  axis_only_theme+theme(plot.title=element_markdown(hjust=0.5))
pplot2

grid.arrange(pplot1,pplot2,nrow=1)

2*(1-pnorm(15010,mean=15000,sd=sqrt(45)))
