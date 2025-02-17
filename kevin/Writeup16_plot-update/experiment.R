# from https://stackoverflow.com/questions/34810857/plotting-a-kde-result-in-ggtern
#Required Libraries
library(ggtern)
library(ggplot2)
library(compositions)
library(MASS)
library(scales)

set.seed(1) #For Reproduceability
mydata <- data.frame(
  x = runif(100, min = 0.25, max = 0.5),
  y = runif(100, min = 0.1, max = 0.4),
  z = runif(100, min = 0.5, max = 0.7)) 

#VARIABLES
nlevels  = 7
npoints  = 200
expand   = 0.5

#Prepare the data, put on isometric logratio basis
df     = data.frame(acomp(mydata)); colnames(df) = colnames(mydata)
data   = data.frame(ilr(df)); colnames(data) = c('x','y')

#Prepare the Density Estimate Data
h.est  = c(MASS::bandwidth.nrd(data$x), MASS::bandwidth.nrd(data$y))
lims   = c(expand_range(range(data$x),expand),expand_range(range(data$y),expand))
dens   = MASS::kde2d(data$x,data$y,h=h.est,n=npoints,lims=lims)

#-------------------------------------------------------------
#<<<<< Presumably OP has data at this point, 
#      and so the following should achieve solution
#-------------------------------------------------------------

#Generate the contours via ggplot2's non-exported function
lines  = ggplot2:::contour_lines(data.frame(expand.grid(x = dens$x, y = dens$y),
                                            z=as.vector(dens$z),group=1),
                                 breaks=pretty(dens$z,n=nlevels))

#Transform back to ternary space
lines[,names(mydata)] = data.frame(ilrInv(lines[,names(data)]))

#Render the plot
ggtern(data=lines,aes(x,y,z)) +
  theme_dark() + 
  theme_legend_position('topleft') + 
  geom_polygon(aes(group=group,fill=level),colour='grey50') +
  scale_fill_gradient(low='green',high='red') + 
  labs(fill  = "Density",
       title = "Example Manual Contours from Density Estimate Data")