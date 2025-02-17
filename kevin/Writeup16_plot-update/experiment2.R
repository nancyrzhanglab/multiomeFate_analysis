# from https://rpubs.com/KDVdecisions/triadtutorial1
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

#Notice that you could leave out the points, if you prefer
ggtern(data=mydata, aes(x=x, y=y, z=z)) +      #define data sources
  geom_density_tern()                                 #define a data geometery

#Or you can apply a color gradient to space between the contour lines
ggtern(data=mydata, aes(x=x, y=y, z=z)) +                          #define data sources
  stat_density_tern(aes(fill=..level.., alpha=..level..),geom='polygon') +#now you need to use stat_density_tern
  scale_fill_gradient2(high = "red") +                                    #define the fill color
  guides(color = "none", fill = "none", alpha = "none")                   #we don't want to display legend items
