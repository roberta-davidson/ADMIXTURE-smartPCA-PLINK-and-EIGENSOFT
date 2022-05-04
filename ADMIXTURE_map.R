#Initialize packages
library(maps)
library(plotrix)
library(dplyr)
library(tidyverse)
library(janitor)
library(plotly)

setwd("/Users/robertadavidson/Box Sync/Robbi_PhD/Bioinf/Sequencing/Feb2022_SAM/ADMIXTURE_PONG/America_data18")

#read in admix output matrix
admix <- read.table("/Users/robertadavidson/Box Sync/Robbi_PhD/Bioinf/Sequencing/Feb2022_SAM/ADMIXTURE_PONG/America_data18/data18.2.6.Q",
                       col.names=c("K1", "K2","K3","K4","K5","K6"))
admix <- as.data.frame(admix)
#read in gps coordinates, rows in same sample order as admixture output matric
gps <- read.table("/Users/robertadavidson/Box Sync/Robbi_PhD/Bioinf/Sequencing/Feb2022_SAM/ADMIXTURE_PONG/America_data18/data18.coord", 
                  col.names=c("Lat","Lon")) 
gps <- as.data.frame(gps) 
gps$ID = as.numeric(as.character(gps$ID))

#plot map
map(xlim=c(-90,-30), ylim=c(-60,0)) #Plot maps
map.axes() #Add axes
#add points
points(gps$Lon, gps$Lat,
       cex = 1, col="red", pch=19) #To add just points
#To add admixture pie plots – here I used K = 6
for (x in 1:nrow(gps)) {
    floating.pie(gps$Lon[x], gps$Lat[x], edges=800,
      c(admix$K1[x], admix$K2[x], admix$K3[x], admix$K4[x],admix$K5[x], admix$K6[x]), 
      radius=.6, col=c("red","blue", "green", "purple", "orange", "yellow"), alphha=0.6)
      title(main="K=6", xlab="Longitude", ylab="Latitude") }
