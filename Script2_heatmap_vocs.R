# Making a Heatmap with ggplot of volatile compounds
# We will make a simple heatmap of identified VOCS
# Samples: Trichoderma strains grown in different media
# You can use this Script for any volatile compounds identified in any samples!

# Version 2.0 (Karina Rules!)
# For more info go to: http://rpubs.com/htejero/212365

# First set in your working folder
setwd("~/R/VOCS_2019")

#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("gridExtra")) {
  install.packages("gridExtra", dependencies = TRUE)
  library(gridExtra)
}

#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

data <- read.csv("230519_Log2_invitromatrix_NA.csv", comment.char="#")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
vocs <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(vocs) <- rnames                  # assign row names

#########################################################
### C) Customizing and plotting the heat map
#########################################################

if (!require("reshape")) {
  install.packages("reshape", dependencies = TRUE)
  library(reshape)
}

mvocs <- melt(vocs)
head(mvocs)

#### PLOTTING ###
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

if (!require("forcats")) {
  install.packages("forcats", dependencies = TRUE)
  library(forcats)
}

ggplot(data = data.frame(mvocs), aes(x = fct_inorder(X2), y = fct_rev(fct_inorder(X1)), fill= value)) + 
  geom_tile() +
  ylab("Volatile Compounds") +
  xlab("Samples") +
  theme(text = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "SCALE: LOG2 RATIO SIGNAL")


## Changing colors
ggplot(data = data.frame(mvocs), aes(x = fct_inorder(X2), y = fct_rev(fct_inorder(X1)), fill= value)) + 
  geom_tile() +
  scale_fill_gradient2() +
    ylab("Volatile Compounds") +
  xlab("Samples") +
  theme(text = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "LOG2 RATIO SIGNAL")

#Option 3 <----------------
ggplot(data = data.frame(mvocs), aes(x = fct_inorder(X2), y = fct_rev(fct_inorder(X1)), fill= value)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Paired", direction = 1) +
  ylab("Volatile Compounds") +  xlab("Samples") +
  theme(text = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Abundance: Peak Area")

# OPTION 4 # I like this for in vitro cuz has low peak areas
ggplot(data = data.frame(mvocs), aes(x = fct_inorder(X2), y = fct_rev(fct_inorder(X1)), fill= value)) + 
  geom_tile() +
  scale_fill_gradient2(low="white", mid="blue", high="red", 
                       midpoint=7e+07) +
  ylab("Volatile Compounds") +  xlab("Samples") +
  theme(text = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=14),
        axis.title=element_text(size=14,face="bold"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Scale")


# Let's add broders to the square in OPTION 4 <---- I LIKE THIS =)
ggplot(data = data.frame(mvocs), aes(x = fct_inorder(X2), y = fct_rev(fct_inorder(X1)), fill= value)) + 
  geom_tile(colour = "gray") +
  scale_fill_gradient2(low="green", mid="black", high="red", 
                       midpoint=0) +
  ylab("Volatile Compounds") +  xlab("Samples") +
  theme(text = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=14),
        axis.title=element_text(size=14,face="bold"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Scale")


# OPTION 5 THIS ONE FOR COLORPALETTE #030519 with borders inside 230519
ggplot(data = data.frame(mvocs), aes(x = fct_inorder(X2), y = fct_rev(fct_inorder(X1)), fill= value)) + 
  geom_tile(colour = "gray") +
  scale_fill_distiller(palette = "PuRd", direction = 1) +
  ylab("Volatile Compounds") +  xlab("Samples") +
  theme(text = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=14),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Scale for Peak Area")


# Changing Labels and Scale format
# https://www.littlemissdata.com/blog/heatmaps

ggplot(data = data.frame(mvocs), aes(x = fct_inorder(X2), y = fct_rev(fct_inorder(X1)), fill= value), colour = "white", na.rm = TRUE) + 
  geom_tile(colour = "gray") +
  scale_fill_gradient2(low="green", mid="black", high="red", 
                       midpoint=0) +
  guides(fill=guide_legend(title="Scale")) +
  ylab("Volatile Compounds") +  xlab("Samples") +
  theme(text = element_text(size=16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=14),
        axis.title=element_text(size=14,face="bold"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="",y="", title = "Identified volatile compounds")


# Another helpfull pages: 
# https://stackoverflow.com/questions/39361430/r-program-how-to-avoid-ggplot-re-order-in-both-x-axis-and-y-axis
# http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually#default-colors
