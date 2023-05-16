# Barplots for any set of group of treatments
# Grouped barplot with ggplot2
# From: https://www.r-graph-gallery.com/48-grouped-barplot-with-ggplot2/
# Date: January 30th 2019 - February 1th 2019, modified: August 6th 2019
# Last edition: August 25th 2020

# Set on your working directory
setwd("~/R/VOCS_2019/060819_VOCS")

# Create a dataset

# Reading the table for this experiment
growth<-read.table("060819_GrowthKinetics.txt", sep="\t", header=T)
head(growth)

# We will work with Spores/ml
Trichoderma<-growth$Strain
Spores<-growth$Spores
Media<-growth$Media
growth.new <- data.frame(Media, Trichoderma, Spores)

# Library
# Install if necessary
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
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

# Grouped barplot in ggplot2
# Example:
# In this case we can set the y axis scale or range because we know which is the limit in 
# the highest set (experimental set)

ggplot(growth.new, aes(fill=Media, y=Spores, x=Trichoderma)) +
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Spores/ml") +
  scale_y_continuous(limits = c(0, 1.0E+08)) +
  ylab("Spores/ml") +
  xlab("Trichoderma strain") +
  theme(text = element_text(size=16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size=16, face="italic"), # For Colletotrichum italic
        axis.title=element_text(size=16),
        axis.text.x = element_text(angle = 0, hjust = 1)) +
  labs(fill = "Media")

## How to add standard error bars to the plot?
## VISIT: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

if (!require("Rmisc")) {
  install.packages("Rmisc", dependencies = TRUE)
  library(Rmisc)
}

if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE)
  library(plyr)
}

# summarySE provides the standard deviation, standard error of the mean, 
# and a (default 95%) confidence interval
tgc <- summarySE(growth.new, measurevar="Spores", groupvars=c("Media","Trichoderma"), na.rm=TRUE)
tgc

## na.rm=TRUE argument to the  end of the argument list should clear it if theres NAs**

## Then --> 
ggplot(tgc, aes(fill=Media, y=Spores, x=Trichoderma)) +
  geom_bar(position="dodge", stat="identity") +
  ggtitle("Spores production") +
  scale_fill_manual(values=c("#99FFFF", "#FFFF99", "#FF99FF")) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1.15E+08)) +
  ylab("Spores/ml") +
  xlab("Trichoderma strain") +
  theme(text = element_text(size=16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 22),
        plot.title = element_blank(), # To hide plot title
        axis.title.y = element_text(size=22),
        axis.title.x = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1, size=24, color="black"),
        axis.text.x = element_text(angle = 0, hjust = 1, size=24, color="black")) +
  labs(fill = "Media") +
  geom_errorbar(aes(ymin=Spores-sd, ymax=Spores+sd), # dont forget to change variable here
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))

# TO SEE COLOR PALLETES AND HEXADECIMAL COLOR CODE CHART VISIT:
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

# I´m not an expert programming in R, so it costed to me one and half of the other to make
# this Script
# Trust me!! :v


# STATISTICAL ANALYSIS FOR REPORT
# Statistical Analyses of "SPORES"
cfu.aov <- aov(Spores ~ Media * Strain, data=growth)
anova(cfu.aov)

model.tables(cfu.aov, type="means", se=T)

TukeyHSD(cfu.aov, which="Media")

TukeyHSD(cfu.aov, which="Strain")

TukeyHSD(cfu.aov, which="Media:Strain")
