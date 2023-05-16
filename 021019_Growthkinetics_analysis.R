# karina.gutierrez@cinvestav.mx
# Karina Gutierrez Moreno
# PhD Student on Integrative Biology (Cinvestav Irapuato Mexico)
# Supervisor: Martin Heil
# Statistical analysis of 2019 field experiments
# August 6th 2019 
# Edited: September 6th 2019

# Set on your working directory
setwd("~/R/VOCS_2019/060819_VOCS")

# Create a dataset
# Reading the table for this experiment
growth <- read.table("021019_growthrate_diameter.txt", sep="\t", header=T)
head(growth)
growth$Medium <- as.factor(growth$Medium)
summary(growth$Medium)

# Segregating the data of the different media
# Segregating each medium into different data frames.
mea <- growth[as.factor(growth$Medium)=="MEA",] # Only MEA
pda <- growth[as.factor(growth$Medium)=="PDA",] # Only PDA
sea <- growth[as.factor(growth$Medium)=="SEA",] # Only SEA

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

if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
}

if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE)
  library(plyr)
}

if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

if (!require("lme4")) {
  install.packages("lme4", dependencies = TRUE)
  library(lme4)
}

# Visit http://www.danmirman.org/gca for more details
# Working first with MEA data
summary(mea)
head(mea)

# For statistical analysis
# The first step is to create a third-order polynomial in the range of Replicate.
t <- poly((unique(mea$Replicate)), 4) 

# The next step is to create variables ot1, ot2, ot3 corresponding to the orthogonal 
# polynomial time terms and populate their values with the Replicate-appropriate orthogonal 
# polynomial values:
mea[,paste("ot", 1:4, sep="")] <- t[mea$Replicate, 1:4]

# Since this is a simple case with just one within-subjects fixed effect that has only 
# two levels, we can skip to the full model and examine its parameter estimates:
m.full <- lmer(Diameter ~ (ot1+ot2+ot3+ot4)* Strain + 
                 (ot1+ot2+ot3+ot4 | Medium) + 
                 (ot1+ot2+ot3 | Medium:Strain), 
               control = lmerControl(optimizer="bobyqa"),
               data=mea, REML=F)
coef(summary(m.full))

p <- ggplot(mea, aes(x=Time.H, y=Diameter, color=Strain)) + 
  stat_summary(fun.data=mean_se, geom="pointrange") + 
  ggtitle("Growth kinetics on MEA") +
  labs(y="Mycelium diameter (cm)", x="Time after inoculation (H)") + 
  theme_bw() + 
  theme(text=element_text(size=14),
        legend.title=element_text(size = 14),
        legend.text=element_text(size = 14),
        plot.title=element_text(size=14),
        axis.title=element_text(size=14),
        axis.text.y=element_text(angle=0, hjust=1, size=14), 
        axis.text.x=element_text(angle=0, hjust=1, size=14)) +
  scale_color_manual(values=c("darkorchid4","deepskyblue3","forestgreen","orange"))

p 

last_plot() + stat_summary(aes(y=fitted(m.full)), fun.y=mean, geom="line")


