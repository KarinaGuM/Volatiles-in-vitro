# Making a Heatmap with ggplot of volatile compounds
# We will make a simple heatmap of identified VOCS
# Samples: Trichoderma strains grown in different media

# First set in your working folder
setwd("~/R/VOCS_2019")

# If necessary wi will download and instal what we need

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

library(ggplot2)

# Reading my .csv archive from excel
# Asegurarse de estar en la misma area, carpeta de trabajo donde esta el .csv
read.csv("28ene19_Invitro.csv")
datavocs <- read.csv("28ene19_Invitro.csv")

# Making the heatmap (white - steelblue)

ggplot(datavocs, aes(Sample, Compounds)) +
  geom_tile(aes(fill = Abundance), color = "white") +
  scale_fill_gradient2(low="white", mid="navy", high="black", 
                       midpoint=3e+08, limits=range(datavocs$Abundance)) +
  theme_classic() +
  ylab("Volatile Compounds") +
  xlab("Samples") +
  theme(text = element_text(size=17),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=17),
        axis.title=element_text(size=18,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Relative Abundance")