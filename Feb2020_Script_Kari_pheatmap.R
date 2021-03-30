# Script by Frank Medina and Kamil Slowkowski 
# Edited by Karina G. Moreno
# Cinvestav Unidad Irapuato
# February 2020
# LAST EDITION: AUGUST 20TH 2020

library(gplots)
library(pheatmap)
library(RColorBrewer)
library(viridis) 

setwd("~/R/Heatmaps SPME VOCs 2020")

# Reading the data
ext_todos <- read.csv("100320_MIX_vocs.csv")

# Converting to matrix
ext_matrix <- data.matrix(ext_todos)
rownames(ext_matrix)<- ext_todos[,1] # Asigning row names (compounds or genes)
ext_matrix[is.na(ext_matrix)] <- 0 # Converting NA or Blanks to Zero

# Making the Heatmap with pheatmap
pheatmap(ext_matrix[,2:26], main="Petri dish",fontsize_row= 10, show_rownames = T, display_numbers = F, fontsize = c(12), number_color = "black",number_format = "%.0f", color = greenred(75))

# Changing colors and clusters
pheatmap(ext_matrix[,2:13], main="Petri dish",fontsize_row= 11, show_rownames = T, display_numbers = F, cluster_rows = TRUE, cluster_cols = FALSE, fontsize = c(10), number_color = "black", number_format = "%.0f", clustering_method = "average", color = cm.colors(5)) # Aqui el numero indica que tan degradado estara la escala
pheatmap(ext_matrix[,2:13], main="Petri dish",fontsize_row= 11, 
         show_rownames = T, display_numbers = F, cluster_rows = TRUE, 
         cluster_cols = FALSE, fontsize = c(10), number_color = "black", 
         number_format = "%.0f", clustering_method = "average", 
         color = colorRampPalette(brewer.pal(7, "YlGnBu"))(150))

# Pallets of colors for RColorBrewer
# https://www.r-graph-gallery.com/38-rcolorbrewers-palettes.html

pheatmap(ext_matrix[,2:26], main="MIX",fontsize_row= 11, 
         show_rownames = T, display_numbers = F, cluster_rows = TRUE, 
         cluster_cols = FALSE, fontsize = c(10), number_color = "black", 
         number_format = "%.0f", color = cividis(50))


pheatmap(ext_matrix[,2:13], main="Petri dish",fontsize_row= 11, 
         show_rownames = T, 
         display_numbers = F, 
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         clustering_distance_rows = "correlation",#Pearson's
         clustering_method = "average",
         gaps_col=c(4,8),
         fontsize = c(11), 
         color = colorRampPalette(brewer.pal(7, "YlGnBu"))(150))

pheatmap(ext_matrix[,2:13], main="Petri dish",fontsize_row= 11, 
         show_rownames = T, display_numbers = F, cluster_rows = TRUE, 
         cluster_cols = FALSE, fontsize = c(10), number_color = "black",
         clustering_distance_rows = "correlation",#Pearson's
         clustering_method = "average",
         gaps_col=c(4,8),
         color = colorRampPalette(brewer.pal(9, "Spectral"))(150))

pheatmap(ext_matrix[,2:26], main="Sterile Soil",fontsize_row= 12, 
         show_rownames = T, 
         display_numbers = F, 
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         clustering_distance_rows = "correlation",#Pearson's
         clustering_method = "average",
         gaps_col=c(5,10,15,20),
         angle_col= 45,
         fontsize = c(14), 
         color = viridis(250)) # Numero menor para Mix

pheatmap(ext_matrix[,2:13], main="Petri dish",fontsize_row= 11, 
         show_rownames = T, 
         display_numbers = F, 
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         clustering_distance_rows = "correlation", #Pearson's
         clustering_method = "average",
         gaps_col=c(4,8), #change for your gap_col by treatments (e.j.every 5)
         fontsize = c(11), 
         color = viridis(350, direction = -1)) # change for any viridis colorpallette

# https://slowkow.com/notes/pheatmap-tutorial/#uniform-breaks

library(ggplot2)

# Set the theme for all the following plots.
theme_set(theme_bw(base_size = 16))

dat <- data.frame(values = as.numeric(ext_matrix))
ggplot(dat, aes(values)) + 
  geom_density(bw = "SJ", fill = "salmon", colour= NA) +
  geom_rug(colour = "salmon") +
  labs(x = "Abundances", y = "Density") +
  theme_minimal()


# Uniform breaks
# We can visualize the unequal proportions of data represented by each color:
mat_breaks <- seq(min(ext_matrix), max(ext_matrix), length.out = 10)


## ----uniform-color-breaks-detail, fig.height=2, echo=FALSE---------------
dat_colors <- data.frame(
  xmin = mat_breaks[1:(length(mat_breaks)-1)],
  xmax = mat_breaks[2:length(mat_breaks)],
  ymin = 0,
  ymax = max(density(ext_matrix, bw = "SJ")$y),
  fill = rev(inferno(length(mat_breaks) - 1)), # you can change for any viridis
  stringsAsFactors = FALSE
)
ggplot() +
  geom_rect(
    data = dat_colors,
    mapping = aes(
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill
    )
  ) +
  geom_density(
    data = dat,
    mapping = aes(values),
    bw = "SJ", color = "cyan"
  ) +
  scale_fill_manual(values = dat_colors$fill) +
  theme(legend.position = "none") +
  labs(title = "Uniform breaks")

## ----uniform-color-breaks-bars, fig.height=3, echo=FALSE-----------------
dat2 <- as.data.frame(table(cut(
  ext_matrix, mat_breaks
)))
dat2$fill <- inferno(nrow(dat2)) #change for any viridis colorpallette
ggplot() +
  geom_bar(
    data = dat2,
    mapping = aes(x = Var1, weight = Freq, fill = Var1),
    color = "black", size = 0.1
  ) +
  coord_flip() +
  scale_fill_manual(values = dat2$fill) +
  theme(legend.position = "none") +
  labs(y = "data points", x = "breaks",
       title = "Number of data points per color")


# With our uniform breaks and non-uniformly distributed data, we represent 
# 86.5% of the data with a single color.
# On the other hand, 6 data points greater than or equal to 100 are represented 
# with 4 different colors.


# Quantile breaks
# If we reposition the breaks at the quantiles of the data, then each color will 
# represent an equal proportion of the data:
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(ext_matrix, n = 11)


#### Taken from pheatmap-tutorial.R of Kamil Slowikowski: pheatmap R 
# package by Raivo Kolde
# https://slowkow.com/notes/pheatmap-tutorial/

## ----quantile-color-breaks-detail, fig.height=2, echo=FALSE--------------
dat_colors <- data.frame(
  xmin = mat_breaks[1:(length(mat_breaks)-1)],
  xmax = mat_breaks[2:length(mat_breaks)],
  ymin = 0,
  ymax = max(density(ext_matrix, bw = "SJ")$y),
  fill = rev(viridis(length(mat_breaks) - 1)),
  stringsAsFactors = FALSE
)
ggplot() +
  geom_rect(
    data = dat_colors,
    mapping = aes(
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill
    )
  ) +
  geom_density(
    data = dat,
    mapping = aes(values),
    bw = "SJ", color = "cyan"
  ) +
  scale_fill_manual(values = dat_colors$fill) +
  theme(legend.position = "none") +
  labs(title = "Quantile breaks")

## ----quantile-color-breaks-bars, fig.height=3, echo=FALSE----------------
dat2 <- as.data.frame(table(cut(
  ext_matrix, mat_breaks
)))
dat2$fill <- viridis(nrow(dat2))
ggplot() +
  geom_bar(
    data = dat2,
    mapping = aes(x = Var1, weight = Freq, fill = Var1),
    color = "black", size = 0.1
  ) +
  coord_flip() +
  scale_fill_manual(values = dat2$fill) +
  theme(legend.position = "none") + # Edited 20/08/20
  theme(axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        axis.title=element_text(size=22,face="bold"),
        plot.title = element_text(size=22)) +
  labs(y = "data points", x = "breaks",
       title = "Number of data points per color")


############################################################################

# Then:
#' When we use quantile breaks in the heatmap, we can clearly see that
#' group 1 values are much larger than values in groups 2 and 3, and we can
#' also distinguish different values within groups 2 and 3:
#' 
#' 
pheatmap(ext_matrix[,2:13], main="Petri dish",fontsize_row= 11, 
         show_rownames = T, display_numbers = F, cluster_rows = TRUE, 
         cluster_cols = FALSE, fontsize = c(10), number_color = "black", 
         number_format = "%.0f", color = plasma(length(mat_breaks)-1), breaks= mat_breaks)

pheatmap(ext_matrix[,2:13], main="Petri dish",fontsize_row= 11, 
         show_rownames = T, display_numbers = F, row_dend_width = unit(4, "cm"), 
         cluster_cols = FALSE, fontsize = c(10), number_color = "black", 
         number_format = "%.0f", color = magma(length(mat_breaks)-1), breaks= mat_breaks)


#########################################################################

# FINAL PHEATMAP SCRIPT
# DON´T FORGET TO MODIFY WITH GAPS_COL BY TREATMENTS
# CHANGE ALSO ,2:NUMBER OF TOTAL COLUMNS IN YOUR DATA
pheatmap(ext_matrix[,2:26], main="Sterile soil",fontsize_row= 11, 
         show_rownames = T, 
         display_numbers = F, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         fontsize = c(12), 
         clustering_distance_rows = "correlation",#Pearson's
         clustering_method = "average",
         gaps_col=c(5,10,15,20),
         angle_col= 45,
         color = viridis(length(mat_breaks)-1), 
         breaks= mat_breaks)