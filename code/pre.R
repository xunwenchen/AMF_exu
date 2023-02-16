# preparation ----

# ~~~~ loading packages and own functions ----

my_packages <- c('readxl', 'phyloseq', 'ape', 'Biostrings', 'tidyverse', 'dplyr', 'pairwiseAdonis', 'reshape2', 'ggplot2', 'ggpubr', 'DESeq2', 'GGally', 'igraph', 'vegan', 'ecodist', 'agricolae', 'stats')

lapply(my_packages, library, character.only = TRUE) # load my_packages

theme_set(theme_bw()) # Set plot theme of ggplot2

source('code/fun.R') # load own function(s)

# ~~~~ loading data   ----
