#####################################################
# Name: kmer_plots_type_5.R
# Description: Creates a confusion matrix heat map for
#              Experiment 5
#
# Date: June 28, 2022
#
#####################################################

library(ggplot2)
library(reshape)
library(tidyr)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################
confusion_matrix_data_path <- "/Users/mwche/Downloads/exp5_salmonella_ecoli/confusion_matrix.csv" # FILL IN

#dataset_names <- c("B. Cereus", "B. Anthracis", "B. Thuringiensis", "B. Weihenstephanensis ")
dataset_names <- c("E. Coli", "Salmonella")

working_dir <- "/Users/mwche/Downloads/exp5_salmonella_ecoli/"

########################################################################
# Methods for generating accuracy plot
########################################################################

make_confusion_heatmap <- function(data_melt,k){
  colors = c("brown3","darkgoldenrod1","chartreuse2")
  plot <- ggplot(data_melt, aes(Predicted, True)) +                           
    geom_tile(aes(fill = Value)) + 
    geom_text(aes(label = Value), color = "white", size = 4) +
    scale_fill_gradientn(colors = colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust=1),
          plot.title=element_text(hjust = 0.5, size=14, face="bold"),
          axis.title.x=element_text(size =14),
          axis.title.y=element_text(size=14),
          legend.title=element_text(size=12),
          axis.text=element_text(size=12, color="black")) +
    coord_fixed() + 
    guides(fill = guide_colourbar(title = "% Kmers")) 
  return (plot)
}

#########################################################################
# Start of the main method of code ...
#########################################################################

# Generate Confusion Heat map for k
confusion <- read.csv(confusion_matrix_data_path, header=FALSE)

confusion_matrix = as.matrix(confusion)

# Normalize confusion matrix values to # hm
group_sums = rowSums(confusion_matrix)
norm_matrix <-matrix(nrow=nrow(confusion_matrix),ncol=ncol(confusion_matrix))
for(i in 1:nrow(norm_matrix)){
  for(j in 1:ncol(norm_matrix)){
    norm_matrix[i,j] = round(confusion_matrix[i,j]/group_sums[i], 2)
  }
}

# Modify row and col names
rownames(norm_matrix) <- dataset_names
colnames(norm_matrix) <- dataset_names


# Transform normalized matrix
data_melt <- melt(norm_matrix)
colnames(data_melt) <- c("True","Predicted","Value")
my_plot <- make_confusion_heatmap(data_melt)   
print(my_plot)
output_name <- paste(working_dir, "confusion_matrix_normalized.png", sep="")
print("output_name")
ggsave(output_name, plot=my_plot, dpi=800, device="jpeg", width=7, height=6)







