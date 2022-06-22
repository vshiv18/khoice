#####################################################
# Name: kmer_plots_type_4.R
# Description: Plots accuracy across values of k for 
# a given set of groups
#
# Date: June 17, 2022
#
#####################################################

library(ggplot2)
library(data.table)
library(tidyr)


########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################
accuracy_score_data_path <- "/Users/mwche/Downloads/exp4_salmonella_ecoli/accuracy_scores.csv"
confusion_matrix_data_path <- "/Users/mwche/Downloads/exp4_salmonella_ecoli/k_31_confusion_matrix.csv"

#dataset_names <- c("B. Cereus", "B. Anthracis", "B. Thuringiensis", "B. Weihenstephanensis ")
dataset_names <- c("E. Coli", "Salmonella")
k = "31"
working_dir <- "/Users/mwche/Downloads/exp4_salmonella_ecoli/"

########################################################################
# Methods for generating accuracy plot
########################################################################
make_accuracy_plot <- function(accuracy_df, dataset_names){
  plot <- ggplot(data=accuracy_df, aes(x=k, y= accuracy_score, color = dataset)) +
    geom_line()+
    geom_point() + 
    theme_bw() +
    theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
          axis.title.x=element_text(size =14),
          axis.title.y=element_text(size=14),
          legend.position = "bottom", 
          legend.text=element_text(size=12),
          legend.box="vertical",
          legend.title=element_text(size=12),
          axis.text=element_text(size=12, color="black")) +
    scale_y_continuous(breaks=seq(0, 1.0, 0.1)) +
    scale_x_continuous(breaks = c(10, 15, 20, 25, 30, 35)) +
    labs(x="Kmer Length (k)",
         y="Accuracy",
         title="Kmer Classification Accuracy Using LCA Across Groups") +
    scale_color_discrete(name = "", labels = dataset_names)
  return (plot)
}

make_confusion_heatmap <- function(data_melt,k){
  graph_title <- paste("Confusion Matrix for k =",k)
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
    labs(title= graph_title) +
    coord_fixed() + 
    guides(fill = guide_colourbar(title = "% Kmers")) 
  return (plot)
}
  

#########################################################################
# Start of the main method of code ...
########################################################################

# Generate multiline plot
accuracy_df <- read.csv(accuracy_score_data_path, header = FALSE)

# Wrangle df to generate a single accuracy column
name_id = 2
colnames(accuracy_df)[1] <- "k"
for (name in dataset_names){
  colnames(accuracy_df)[name_id] <- paste("Dataset: ", name, sep=" ")
  name_id = name_id + 1
}
long_accuracy_df <- pivot_longer(accuracy_df, cols = contains("Dataset: "), 
                                 names_to ="dataset",
                                values_to = "accuracy_score")
my_plot <- make_accuracy_plot(long_accuracy_df, dataset_names)
output_name <- paste(working_dir, "accuracy_plot.png", sep="")
ggsave(output_name, plot=my_plot, dpi=800, device="jpeg", width=7, height=6)


# Generate Confusion Heatmap for k
confusion <- read.csv(confusion_matrix_data_path, header=FALSE)

confusion_matrix = as.matrix(confusion)

group_sums = rowSums(confusion_matrix)
norm_matrix <-matrix(nrow=nrow(confusion_matrix),ncol=ncol(confusion_matrix))
for(i in 1:nrow(norm_matrix)){
  for(j in 1:ncol(norm_matrix)){
    norm_matrix[i,j] = round(confusion_matrix[i,j]/group_sums[i], 2)
  }
}

rownames(norm_matrix) <- dataset_names

cols <- append(dataset_names,"Unidentified")

colnames(norm_matrix) <- cols

data_melt <- melt(norm_matrix)
colnames(data_melt) <- c("True","Predicted","Value")
my_plot <- make_confusion_heatmap(data_melt,k)   
print(my_plot)
output_name <- paste(working_dir, "k_",k,"_confusion_matrix_normalized.png", sep="")
ggsave(output_name, plot=my_plot, dpi=800, device="jpeg", width=7, height=6)

