#####################################################
# Name: kmer_plots_type_3.R
# Description: Generates grouped bar charts for the
#              type 3 experiment from the Snakemake
#              workflow.
#
# Date: February 10, 2022
#####################################################

library(ggplot2)
library(data.table)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################
data_path <- "/Users/omarahmed/Downloads/final_analysis_type3.csv"

dataset_names <- c("Escherichia coli","Salmonella enterica")
seq_names <- c("Illumina", "ONT")
working_dir <- "/Users/omarahmed/downloads/current_research/khoice_exps/type_3/exp_1/"

########################################################################
# Methods for generating the grouped bar chart
########################################################################

make_group_bar_plot <- function(input_df, read_type, pivot_name) {
  # Create the grouped bar-chart
  input_df$dataset_num <- as.factor(input_df$dataset_num)
  input_df$k <- as.factor(input_df$k)
  
  graph_title <- paste("Kmer Overlap % of", read_type, "Reads Simulated from", pivot_name, "Varying k")
  plot <- ggplot(input_df, aes(fill=dataset_num, x=k, y=intersection_percent)) + 
    geom_bar(position="dodge", stat="identity")+
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
    labs(x="Kmer Length (k)",
         y="% of Pivot's Kmers Present in Dataset",
         title=graph_title) +
    scale_fill_discrete(name="Species", labels = dataset_names)
  return (plot)
}

#########################################################################
# Start of the main method of code ...
########################################################################

all_data_df <- read.csv(data_path, header=TRUE)

seq_pos <- 1
for (seq_type in unique(all_data_df$read_type)) { 
  pivot_pos <- 1
  for (pivot_id in unique(all_data_df$pivot_num)){

    curr_df <- subset(all_data_df, read_type==seq_type & pivot_num == pivot_id)
    read_type <- seq_names[seq_pos]
    pivot_name <- dataset_names[pivot_pos]
    curr_plot <- make_group_bar_plot(curr_df, read_type, pivot_name)
    
    output_name <- paste(working_dir, "kmer_overlap_", read_type, "_", pivot_pos, "_kmer_analysis.png", sep="")
    ggsave(output_name, plot=curr_plot, dpi=800, device="jpeg", width=9.5, height=7)
    pivot_pos <- pivot_pos + 1
  }
  seq_pos <- seq_pos + 1
}
