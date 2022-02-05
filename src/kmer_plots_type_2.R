#####################################################
# Name: kmer_plots.R
# Description: Generates stacked bar charts using
#              results produced by the snakemake
#              workflow.
#
# Date: January 14, 2022
#
#####################################################

library(ggplot2)
library(data.table)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
########################################################################
within_groups_data_path <- "/Users/omarahmed/Downloads/current_research/khoice_exps/type_2/exp_1a/within_dataset_analysis.csv"
across_groups_data_path <- "/Users/omarahmed/Downloads/current_research/khoice_exps/type_2/exp_1a/across_dataset_analysis.csv"

dataset_names <- c("Bacillus cereus", "Bacillus anthracis", "Bacillus thuringiensis", "Bacillus weihenstephanensis")
working_dir <- "/Users/omarahmed/Downloads/current_research/khoice_exps/type_2/exp_1a/"

########################################################################
# Methods for generating the two types of plots: within groups, and 
# across of multiple groups
########################################################################

make_within_groups_plot <- function(group_df, name) {
  # Convert to data.table, melt, and then revert to data.frame
  setDT(group_df)
  group_df_melt <- melt(group_df, 
                        measure.vars = c("percent_75_or_more","percent_25_to_75", "percent_25_or_less", "percent_1_occ"),
                        variable.name = "range", 
                        value.name = "percent")
  setDF(group_df_melt)
  
  # Build the stacked bar chart ...
  graph_title <- paste("Uniqueness of Kmers Across", name, "Genomes w.r.t a Pivot Genome")
  plot <- ggplot(group_df_melt, aes(fill=range, x=factor(k), y=percent)) + 
          geom_bar(position="fill", stat="identity") +
          theme_classic() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =14),
                axis.title.y=element_text(size=14),
                legend.position = "bottom", 
                legend.text=element_text(size=10),
                legend.box="vertical",
                legend.title=element_text(size=10),
                axis.text=element_text(size=14, color="black")) +
        labs(x="Kmer Length (k)",
             y="Ratio of Unique Kmers",
             title=graph_title) +
        scale_fill_discrete(name = "Number of Genomes the Kmers Occur In:",
                            labels = c(">75% of Genomes", "25% to 75% of Genomes","<25% of Genomes", "Only in Pivot")) +
        guides(fill= guide_legend(nrow = 2))
  return (plot)
}

make_uniqueness_plot <- function(group_df) {
  graph_title <- "Uniqueness statistic as k increases for each dataset"
  plot <- ggplot(group_df, aes(x=k, y=unique_stat)) + 
          geom_line(aes(color=group_num)) +
          geom_point(aes(color=group_num)) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =14),
                axis.title.y=element_text(size=14),
                legend.position = "bottom", 
                legend.text=element_text(size=10),
                legend.box="vertical",
                legend.title=element_text(size=10),
                axis.text=element_text(size=14, color="black")) +
          labs(x="Kmer Length (k)",
               y="Uniqueness Statistic",
               title=graph_title) +
          scale_color_discrete(name="Dataset",
                               labels=dataset_names) +
          #scale_shape_discrete(name="Dataset",labels=dataset_names) +
          scale_x_continuous(breaks=seq(0, 100, 5)) +
          geom_hline(yintercept=1, linetype="dashed", color = "red", size=1) +
          guides(fill= guide_legend(nrow = 3))
          
  return (plot)
}

make_across_groups_plot <- function(group_df, name) {
  # Convert to data.table, melt, and then revert to data.frame
  setDT(group_df)
  group_df_melt <- melt(group_df, 
                        measure.vars = c("percent_4_to_8","percent_2_to_3", "percent_1_occ"),
                        variable.name = "range", 
                        value.name = "percent")
  setDF(group_df_melt)
  
  # Build the stacked bar chart ...
  graph_title <- paste("Uniqueness of Kmers Across All", length(dataset_names), "Groups w.r.t", name, "Pivot")
  plot <- ggplot(group_df_melt, aes(fill=range, x=factor(k), y=percent)) + 
          geom_bar(position="fill", stat="identity") +
          theme_classic() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =14),
                axis.title.y=element_text(size=14),
                legend.position = "bottom", 
                legend.text=element_text(size=10),
                legend.box="vertical",
                legend.title=element_text(size=10),
                axis.text=element_text(size=14, color="black")) +
          labs(x="Kmer Length (k)",
               y="Ratio of Unique Kmers",
               title=graph_title) +
          scale_fill_discrete(name = "Number of Groups the Kmers Occur In:",
                        labels = c("4 Groups", "2 to 3 Groups", "Only Pivot"))
  return (plot)
}

make_uniqueness_plot_across_groups <- function(group_df) {
  graph_title <- paste("Uniqueness statistic across all", length(dataset_names),"datasets as k increases")
  plot <- ggplot(group_df, aes(x=k, y=unique_stat, type=group_num, color=group_num)) + 
          geom_line() +
          geom_point() +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =14),
                axis.title.y=element_text(size=14),
                legend.position = "bottom", 
                legend.text=element_text(size=10),
                legend.box="vertical",
                legend.title=element_text(size=10),
                axis.text=element_text(size=14, color="black")) +
          labs(x="Kmer Length (k)",
               y="Uniqueness Statistic",
               title=graph_title) +
          scale_x_continuous(breaks=seq(0, 100, 5)) +
          geom_hline(yintercept=1, linetype="dashed", color = "red", size=1) +
          scale_color_discrete(name="Pivot Group", labels=dataset_names)
  
  return (plot)
}

#########################################################################
# Start of the "main" method of code ...
########################################################################

## SECTION 1 - Looking at data regarding within groups
within_df <- read.csv(within_groups_data_path, header=TRUE)

# Creates the individual plots for each group of data ...
pos <- 1
for (dataset_num in unique(within_df$group_num)) { 
  within_df_subset <- subset(within_df, within_df$group_num == dataset_num)
  curr_plot <- make_within_groups_plot(within_df_subset, dataset_names[pos])
  
  output_name <- paste(working_dir, "within_dataset_", pos, "_kmer_analysis.png", sep="")
  ggsave(output_name, plot=curr_plot, dpi=800, device="jpeg", width=9.5, height=7)
  pos <- pos + 1
}

# Creates the plot looking at the uniqueness measure
curr_plot <- make_uniqueness_plot(within_df)
output_name <- paste(working_dir, "within_dataset_unique_stat.png", sep="")
ggsave(output_name, plot=curr_plot, dpi=800, device="jpeg", width=9, height=6)


## SECTION 2 - Looking at data regarding across groups
across_df <- read.csv(across_groups_data_path, header=TRUE)

# Creates the individual plots for each pivot against all other groups
pos <- 1
for (dataset_num in unique(across_df$group_num)) {
  across_df_subset <- subset(across_df, across_df$group_num == dataset_num)
  curr_plot <- make_across_groups_plot(across_df_subset, dataset_names[pos])
  
  output_name <- paste(working_dir, "across_datasets_", pos, "_kmer_analysis.png", sep="")
  ggsave(output_name, plot=curr_plot, dpi=800, device="jpeg", width=9.5, height=7)
  pos <- pos + 1
}

# Next, make a plot for the uniqueness statistic across groups ...
curr_plot <- make_uniqueness_plot_across_groups(across_df)
output_name <- paste(working_dir, "across_dataset_unique_stat.png", sep="")
ggsave(output_name, plot=curr_plot, dpi=800, device="jpeg", width=9, height=6)
