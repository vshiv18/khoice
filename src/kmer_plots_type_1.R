#####################################################
# Name: kmer_plots_type1.R
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
within_groups_data_path <- "/Users/omarahmed/Downloads/current_research/khoice_exps/type_1/test_1/within_datasets_analysis.csv"
across_groups_data_path <- "/Users/omarahmed/Downloads/current_research/khoice_exps/type_1/test_1/across_datasets_analysis.csv"

dataset_names <- c("Escherichia coli", "Salmonella enterica")
working_dir <- "/Users/omarahmed/Downloads/current_research/khoice_exps/type_1/test_1/"

########################################################################
# Methods for generating the two types of plots: within groups, and 
# across of multiple groups
########################################################################

make_within_groups_plot <- function(group_df, name) {
  # Convert to data.table, melt, and then revert to data.frame
  setDT(group_df)
  group_df_melt <- melt(group_df, 
                        measure.vars = c("percent_75_or_more","percent_25_to_75","percent_25_or_less","percent_1_occ"),
                        variable.name = "range", 
                        value.name = "percent")
  setDF(group_df_melt)
  
  # Build the stacked bar chart ...
  graph_title <- paste("Uniqueness of Kmers Across", name, "Genomes")
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
                            labels = c(">75% of Genomes","25% to 75% of Genomes","< 25% of Genomes","One Genome")) +
        guides(fill= guide_legend(nrow = 2))
  return (plot)
}

make_uniqueness_plot <- function(group_df) {
  # Determine the coefficient used for scaling second axis
  coeff <- max(group_df$delta_frac)/max(group_df$unique_stat)
  
  # Make the plot
  graph_title <- "Uniqueness statistic as k increases for each dataset"
  plot <- ggplot(group_df, aes(x=k)) + 
          geom_line(aes(y=unique_stat, color=group_num)) +
          geom_point(aes(y=unique_stat, color=group_num)) +
          geom_line(aes(y=delta_frac/coeff, color=group_num), linetype="dashed")+
          geom_point(aes(y=delta_frac/coeff, color=group_num))+
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
          geom_hline(yintercept=1, color = "red", size=0.5) +
          scale_y_continuous(name="Uniqueness Statistic", 
                       sec.axis = sec_axis(~.*coeff, name="Kmer Cardinality/k", 
                                           labels = function(x) format(x, scientific = TRUE)))+
          guides(fill= guide_legend(nrow = 3))
          
  return (plot)
}

make_across_groups_plot <- function(group_df) {
  # Convert to data.table, melt, and then revert to data.frame
  setDT(group_df)
  group_df_melt <- melt(group_df, 
                        measure.vars = c("percent_5_to_20","percent_2_to_5","percent_1_occ"),
                        variable.name = "range", 
                        value.name = "percent")
  setDF(group_df_melt)
  
  # Build the stacked bar chart ...
  graph_title <- paste("Uniqueness of Kmers Across the", length(dataset_names), "Groups of Genomes")
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
                        labels = c("5 to 8 Groups","2 to 4 Groups","One Group"))
  return (plot)
}

make_uniqueness_plot_across_groups <- function(group_df) {
  # Determine the coefficient used for scaling second axis
  coeff <- max(group_df$delta_frac)/max(group_df$unique_stat)
  
  # Make the plot
  graph_title <- paste("Uniqueness statistic across all", length(dataset_names),"datasets as k increases")
  plot <- ggplot(group_df, aes(x=k)) + 
          geom_line(aes(y=unique_stat), color="dodgerblue") +
          geom_line(aes(y=delta_frac/coeff), color="dodgerblue", linetype="dashed")+
          geom_point(aes(y=unique_stat), color="dodgerblue") +
          geom_point(aes(y=delta_frac/coeff), color="dodgerblue") +
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
          scale_y_continuous(name="Uniqueness Statistic", 
                             sec.axis = sec_axis(~.*coeff, name="Kmer Cardinality/k", 
                                                 labels = function(x) format(x, scientific = TRUE)))+
          geom_hline(yintercept=1, color = "red", size=0.5)
  
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
  ggsave(output_name, plot=curr_plot, dpi=800, device="jpeg", width=8, height=7)
  pos <- pos + 1
}

# Creates the plot looking at the uniqueness measure
curr_plot <- make_uniqueness_plot(within_df)
output_name <- paste(working_dir, "within_dataset_unique_stat.png", sep="")
ggsave(output_name, plot=curr_plot, dpi=800, device="jpeg", width=9, height=6)


## SECTION 2 - Looking at data regarding across groups
across_df <- read.csv(across_groups_data_path, header=TRUE)

# Start again, making the stacked bar charts ...
curr_plot <- make_across_groups_plot(across_df)
output_name <- paste(working_dir, "across_dataset_kmer_analysis.png", sep="")
ggsave(output_name, plot=curr_plot, dpi=800, device="jpeg", width=8, height=7)

# Next, make a plot for the uniqueness statistic across groups ...
curr_plot <- make_uniqueness_plot_across_groups(across_df)
output_name <- paste(working_dir, "across_dataset_unique_stat.png", sep="")
ggsave(output_name, plot=curr_plot, dpi=800, device="jpeg", width=7, height=6)
