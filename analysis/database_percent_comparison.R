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

dataset_names <- c("B. Cereus", "B. Anthracis", "B. Thuringiensis", "B. Weihenstephanensis")
working_dir <- "/Users/mwche/Downloads/b_p_data/"

########################################################################
# Methods for generating the two types of plots: within groups, and 
# across of multiple groups
########################################################################
make_subset_across_plot <- function(group_df, name) {
  # Convert to data.table, melt, and then revert to data.frame
  setDT(group_df)
  group_df_melt <- melt(group_df, 
                        measure.vars = c("percent_4_to_8","percent_2_to_3", "percent_1_occ"),
                        variable.name = "range", 
                        value.name = "percent")
  setDF(group_df_melt)
  
  # Build the stacked bar chart ...
  graph_title <- paste("Uniqueness of Kmers Across All", length(dataset_names), "Groups w.r.t", name, "Pivot")
  plot <- ggplot(group_df_melt, aes(fill=range, x=database_size, y=percent)) + 
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
    labs(x="Database Size",
         y="Ratio of Unique Kmers",
         title=graph_title) +
    scale_fill_discrete(name = "Number of Groups the Kmers Occur In:",
                        labels = c("4 Groups", "2 to 3 Groups", "Only Pivot"))
  return (plot)
}

############################################################################
# Start of the main method of code ...
############################################################################

#######################################################################
# WRANGLE KMER DATA
#######################################################################

percent_files <- list.files(path=working_dir,pattern = "\\.csv$")
col_names <- c("k","dataset","TP","TN","FP","FN")

percent_df_list = list()

for(i in seq_along(percent_files)){
  print(percent_files[i])
  temp_df = read.csv(file = paste(working_dir,percent_files[i], sep = ""),header = TRUE)
  percent_df_list[[percent_files[i]]] <- temp_df
}

percentages = c("20","40","60","80","100")
print(strsplit("trial_100.csv", "_")[1])
strsplit(strsplit("trial_100.csv", "_")[[1]],".csv")[[2]]

anthracis_subset_list = list()
for(i in seq_along(percent_files)){
  temp_df <- percent_df_list[[percent_files[i]]]
  temp_df <- temp_df[temp_df$group_num == 'group_2', ]
  temp_df <- temp_df[temp_df$k == '31', ]
  temp_df[,"database_size"] <- as.factor(strsplit(strsplit(percent_files[i], "_")[[1]],".csv")[[2]])
  anthracis_subset_list[[percent_files[i]]] <- temp_df
}
anthracis_subset <- Reduce(function(x, y) merge(x, y, all=TRUE), anthracis_subset_list)

make_subset_across_plot(anthracis_subset, "B. Anthracis")


