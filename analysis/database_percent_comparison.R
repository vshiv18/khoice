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
library(ggpubr)

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
  plot <- ggplot(group_df_melt, aes(fill=range, x=database_size, y=percent)) + 
    geom_bar(position="fill", stat="identity") +
    theme_classic() +
    theme(plot.title=element_text(hjust = 0.5, size=14),
          axis.title.x=element_text(size =14),
          axis.title.y=element_text(size=14),
          legend.position = "bottom", 
          legend.text=element_text(size=14),
          legend.box="horizontal",
          legend.title=element_text(size=14),
          axis.text=element_text(size=12, color="black")) +
    scale_x_continuous(breaks = seq(20, 100, 20)) +
    labs(x="Database Size (% RefSeq)",
         y="Ratio of Unique Kmers",
         title=name) +
    scale_fill_discrete(name = "Number of Groups the Kmers Occur In:",
                        labels = c("4 Groups", "2 to 3 Groups", "Only Pivot"))
  return (plot)
}

make_subset_unique_across_plot <- function(group_df, name) {
  # Convert to data.table, melt, and then revert to data.frame
  
  # Build the stacked bar chart ...
  plot <- ggplot(group_df, aes(x=database_size, y=percent_1_occ)) + 
    geom_bar(stat="identity",fill="steelblue") +
    theme_classic() +
    theme(plot.title=element_text(hjust = 0.5, size=14),
          axis.title.x=element_text(size =14),
          axis.title.y=element_text(size=14),
          legend.position = "bottom", 
          legend.text=element_text(size=14),
          legend.box="horizontal",
          legend.title=element_text(size=14),
          axis.text=element_text(size=12, color="black")) +
    scale_x_continuous(breaks = seq(20, 100, 20)) +
    labs(x="Database Size (% RefSeq)",
         y="% Kmers Unique to Pivot",
         title=name) 
  return (plot)
}

make_across_groups_plot <- function(group_df, name) {
  # Convert to data.table, melt, and then revert to data.frame
  setDT(group_df)
  group_df_melt <- melt(group_df, 
                        measure.vars = c("percent_4_to_8","percent_2_to_3", "percent_1_occ"),
                        variable.name = "range", 
                        value.name = "percent")
  print(tail(group_df_melt))
  setDF(group_df_melt)
  
  # Build the stacked bar chart ...
  plot <- ggplot(group_df_melt, aes(fill=range, x=factor(k), y=percent)) + 
    geom_bar(position="fill", stat="identity") +
    theme_classic() +
    theme(plot.title=element_text(hjust = 0.5, size=14),
          axis.title.x=element_text(size =14),
          axis.title.y=element_text(size=14),
          legend.position = "bottom", 
          legend.text=element_text(size=14),
          legend.box="horizontal",
          legend.title=element_text(size=14),
          axis.text=element_text(size=12, color="black")) +
    labs(x="Kmer Length (k)",
         y="Ratio of Unique Kmers",
         title=name) +
    scale_fill_discrete(name = "Number of Groups the Kmers Occur In: ",
                        labels = c("4 Groups", "2 to 3 Groups", "Only Pivot")) +
    scale_x_discrete(breaks = seq(7, 35, 4)) 
  return (plot)
}

############################################################################
# Start of the main method of code ...
############################################################################

#######################################################################
# WRANGLE KMER DATA
#######################################################################

percent_files <- list.files(path=working_dir,pattern = "\\.csv$")

percent_df_list = list()

for(i in seq_along(percent_files)){
  print(percent_files[i])
  temp_df = read.csv(file = paste(working_dir,percent_files[i], sep = ""),header = TRUE)
  temp_df[,"database_size"] <- strtoi(strsplit(strsplit(percent_files[i], "_")[[1]],".csv")[[2]])
  percent_df_list[[percent_files[i]]] <- temp_df
}
merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), percent_df_list)

###################################
# GENERATE PLOTS 
###################################

### Arrange 4 species across k at 20 and 100%
percent_db = 20

plot_list = list()
for(dataset in 0:(length(dataset_names)-1)){
  group = paste("group_",dataset+1, sep = "")
  temp_df = merge_df[merge_df$group_num == group, ]
  temp_df = temp_df[temp_df$database_size == percent_db, ]
  temp_plot <-make_across_groups_plot(temp_df,dataset_names[dataset+1])
  plot_list[[dataset + 1]]<- temp_plot
}

final_plot <- ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
                        labels = c("a", "b", "c", "d"), ncol = 2, nrow = 2,
                        common.legend = TRUE, legend = "bottom")


print(final_plot)
output_name <- paste(working_dir, "bacillus_database_20.jpeg", sep="")
ggsave(output_name, plot=final_plot, dpi=800, device="jpeg", width=10, height=8)

output_name <- paste(working_dir, "bacillus_database_20.pdf", sep="")
ggsave(output_name, plot=final_plot, dpi=800, device="pdf", width=10, height=8)

### Subset for species and k over sizes

subset_list = list()
for(i in seq_along(percent_files)){
  temp_df <- percent_df_list[[percent_files[i]]]
  temp_df <- temp_df[temp_df$group_num == 'group_1', ]
  temp_df <- temp_df[temp_df$k == '31', ]
  subset_list[[percent_files[i]]] <- temp_df
}
subset <- Reduce(function(x, y) merge(x, y, all=TRUE), subset_list)


p_range <- make_subset_across_plot(subset, "B. Cereus")
p_unique <-make_subset_unique_across_plot(subset, "B. Cereus")
print(p_unique)
output_name <- paste(working_dir, "cereus_unique.jpeg", sep="")
ggsave(output_name, plot=p_unique, dpi=800, device="jpeg", width=10, height=8)

output_name <- paste(working_dir, "cereus_unique.pdf", sep="")
ggsave(output_name, plot=p_unique, dpi=800, device="pdf", width=10, height=8)


