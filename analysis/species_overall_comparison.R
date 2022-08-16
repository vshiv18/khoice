#####################################################
# Name: species_overall_comparison.R
# Description: Generates line plot for short/long read
#              and in/out genome kmer f1 across k
#              Originally planned for Section 2
#
# Date: August 2, 2022
#
#####################################################

library(ggplot2)
library(ggpubr)
library(data.table)
library(tidyr)
library(dplyr)


########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
#######################################################################

dataset_names <- c("B. Cereus", "B. Anthracis", "B. Thuringiensis", "B. Weihenstephanensis")

working_dir <- "/Users/mwche/Downloads/genome_trials/"

# read working dir should contain csv files with "long"/"short" in title (2 per trial)
# genome working dir should contain csv files with "in"/"out" in title (1 per trial)
# ex. trial_1_in.csv

########################################################################
# Methods for generating accuracy plot
########################################################################

make_f1_ribbon_sd <- function(read_df, genome_df, dataset_name){
  # Makes lineplot with 4 lines (out-short, out-long, out-genome, in-genome)
  
  # Create separate datasets for each
  in_genome = genome_df[genome_df$pivot_type == 'in', ]
  out_genome = genome_df[genome_df$pivot_type == 'out', ]
  short = read_df[read_df$read_type == 'short', ]
  long = read_df[read_df$read_type == 'long', ]
  
  # Find value for vertical line
  max_long = long$k[long$avg == max(long$avg)][1]
  
  plot <- ggplot() +
    # Creates ribbon and line for each
    geom_ribbon(data = short, aes(x=k, ymin=minus_sd, ymax=plus_sd, fill = "Out-Short"),
                alpha=0.2, colour = NA) +
    geom_ribbon(data = long, aes(x=k, ymin=minus_sd, ymax=plus_sd, fill = "Out-Long"),
                alpha=0.2, colour = NA) +
    geom_ribbon(data = out_genome, aes(x=k, ymin=minus_sd, ymax=plus_sd, fill = "Out-Genome"),
                alpha=0.2, colour = NA) +
    geom_ribbon(data = in_genome, aes(x=k, ymin=minus_sd, ymax=plus_sd, fill = "In-Genome"),
                alpha=0.2, colour = NA) +
    geom_line(data = short, aes(x = k, y = avg, color = "Out-Short"), size=1.0) +
    geom_line(data = long, aes(x = k, y = avg, color = "Out-Long"), size=1.0) +
    geom_line(data = out_genome, aes(x = k, y = avg, color = "Out-Genome"), size=1.0) +
    geom_line(data = in_genome, aes(x = k, y = avg, color = "In-Genome"), size=1.0) +
    # Creates long read vertical line
    geom_vline(xintercept=max_long, linetype="dashed", color = "#F8766D", size = 1.0) +
    theme_bw() +
    theme(plot.title=element_text(hjust = 0.5, size=14),
          axis.title.x=element_text(size =14),
          axis.title.y=element_text(size=14),
          legend.position = "none",
          legend.text=element_text(size=14),
          legend.box="horizontal",
          legend.title=element_text(size=12),
          axis.text=element_text(size=12, color="black")) +
    scale_x_continuous(breaks = seq(5, 50, 5)) +
    labs(title = dataset_name,
         x = "Kmer length (k)",
         y = "F1") +
    scale_color_manual(name = "",
                          breaks = c("Out-Short","Out-Long", "Out-Genome","In-Genome"), 
                          values = c("Out-Short"="#00BFC4","Out-Long"= "#F8766D", "Out-Genome"= "#7CAE00","In-Genome"="#C77CFF")) +
    scale_fill_manual(name = "",
                       breaks = c("Out-Short","Out-Long", "Out-Genome"), 
                       values = c("Out-Short"="#00BFC4","Out-Long"= "#F8766D", "Out-Genome"= "#7CAE00", "In-Genome"="#C77CFF")) +
    guides(fill = "none", color = guide_legend(nrow = 1))  
  return(plot)
}



############################################################################
# Start of the main method of code ...
############################################################################

#######################################################################
# WRANGLE READ DATA
#######################################################################

# Get all csv files from working directory
read_working_dir <- "/Users/mwche/Downloads/b_subset_trials/"

values_files <- list.files(path=read_working_dir, pattern = "\\.csv$")
col_names <- c("k","dataset","TP","TN","FP","FN")

value_df_list = list()

for(i in seq_along(values_files)){
  print(values_files[i])
  
  # Read in csv and annotate with name and dataset
  temp_df = read.csv(file = paste(read_working_dir,values_files[i], sep = ""),header = FALSE)
  colnames(temp_df) <- col_names
  temp_df[,"dataset"] <- as.factor(temp_df[,"dataset"])
  temp_df[,'source_csv'] <-values_files[[i]]
  
  # Fill columns with read type
  if(grepl("short", values_files[i],fixed = TRUE)){
    temp_df[,"read_type"] <-"short"
  }else{
    temp_df[,"read_type"] <-"long"
  }
  
  # Calculate recall, precision and f1
  temp_df["recall"] = (temp_df["TP"])/(temp_df["TP"] + temp_df["FN"])
  temp_df["precision"] = (temp_df["TP"])/(temp_df["TP"] + temp_df["FP"])
  temp_df["f1"] = (2 * temp_df["recall"] * temp_df["precision"] / (temp_df["precision"] + temp_df["recall"]))
  value_df_list[[values_files[i]]] <- temp_df
}

merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), value_df_list)

# Aggregate mean, min and max of precision for each k

f1_df <- merge_df %>%
  group_by(k,read_type,dataset) %>%
  summarise_at(vars(f1), list(avg = mean, min = min, max = max, sd = sd))  %>%
  mutate(minus_sd = if_else(avg - sd < 0, 0, avg - sd),
         plus_sd = if_else(avg + sd > 1, 1, avg + sd))

#######################################################################
# WRANGLE GENOME DATA
#######################################################################

genome_working_dir <- "/Users/mwche/Downloads/genome_trials/"
values_files <- list.files(path=genome_working_dir, pattern = "\\.csv$")
col_names <- c("k","dataset","TP","TN","FP","FN")

value_df_list = list()

for(i in seq_along(values_files)){
  print(values_files[i])
  
  # Read in dataframe, and annotate with name
  temp_df = read.csv(file = paste(genome_working_dir,values_files[i], sep = ""),header = FALSE)
  colnames(temp_df) <- col_names
  temp_df[,"dataset"] <- as.factor(temp_df[,"dataset"])
  temp_df[,'source_csv'] <-values_files[[i]]
  
  # Fill columns with read type
  if(grepl("in", values_files[i],fixed = TRUE)){
    temp_df[,"pivot_type"] <-"in"
  }else{
    temp_df[,"pivot_type"] <-"out"
  }
  
  # Calculate recall and precision
  temp_df["recall"] = (temp_df["TP"])/(temp_df["TP"] + temp_df["FN"])
  temp_df["precision"] = (temp_df["TP"])/(temp_df["TP"] + temp_df["FP"])
  temp_df["f1"] = (2 * temp_df["recall"] * temp_df["precision"] / (temp_df["precision"] + temp_df["recall"]))
  value_df_list[[values_files[i]]] <- temp_df
}

merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), value_df_list)

genome_f1_df <- merge_df %>%
  group_by(k,pivot_type,dataset) %>%
  summarise_at(vars(f1), list(avg = mean, min = min, max = max, sd = sd))  %>%
  mutate(minus_sd = if_else(avg - sd < 0, 0, avg - sd),
         plus_sd = if_else(avg + sd > 1, 1, avg + sd))

#######################################################################
# GENERATE PLOTS
#######################################################################

# Create list of plots
plot_list = list()

# Create plot for each dataset
for(dataset in 0:(length(dataset_names)-1)){
  read_dataset_df = read_f1_df[read_f1_df$dataset == as.character(dataset),]
  genome_dataset_df = genome_f1_df[genome_f1_df$dataset == as.character(dataset),]
  temp_plot <- make_f1_ribbon_sd(read_dataset_df, genome_dataset_df,dataset_names[dataset+1])
  mylegend <- get_legend(temp_plot)
  plot_list[[dataset + 1]] <- temp_plot
}


# Combining plots
final_plot <- ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
                        labels = c("a", "b", "c", "d"), ncol = 1, nrow = 4,
                        common.legend = TRUE, legend = "bottom")
print(final_plot)

# Saving plots: a vector and non-vector graphic
output_name <- paste(working_dir, "combined_plot_v2.png", sep="")
ggsave(output_name, plot=final_plot, dpi=800, device="jpeg", width=6, height=12)


output_name <- paste(working_dir, "combined_plot_v2.pdf", sep="")
ggsave(output_name, plot=final_plot, dpi=800, device="pdf", width=11, height=8)


