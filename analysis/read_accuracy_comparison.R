#####################################################
# Name: kmer_plots_type_6.R
# Description: Plots accuracy across values of k for
# short and long simulated reads
#
# Date: July 5, 2022
#
#####################################################
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
#######################################################################

dataset_names <- c("B. Cereus", "B. Anthracis", "B. Thuringiensis", "B. Weihenstephanensis ")
#dataset_names <- c("E. Coli", "Salmonella")

#working_dir <- "/Users/mwche/Downloads/exp6_bacillus/"
working_dir <- "/Users/mwche/Downloads/b_subset_trials/"

########################################################################
# Methods for generating accuracy plot
########################################################################

make_recall_ribbon_sd <- function(recall_df, dataset_names){
  plot <- ggplot(recall_df, color = dataset) +
    geom_ribbon(aes(x=k, ymin=minus_sd, ymax=plus_sd, fill = dataset),
                alpha=0.2, colour = NA) +
    geom_line(aes(x = k, y = avg, color = dataset)) +
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
    scale_x_continuous(breaks = seq(5, 50, 5)) +
    labs(x = "Kmer length (k)",
         y = "Recall") +
    scale_fill_discrete(name = "+/- 1 Standard Deviation", labels = dataset_names) +
    guides(colour = 'none') 
  return(plot)
}


############################################################################
# Start of the main method of code ...
############################################################################

#######################################################################
# WRANGLE KMER DATA
#######################################################################

values_files <- list.files(path=working_dir,pattern = "\\.csv$")
col_names <- c("k","dataset","TP","TN","FP","FN")

value_df_list = list()

for(i in seq_along(values_files)){
  print(values_files[i])
  temp_df = read.csv(file = paste(working_dir,values_files[i], sep = ""),header = FALSE)
  colnames(temp_df) <- col_names
  temp_df[,"dataset"] <- as.factor(temp_df[,"dataset"])
  temp_df[,'source_csv'] <-values_files[[i]]
  # Fill columns with read type
  if(grepl("short", values_files[i],fixed = TRUE)){
    temp_df[,"read_type"] <-"short"
  }else{
    temp_df[,"read_type"] <-"long"
  }
  # Calculate recall and precision
  temp_df["recall"] = (temp_df["TP"])/(temp_df["TP"] + temp_df["FN"])
  temp_df["precision"] = (temp_df["TP"])/(temp_df["TP"] + temp_df["FP"])
  value_df_list[[values_files[i]]] <- temp_df
}

merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), value_df_list)

write.csv(merge_df,paste(working_dir,"/merge_df.csv", sep = ""), row.names = FALSE)

# Aggregate mean, min and max of precision for each k

recall_df <- merge_df %>%
  group_by(k,read_type,dataset) %>%
  summarise_at(vars(recall), list(avg = mean, min = min, max = max, sd = sd)) 
recall_df["minus_sd"] = recall_df["avg"] - recall_df["sd"]
recall_df["plus_sd"] = recall_df["avg"] + recall_df["sd"]

#######################################################################
# WRANGLE MEM DATA
#######################################################################

mem_working_dir <- "/Users/mwche/Downloads/b_mem_trials/"
mem_values_files <- list.files(path=mem_working_dir,pattern = "\\.csv$")
mem_col_names <- c("dataset","TP","TN","FP","FN")

mem_value_df_list = list()

for(i in seq_along(mem_values_files)){
  print(mem_values_files[i])
  temp_df = read.csv(file = paste(mem_working_dir,mem_values_files[i], sep = ""),header = FALSE)
  colnames(temp_df) <- mem_col_names
  temp_df[,"dataset"] <- as.factor(temp_df[,"dataset"])
  temp_df[,'source_csv'] <-mem_values_files[[i]]
  # Fill columns with read type
  if(grepl("short", mem_values_files[i],fixed = TRUE)){
    temp_df[,"read_type"] <-"short"
  }else{
    temp_df[,"read_type"] <-"long"
  }
  if(grepl("hm", mem_values_files[i],fixed = TRUE)){
    temp_df[,"mem_type"] <-"half_mem"
  }else{
    temp_df[,"mem_type"] <-"mem"
  }
  # Calculate recall and precision
  temp_df["recall"] = (temp_df["TP"])/(temp_df["TP"] + temp_df["FN"])
  temp_df["precision"] = (temp_df["TP"])/(temp_df["TP"] + temp_df["FP"])
  mem_value_df_list[[mem_values_files[i]]] <- temp_df
}

mem_merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), mem_value_df_list)

# Aggregate mean, min and max of precision for mem and read type
mem_precision_df <- mem_merge_df %>%
  group_by(mem_type,read_type,dataset) %>%
  summarise_at(vars(precision), list(avg = mean, min = min, max = max, sd = sd)) 
mem_precision_df["minus_sd"] = mem_precision_df["avg"] - mem_precision_df["sd"]
mem_precision_df["plus_sd"] = mem_precision_df["avg"] + mem_precision_df["sd"]

mem_recall_df <- mem_merge_df %>%
  group_by(mem_type,read_type,dataset) %>%
  summarise_at(vars(recall), list(avg = mean, min = min, max = max, sd = sd)) 
mem_recall_df["minus_sd"] = mem_recall_df["avg"] - mem_recall_df["sd"]
mem_recall_df["plus_sd"] = mem_recall_df["avg"] + mem_recall_df["sd"]

#######################################################################
# GENERATE PLOTS
#######################################################################

short_k_df = recall_df[recall_df$read_type == 'short',]
long_k_df = recall_df[recall_df$read_type == 'long',]

short_plot = make_recall_ribbon_sd(short_k_df, dataset_names)
long_plot = make_recall_ribbon_sd(long_k_df, dataset_names)

half_mem_recall = mem_recall_df[mem_recall_df$mem_type == 'half_mem', ]

long_0 = half_mem_recall[half_mem_recall$read_type == 'long'& half_mem_recall$dataset=='0', ]
long_1 = unlist(half_mem_recall[half_mem_recall$read_type == 'long'& half_mem_recall$dataset=='1', ]['avg'])
long_2 = unlist(half_mem_recall[half_mem_recall$read_type == 'long'& half_mem_recall$dataset=='2', ]['avg'])
long_3 = unlist(half_mem_recall[half_mem_recall$read_type == 'long'& half_mem_recall$dataset=='3', ]['avg'])

short_0 = unlist(half_mem_recall[half_mem_recall$read_type == 'short'& half_mem_recall$dataset=='0', ]['avg'])
short_1 = unlist(half_mem_recall[half_mem_recall$read_type == 'short'& half_mem_recall$dataset=='1', ]['avg'])
short_2 = unlist(half_mem_recall[half_mem_recall$read_type == 'short'& half_mem_recall$dataset=='2', ]['avg'])
short_3 = unlist(half_mem_recall[half_mem_recall$read_type == 'short'& half_mem_recall$dataset=='3', ]['avg'])

comb_long_plot = long_plot + 
  geom_hline(yintercept=unlist(long_0["avg"]), linetype="dashed", color = "#F8766D") +
  #annotate(geom = "rect", xmin = 5, xmax = 50, ymin = unlist(long_0["minus_sd"]), ymax = unlist(long_0["plus_sd"]),
           #fill = "#F8766D", alpha = 0.5) +
  geom_hline(yintercept=long_1, linetype="dashed", color = "#7CAE00") +
  geom_hline(yintercept=long_2, linetype="dashed", color = "#00BFC4") +
  geom_hline(yintercept=long_3, linetype="dashed", color = "#C77CFF") +
  labs(title = "Kmer vs. Mem Recall w.r.t Long-Read Out-Pivot")

comb_short_plot = short_plot + 
  geom_hline(yintercept=short_0, linetype="dashed", color = "#F8766D") +
  geom_hline(yintercept=short_1, linetype="dashed", color = "#7CAE00") +
  geom_hline(yintercept=short_2, linetype="dashed", color = "#00BFC4") +
  geom_hline(yintercept=short_3, linetype="dashed", color = "#C77CFF") +
  labs(title = "Kmer vs. Mem Recall w.r.t Short-Read Out-Pivot")

mem_recall = mem_recall_df[mem_recall_df$mem_type == 'mem', ]

long_0 = unlist(mem_recall[mem_recall$read_type == 'long'& mem_recall$dataset=='0', ]['avg'])
long_1 = unlist(mem_recall[mem_recall$read_type == 'long'& mem_recall$dataset=='1', ]['avg'])
long_2 = unlist(mem_recall[mem_recall$read_type == 'long'& mem_recall$dataset=='2', ]['avg'])
long_3 = unlist(mem_recall[mem_recall$read_type == 'long'& mem_recall$dataset=='3', ]['avg'])

short_0 = unlist(mem_recall[mem_recall$read_type == 'short'& mem_recall$dataset=='0', ]['avg'])
short_1 = unlist(mem_recall[mem_recall$read_type == 'short'& mem_recall$dataset=='1', ]['avg'])
short_2 = unlist(mem_recall[mem_recall$read_type == 'short'& mem_recall$dataset=='2', ]['avg'])
short_3 = unlist(mem_recall[mem_recall$read_type == 'short'& mem_recall$dataset=='3', ]['avg'])

comb_long_plot = comb_long_plot + 
  geom_hline(yintercept=unlist(long_0["avg"]), linetype="dotted", color = "#F8766D") +
  #annotate(geom = "rect", xmin = 5, xmax = 50, ymin = unlist(long_0["minus_sd"]), ymax = unlist(long_0["plus_sd"]),
  #fill = "#F8766D", alpha = 0.5) +
  geom_hline(yintercept=long_1, linetype="dotted", color = "#7CAE00") +
  geom_hline(yintercept=long_2, linetype="dotted", color = "#00BFC4") +
  geom_hline(yintercept=long_3, linetype="dotted", color = "#C77CFF") +
  labs(title = "Kmer vs. Half Mem Recall w.r.t Long-Read Out-Pivot")

print(long_2)

comb_short_plot = comb_short_plot + 
  geom_hline(yintercept=short_0, linetype="dotted", color = "#F8766D") +
  geom_hline(yintercept=short_1, linetype="dotted", color = "#7CAE00") +
  geom_hline(yintercept=short_2, linetype="dotted", color = "#00BFC4") +
  geom_hline(yintercept=short_3, linetype="dotted", color = "#C77CFF") +
  labs(title = "Kmer vs. Half Mem Recall w.r.t Short-Read Out-Pivot")


print(comb_long_plot)
print(comb_short_plot)



########################################################################
# SAVE PLOTS
########################################################################

# Save plots
output_name <- paste(mem_working_dir, "long_recall.png", sep="")
ggsave(output_name, plot=comb_long_plot, dpi=800, device="jpeg", width=9, height=6)

output_name <- paste(mem_working_dir, "short_recall.png", sep="")
ggsave(output_name, plot=comb_short_plot, dpi=800, device="jpeg", width=9, height=6)


