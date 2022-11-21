#####################################################
# Name: species_read_comparison.R
# Description: Generates line/bar plot pairing
#              to compare f1 of kmers/half-MEM/MEMs
#              Originally planned for Section 3
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
kmer_working_dir <- "/Users/omarahmed/downloads/current_research/khoice_exps/results/section_3_plots/trial_6/kmer_data/"

# Kmer working dir should contain csv files with "long"/"short" in title (2 per trial)
# MEM working dir should contain csv files with "long"/"short" AND "hm"/"m" in title (4 per trial)

########################################################################
# Methods for generating accuracy plot
########################################################################

make_f1_ribbon_sd <- function(kmer_df, dataset_name, mem_df){
  # Makes line plot of kmer read f1 across k, FOR A SINGULAR SPECIES DATASET
  
  # Commented out plots MEM/half-MEMS as horizontal lines
  #hm = mem_df[mem_df$mem_type == 'half_mem', ]
  #m = mem_df[mem_df$mem_type == 'mem', ]
  
  plot_title <- dataset_name
  plot <- ggplot(kmer_df, aes(color = read_type)) +
          geom_ribbon(aes(x=k, ymin=f1_avg+f1_sd, ymax=f1_avg-f1_sd, fill = read_type),
                alpha=0.2, colour = NA) +
          geom_line(aes(x = k, y = f1_avg, color = read_type, linetype="solid"), size=1.0) +
          geom_ribbon(aes(x=k, ymin=f12_avg+f12_sd, ymax=f12_avg-f12_sd, fill = read_type),
                alpha=0.2, colour = NA) +
          geom_line(aes(x = k, y = f12_avg, color = read_type, linetype="dashed"), size=1.0) +
          #geom_hline(aes(yintercept=avg, linetype="Half-MEM", color=read_type), data = hm, size=1.0) +
          #geom_hline(aes(yintercept=avg, linetype="MEM", color=read_type), data = m, size=1.0) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=20),
          axis.title.x=element_text(size =20),
          axis.title.y=element_text(size=20),
          legend.position = "none",
          legend.text=element_text(size=20),
          legend.box="horizontal",
          legend.title=element_text(size=20),
          axis.text=element_text(size=14, color="black")) +
          scale_x_continuous(breaks = seq(5, 50, 5)) +
          labs(title = plot_title,
               x = "Kmer length (k)",
               y = "F1") +
          #scale_linetype_manual(name = "",
          #                      breaks = c("Kmer","Half-MEM", "MEM"), 
          #                      values = c("Kmer"="solid","Half-MEM"= "dashed", "MEM"= "dotted")) +
          scale_color_discrete(name = "", labels = c("Long Read", "Short Read")) +
          guides(fill = "none", color = guide_legend(nrow = 1))  

  return(plot)
}

make_f1 <- function(kmer_df, dataset_names){
  # Makes line plot of kmer read f1 across k, FOR ALL DATASETS ON ONE 
  
  plot <- ggplot(kmer_df, aes(group = interaction(read_type, dataset), color = dataset)) +
          geom_line(aes(x =k, y = f1, linetype = read_type), size=1.0) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=14, face="bold"),
                axis.title.x=element_text(size =14),
                axis.title.y=element_text(size=14),
                legend.position = "bottom", 
                legend.text=element_text(size=12),
                legend.box="vertical",
                legend.title=element_text(size=12),
                axis.text=element_text(size=12, color="black")) +
          scale_x_continuous(breaks = seq(5, 50, 5)) +
          labs(x = "Kmer length (k)",
               y = "f1") +
          scale_color_discrete(name = "", labels = dataset_names)
  return(plot)
}


make_f1_bar_chart_sd <- function(dataset_name, mem_df) {
  # Makes bar chart of half-MEM/MEM f1 for read types
  
  plot <- ggplot(mem_df, aes(fill=read_type, y=avg, x=mem_type)) + 
          geom_bar(position="dodge", stat="identity", color="black", width=0.75, size=0.75) +
          geom_errorbar(aes(ymin=minus_sd, ymax=plus_sd, group=interaction(mem_type, read_type)), width=.2, position=position_dodge(0.75), size=0.5) +
          theme_classic() +
          theme(plot.title=element_text(hjust = 0.5, size=14),
                axis.title.x=element_text(size =14),
                axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.line.y=element_blank(),
                legend.position = "bottom", 
                legend.text=element_text(size=14),
                legend.box="horizontal",
                legend.title=element_text(size=12),
                axis.text=element_text(size=12, color="black")) +           
          labs(title = "", x = "", y = "Recall") + 
          scale_fill_discrete(name = "", labels = c("Long Read", "Short Read")) +
          scale_x_discrete(name="", labels=c("Half-MEM", "MEM"))
    
  return(plot)
}


#######################################################################
# Step 1: Wrangle k-mer data into one data-frame ...
#######################################################################

# Get all csv files from working directory
kmer_values_files <- list.files(path=kmer_working_dir, pattern = "\\.csv$")
#col_names <- c("k","dataset","TP","TN","FP","FN")
col_names <- c("k","dataset","TP","TN","FP","FN","TP-U","TN-U","FP-U","FN-U")

value_df_list = list()

# Wrangle each csv and add to list
for(i in seq_along(kmer_values_files)){
  print(kmer_values_files[i])
  
  # Read in csv and annotate with name and dataset
  temp_df = read.csv(file = paste(kmer_working_dir,kmer_values_files[i], sep = ""),header=TRUE)
  colnames(temp_df) <- col_names
  temp_df[,"dataset"] <- as.factor(temp_df[,"dataset"])
  temp_df[,'source_csv'] <- kmer_values_files[[i]]
  
  # Fill columns with read type
  if(grepl("short", kmer_values_files[i],fixed = TRUE)){
    temp_df[,"read_type"] <-"short"
  }else{
    temp_df[,"read_type"] <-"long"
  }
  
  print(temp_df)
  
  # Calculate recall, precision and f1 for both confusion matrices...
  temp_df["recall"] = (temp_df["TP"])/(temp_df["TP"] + temp_df["FN"])
  temp_df["precision"] = (temp_df["TP"])/(temp_df["TP"] + temp_df["FP"])
  temp_df["f1"] = (2 * temp_df["recall"] * temp_df["precision"] / (temp_df["precision"] + temp_df["recall"]))

  temp_df["recall2"] = (temp_df["TP-U"])/(temp_df["TP-U"] + temp_df["FN-U"])
  temp_df["precision2"] = (temp_df["TP-U"])/(temp_df["TP-U"] + temp_df["FP-U"])
  temp_df["f12"] = (2 * temp_df["recall2"] * temp_df["precision2"] / (temp_df["precision2"] + temp_df["recall2"]))

  value_df_list[[kmer_values_files[i]]] <- temp_df
}

# Merge list to mega-df
merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), value_df_list)

# Summarize merged df with avg, sd, +/- sd
f1_df <- merge_df %>%
         group_by(k,read_type,dataset) %>%
        summarise_at(c("f1", "f12"), list(avg=mean, min=min, max=max, sd=sd)) 
   #     mutate(minus_sd = if_else(mean - sd < 0, 0, mean - sd),
   #             plus_sd = if_else(mean + sd > 1, 1, mean + sd))

#######################################################################
# Step 2: Wrangle half-MEMs/MEMs data into one data-frame ...
#######################################################################

# Get all csv files from working directory
mem_working_dir <- "/Users/omarahmed/downloads/current_research/khoice_exps/results/section_3_plots/trial_6/mem_data/"
mem_values_files <- list.files(path=mem_working_dir,pattern = "\\.csv$")
mem_col_names <- c("dataset","TP","TN","FP","FN")

mem_value_df_list = list()

# Wrangle each csv and add to list
for(i in seq_along(mem_values_files)){
  print(mem_values_files[i])
  
  # Read in csv and annotate with name and dataset
  temp_df = read.csv(file = paste(mem_working_dir,mem_values_files[i], sep = ""),header = FALSE)
  colnames(temp_df) <- mem_col_names
  temp_df[,"dataset"] <- as.factor(temp_df[,"dataset"])
  temp_df[,'source_csv'] <-mem_values_files[[i]]
  
  # Fill columns with read type and mem type
  if(grepl("illumina", mem_values_files[i],fixed = TRUE)){
    temp_df[,"read_type"] <-"short"
  }else{
    temp_df[,"read_type"] <-"long"
  }
  if(grepl("half_mems", mem_values_files[i],fixed = TRUE)){
    temp_df[,"mem_type"] <-"half_mem"
  }else{
    temp_df[,"mem_type"] <-"mem"
  }
  
  # Calculate recall, precision and f1
  temp_df["recall"] = (temp_df["TP"])/(temp_df["TP"] + temp_df["FN"])
  temp_df["precision"] = (temp_df["TP"])/(temp_df["TP"] + temp_df["FP"])
  temp_df["f1"] = (2 * temp_df["recall"] * temp_df["precision"] / (temp_df["precision"] + temp_df["recall"]))
  mem_value_df_list[[mem_values_files[i]]] <- temp_df
}

# Merge list to mega-df
mem_merge_df <- Reduce(function(x, y) merge(x, y, all=TRUE), mem_value_df_list)

# Summarize merged df with avg, sd, +/- sd
mem_f1_df <- mem_merge_df %>%
  group_by(mem_type,read_type,dataset) %>%
  summarise_at(vars(f1), list(avg = mean, min = min, max = max, sd = sd))  %>%
  mutate(minus_sd = if_else(avg - sd < 0, 0, avg - sd),
         plus_sd = if_else(avg + sd > 1, 1, avg + sd))

#######################################################################
# Step 3: Generate the plots ...
#######################################################################

# Create list of plots
plot_list = list()

for(dataset in 0:(length(dataset_names)-1)){
  print(dataset)
  
  # Build the dataframes for kmers/half-mems and mems
  kmer_dataset_df = f1_df[f1_df$dataset == as.character(dataset),]
  mem_dataset_df = mem_f1_df[mem_f1_df$dataset == as.character(dataset),]

  # Build the ribbon and bar chart
  temp_plot1 <- make_f1_ribbon_sd(kmer_dataset_df, dataset_names[dataset+1], mem_dataset_df)
  temp_plot2 <- make_f1_bar_chart_sd(dataset_names[dataset+1], mem_dataset_df)
  
  mylegend <- get_legend(temp_plot1)

  # Determine the y-axis range we want to keep the plots on same scale
  max_value <- max(max(kmer_dataset_df["f1_avg"]+kmer_dataset_df["f1_sd"]), max(mem_dataset_df["plus_sd"]))
  max_value <- ceiling(max_value/0.10)*0.10
  min_value <- 0.0

  limits <- c(min_value, max_value)
  breaks <- seq(limits[1], limits[2], by=.1)
  
  # Place the same y-axis in both plots
  temp_plot1 <- temp_plot1 + scale_y_continuous(limits=limits, breaks=breaks)
  temp_plot2 <- temp_plot2 + scale_y_continuous(limits=limits, breaks=breaks)
  plot_list[[dataset + 1]] <- ggarrange(temp_plot1, temp_plot2 + theme(legend.position="none"), ncol=2, widths=c(2,0.85))

}


# Combining plots
final_plot <- ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
                        labels = c("a", "b", "c", "d"), ncol = 2, nrow = 2, 
                        legend.grob=mylegend, 
                        legend="bottom",
                        common.legend=TRUE,
                        font.label=list(size = 20))
print(final_plot)

# Saving plots: a vector and non-vector graphic

plot_output_dir <- "/Users/omarahmed/downloads/current_research/khoice_exps/results/section_3_plots/trial_6/plots/"
output_name <- paste(plot_output_dir, "combined_plot.jpeg", sep="")
ggsave(output_name, plot=final_plot, dpi=800, device="jpeg", width=11, height=8)

output_name <- paste(plot_output_dir, "combined_plot.pdf", sep="")
ggsave(output_name, plot=final_plot, dpi=800, device="pdf", width=11, height=8)


