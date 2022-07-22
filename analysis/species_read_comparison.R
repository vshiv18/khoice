#####################################################
# Name: kmer_plots_type_6.R
# Description: Plots accuracy across values of k for
# short and long simulated reads
#
# Date: July 5, 2022
#
#####################################################

library(ggplot2)
library(ggpubr)
library(data.table)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)
install.packages('gridExtra')

########################################################################
# IMPORTANT: Experiment-dependent variables below, need to be set ...
#######################################################################

dataset_names <- c("B. Cereus", "B. Anthracis", "B. Thuringiensis", "B. Weihenstephanensis")
#dataset_names <- c("E. Coli", "Salmonella")

#working_dir <- "/Users/mwche/Downloads/exp6_bacillus/"
working_dir <- "/Users/omarahmed/downloads/raw_data/b_subset_trials/"

########################################################################
# Methods for generating accuracy plot
########################################################################

make_recall_ribbon_sd <- function(kmer_df, dataset_name, mem_df){
  hm = mem_df[mem_df$mem_type == 'half_mem', ]
  m = mem_df[mem_df$mem_type == 'mem', ]
  
  #plot_title = paste(dataset_name,"Pivot")
  plot_title <- dataset_name
  plot <- ggplot(kmer_df, aes(color = read_type)) +
          geom_ribbon(aes(x=k, ymin=minus_sd, ymax=plus_sd, fill = read_type),
                alpha=0.2, colour = NA) +
          geom_line(aes(x = k, y = avg, color = read_type, linetype = "Kmer"), size=1.0) +
          geom_hline(aes(yintercept=avg, linetype="Half-MEM", color=read_type), data = hm, size=1.0) +
          geom_hline(aes(yintercept=avg, linetype="MEM", color=read_type), data = m, size=1.0) +
          theme_bw() +
          theme(plot.title=element_text(hjust = 0.5, size=14),
          axis.title.x=element_text(size =14),
          axis.title.y=element_text(size=14),
          legend.position = "bottom", 
          legend.text=element_text(size=14),
          legend.box="horizontal",
          legend.title=element_text(size=12),
          axis.text=element_text(size=12, color="black")) +
          scale_y_continuous(breaks=seq(0, 1.0, 0.1)) +
          scale_x_continuous(breaks = seq(5, 50, 5)) +
          labs(title = plot_title,
               x = "Kmer length (k)",
               y = "Recall") +
          scale_linetype_manual(name = "",
                                breaks = c("Kmer","Half-MEM", "MEM"), 
                                values = c("Kmer"="solid","Half-MEM"= "dashed", "MEM"= "dotted")) +
          scale_color_discrete(name = "", labels = c("Long Read", "Short Read")) +
          guides(fill = "none", color = guide_legend(nrow = 1))  
          #theme(plot.margin = margin(t = 20,  # Top margin
          #                          r = 50,  # Right margin
          #                          b = 40,  # Bottom margin
          #                          l = 10)) # Left margin)

  return(plot)
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

############################################################################
# Start of the main method of code ...
############################################################################

#######################################################################
# WRANGLE KMER DATA
#######################################################################

values_files <- list.files(path=working_dir, pattern = "\\.csv$")
col_names <- c("k","dataset","TP","TN","FP","FN")

value_df_list = list()

for(i in seq_along(values_files)){
  print(values_files[i])
  
  # Read in dataframe, and annotate with name
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


# Aggregate mean, min and max of precision for each k

recall_df <- merge_df %>%
  group_by(k,read_type,dataset) %>%
  summarise_at(vars(recall), list(avg = mean, min = min, max = max, sd = sd)) 
recall_df["minus_sd"] = recall_df["avg"] - recall_df["sd"]
recall_df["plus_sd"] = recall_df["avg"] + recall_df["sd"]

#######################################################################
# WRANGLE MEM DATA
#######################################################################

mem_working_dir <- "/Users/omarahmed/Downloads/raw_data/b_mem_trials/"
#mem_working_dir <- "/Users/mwche/Downloads/b_mem_trials/"
mem_values_files <- list.files(path=mem_working_dir,pattern = "\\.csv$")
mem_col_names <- c("dataset","TP","TN","FP","FN")

mem_value_df_list = list()

for(i in seq_along(mem_values_files)){
  print(mem_values_files[i])
  
  # Read in dataframe, and annotate with name
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

# Aggregate mean, min and max of recall for mem and read type

mem_recall_df <- mem_merge_df %>%
  group_by(mem_type,read_type,dataset) %>%
  summarise_at(vars(recall), list(avg = mean, min = min, max = max, sd = sd)) 
mem_recall_df["minus_sd"] = mem_recall_df["avg"] - mem_recall_df["sd"]
mem_recall_df["plus_sd"] = mem_recall_df["avg"] + mem_recall_df["sd"]

#######################################################################
# GENERATE PLOTS
#######################################################################

# Create list of plots
plot_list = list()

for(dataset in 0:(length(dataset_names)-1)){
  print(dataset)
  kmer_dataset_df = recall_df[recall_df$dataset == as.character(dataset),]
  mem_dataset_df = mem_recall_df[mem_recall_df$dataset == as.character(dataset),]
  temp_plot <- make_recall_ribbon_sd(kmer_dataset_df, dataset_names[dataset+1], mem_dataset_df)
  plot_list[[dataset + 1]]<- temp_plot
}

# Combining plots
final_plot <- ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
                        labels = c("a", "b", "c", "d"), ncol = 2, nrow = 2,
                        common.legend = TRUE, legend = "bottom")
print(final_plot)

# Saving plots: a vector and non-vector graphic
output_name <- paste(mem_working_dir, "combined_plot_v2.jpeg", sep="")
ggsave(output_name, plot=final_plot, dpi=800, device="jpeg", width=10, height=8)

output_name <- paste(mem_working_dir, "combined_plot_v2.pdf", sep="")
ggsave(output_name, plot=final_plot, dpi=800, device="pdf", width=10, height=8)


# Manually creating plots
#dataset = 0
#kmer_dataset_df = recall_df[recall_df$dataset == as.character(dataset),]
#mem_dataset_df = mem_recall_df[mem_recall_df$dataset == as.character(dataset),]
#p1 <- make_recall_ribbon_sd(kmer_dataset_df, dataset_names[dataset+1], mem_dataset_df)

#print(p4)

# Arranging plots
#mylegend<-g_legend(p1)
#p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
#                         p2 + theme(legend.position="none"),
#                         p3 + theme(legend.position="none"),
#                         p4 + theme(legend.position="none"),
#                         nrow=2),
#                        mylegend, nrow=2,heights=c(10, 1),
#                  top = textGrob("Classification Recall w.r.t Out-Pivot",gp=gpar(fontsize=20)))

#output_name <- paste(mem_working_dir, "combined_plot_v1.png", sep="")
#ggsave(output_name, plot=p, dpi=800, device="jpeg", width=12, height=10)


