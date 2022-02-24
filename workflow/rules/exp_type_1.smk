####################################################
# Name: exp_type_1.smk
# Description: Contains functions and rules for
#              the type of experiment type 1: this
#              type is where we generally look at the
#              kmer set across and within groups.
#
# Date: 2/3/22
####################################################

####################################################
# Section 1: Code that is always run, when this 
#            experiment is called. It genereates the
#            needed files.
# 
#            Writes the input files for KMC complex 
#            operations ...
#
# NOTE: the complex operations requires an input file,
#       and this code generates those files, and saves
#       them in the the complex_ops folder.
####################################################

if exp_type == 1:
    # Make tmp directory for kmc
    if not os.path.isdir("tmp"):
        os.mkdir("tmp")

    # Initialize complex_ops directory to start writing files
    if not os.path.isdir("complex_ops"):
        os.mkdir("complex_ops")
        os.mkdir("complex_ops/within_groups")

    # This loop builds the complex_ops files for within groups ...
    for k in k_values:
        if not os.path.isdir(f"complex_ops/within_groups/k_{k}"):
            os.mkdir(f"complex_ops/within_groups/k_{k}")

        for num in range(1, num_datasets+1):
            kmc_input_files = []

            for data_file in os.listdir(f"data/dataset_{num}"):
                if data_file.endswith(".fna.gz"):
                    base_name = data_file.split(".fna.gz")[0]
                    kmc_input_files.append(f"step_2/k_{k}/dataset_{num}/{base_name}.transformed")

            if not os.path.isdir(f"complex_ops/within_groups/k_{k}/dataset_{num}/"):
                os.mkdir(f"complex_ops/within_groups/k_{k}/dataset_{num}")

            with open(f"complex_ops/within_groups/k_{k}/dataset_{num}/within_dataset_{num}.txt", "w") as fd:
                fd.write("INPUT:\n")
                result_str = "("
                for i, path in enumerate(kmc_input_files):
                    fd.write(f"set{i+1} = {path}\n")
                    result_str += "set{} + ".format(i+1)
                result_str = result_str[:-2] + ")"
                fd.write("OUTPUT:\n")
                fd.write(f"step_3/k_{k}/dataset_{num}/dataset_{num}.transformed.combined = {result_str}\n")
                fd.write("OUTPUT_PARAMS:\n-cs5000\n")

    if not os.path.isdir("complex_ops/across_groups"):
        os.mkdir("complex_ops/across_groups")

    # This loop builds the complex_ops across groups ...
    for k in k_values:
        if not os.path.isdir(f"complex_ops/across_groups/k_{k}"):
            os.mkdir(f"complex_ops/across_groups/k_{k}")
        
        kmc_input_files = []
        for i in range(1, num_datasets+1):
            kmc_input_files.append(f"step_6/k_{k}/dataset_{i}/dataset_{i}.transformed.combined.transformed")
        
        with open(f"complex_ops/across_groups/k_{k}/across_all_datasets.txt", "w") as fd:
            fd.write("INPUT:\n")
            result_str = "("
            for i, path in enumerate(kmc_input_files):
                fd.write(f"set{i+1} = {path}\n")
                result_str += "set{} + ".format(i+1)
            result_str = result_str[:-2] + ")"
            fd.write("OUTPUT:\n")
            fd.write(f"step_7/k_{k}/all_datasets.transformed.combined.transformed.combined = {result_str}\n")
            fd.write("OUTPUT_PARAMS:\n-cs5000\n")


####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_group_set(wildcards):
    """ Based on the group, this function will find all the input files """
    input_files = []
    for data_file in os.listdir(f"data/dataset_{wildcards.num}"):
        if data_file.endswith(".fna.gz"):
            file_name = data_file.split(".fna.gz")[0]
            input_files.append(f"step_2/k_{wildcards.k}/dataset_{wildcards.num}/{file_name}.transformed.kmc_pre")
    return input_files

def get_across_group_union(wildcards):
    """ Returns the input files for this union operation - union across all groups """
    k_value = wildcards.k
    input_files = [f"step_6/k_{k_value}/dataset_{num}/dataset_{num}.transformed.combined.transformed.kmc_pre" for num in range(1, num_datasets+1)]
    return input_files

def get_num_of_dataset_members(dataset_num):
    """ Returns the number of genomes in a particular dataset """
    num = 0
    for data_file in os.listdir(f"data/dataset_{dataset_num}"):
        if data_file.endswith(".fna.gz"):
            num += 1
    return num

def summarize_histogram_type1(hist_counts, num_dataset_members, across_group_analysis, k):
    """ 
        Takes in histogram of kmer occurrences, and returns the metrics for the 
        bar charts summarized below:

            %_1_occ - percentage of unique kmers that only occur in one genome
            %_25_or_less - percentage of unique kmers that occur in multiple genomes, in 25% of the genomes or less
            %_25_to_75 - percentage of unique kmers that occur in multiple genomes, in 25% to 75% of the genomes
            %_75_or_more - percentage of unique kmers that occur in multiple genomes, in 75% or more of the genomes
            unique_stat - weighted sum of kmer occurrents = SUM([occ * %_unique_occ for occ in range(255)]) 
    """
    metrics = [0 for i in range(7)]
    total_unique_kmers = sum(hist_counts)
    
    boundaries = [0.25, 0.75]
    boundary_indices = [max(int(percent * num_dataset_members), 1) for percent in boundaries]

    # Special cases where indices are customized ...
    if across_group_analysis:
        boundary_indices = [5, 20]
    
    metrics[0] = round(hist_counts[0]/total_unique_kmers, 3)
    metrics[1] = round(sum([hist_counts[i] for i in range(1, boundary_indices[0])])/total_unique_kmers, 3)
    metrics[2] = round(sum([hist_counts[i] for i in range(boundary_indices[0], boundary_indices[1])])/total_unique_kmers, 3)
    metrics[3] = round(sum([hist_counts[i] for i in range(boundary_indices[1], len(hist_counts))])/total_unique_kmers, 3)

    rounding_error = abs(sum(metrics[0:4])-1) 
    assert rounding_error < 0.05, "Issue occurred with histogram summarization"

    # Both unnormalized, and normalized uniqueness statistics
    metrics[4] = round(sum([((i+1) * (hist_counts[i]/total_unique_kmers)) for i in range(0, len(hist_counts))]), 4)
    metrics[5] = round(sum([(((i+1)/num_dataset_members) * (hist_counts[i]/total_unique_kmers)) for i in range(0, len(hist_counts))]), 4)

    # Calculate the fraction used by the delta measure
    metrics[6] = round(total_unique_kmers/k, 4)
    return metrics

####################################################
# Section 3: Rules needed for this experiment type
####################################################

rule build_kmc_database_on_genome:
    input:
        "data/dataset_{num}/{genome}.fna.gz"
    output:
        "step_1/k_{k}/dataset_{num}/{genome}.kmc_pre", 
        "step_1/k_{k}/dataset_{num}/{genome}.kmc_suf"
    shell:
        "kmc -fm -m64 -k{wildcards.k} -ci1 {input} step_1/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome} tmp/"

rule transform_genome_to_set:
    input: 
        "step_1/k_{k}/dataset_{num}/{genome}.kmc_pre",
        "step_1/k_{k}/dataset_{num}/{genome}.kmc_suf"
    output:
        "step_2/k_{k}/dataset_{num}/{genome}.transformed.kmc_pre",
        "step_2/k_{k}/dataset_{num}/{genome}.transformed.kmc_suf"
    shell:
        "kmc_tools transform step_1/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome} set_counts 1 step_2/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome}.transformed"

rule within_group_union:
    input:
        get_group_set
    output:
        "step_3/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "step_3/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    shell:
        "kmc_tools complex complex_ops/within_groups/k_{wildcards.k}/dataset_{wildcards.num}/within_dataset_{wildcards.num}.txt"

rule within_group_union_histogram:
    input:
        "step_3/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "step_3/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    output:
        "step_4/k_{k}/dataset_{num}/dataset_{num}_k{k}_hist.txt"
    shell:
        "kmc_tools transform step_3/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined histogram step_4/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}_k{wildcards.k}_hist.txt"

rule within_group_union_analysis:
    input:
        expand("step_4/k_{k_len}/dataset_{num}/dataset_{num}_k{k_len}_hist.txt", k_len=k_values, num=list(range(1, num_datasets+1)))
    output:
        "step_5/within_datasets_analysis.csv"
    run:
        with open(output[0], "w") as out_fd:
            out_fd.write(f"group_num,k,percent_1_occ,percent_25_or_less,percent_25_to_75,"
                          "percent_75_or_more,unique_stat,unique_stat_norm,delta_frac,delta_frac_norm\n")
            all_metrics = []
            for input_file in input:
                parts = input_file.split("/")

                k = parts[1][2:]
                dataset_num = parts[2].split("_")[1]
                num_dataset_members = get_num_of_dataset_members(dataset_num)

                with open(input_file, "r") as input_fd:
                    hist_results = input_fd.readlines()
                    hist_counts = [int(record.split()[1]) for record in hist_results]
                
                metrics = [f"group_{dataset_num}", k] + summarize_histogram_type1(hist_counts, num_dataset_members, False, int(k))
                all_metrics.append(metrics)
                            
            # Generate the normalized cardinality/k values for each group
            for i in range(1, num_datasets+1):
                id_str = f"group_{i}"
                values = [metrics[8] for metrics in all_metrics if metrics[0] == id_str]
                max_ratio = max(values)

                # Divide all the values in that group with the max value
                for metrics in all_metrics:
                    if metrics[0] == id_str:
                        metrics.append(round(metrics[8]/max_ratio, 4))
            
            # Print all the values to csv file
            for metrics in all_metrics:
                metrics_str = ",".join([str(x) for x in metrics])
                out_fd.write(f"{metrics_str}\n")

rule build_group_kmer_set:
    input:
        "step_3/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "step_3/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    output:
        "step_6/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.transformed.kmc_pre",
        "step_6/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.transformed.kmc_suf"
    shell:
        "kmc_tools transform step_3/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined set_counts 1 step_6/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined.transformed"
    
rule across_group_union:
    input:
        get_across_group_union
    output:
        "step_7/k_{k}/all_datasets.transformed.combined.transformed.combined.kmc_pre",
        "step_7/k_{k}/all_datasets.transformed.combined.transformed.combined.kmc_suf"
    shell:
        "kmc_tools complex complex_ops/across_groups/k_{wildcards.k}/across_all_datasets.txt"

rule across_group_union_histogram:
    input:
        "step_7/k_{k}/all_datasets.transformed.combined.transformed.combined.kmc_pre",
        "step_7/k_{k}/all_datasets.transformed.combined.transformed.combined.kmc_suf"
    output:
        "step_8/k_{k}/all_datasets_k{k}_hist.txt"
    shell:
        "kmc_tools transform step_7/k_{wildcards.k}/all_datasets.transformed.combined.transformed.combined histogram step_8/k_{wildcards.k}/all_datasets_k{wildcards.k}_hist.txt"


rule across_group_union_analysis:
    input:
        expand("step_8/k_{k_len}/all_datasets_k{k_len}_hist.txt", k_len=k_values)
    output:
        "step_9/across_datasets_analysis.csv"
    run:
        with open(output[0], "w") as out_fd:
            out_fd.write(f"group_num,k,percent_1_occ,percent_2_to_5,percent_5_to_20,percent_20_more,"
                          "unique_stat,unique_stat_norm,delta_frac,delta_frac_norm\n")
            all_metrics = []
            for input_file in input:
                parts = input_file.split("/")

                k = parts[1][2:]
                dataset_num = parts[2].split("_")[1]
                num_dataset_members = num_datasets

                with open(input_file, "r") as input_fd:
                    hist_results = input_fd.readlines()
                    hist_counts = [int(record.split()[1]) for record in hist_results]
                
                metrics = ["full_group", k] + summarize_histogram_type1(hist_counts, num_dataset_members, True, int(k))
                all_metrics.append(metrics)

            # Determine the maximum delta fraction, and use it to normalize
            values = [metrics[8] for metrics in all_metrics]
            max_ratio = max(values)

            # Divide all the values in that group with the max value
            for metrics in all_metrics:
                metrics.append(round(metrics[8]/max_ratio, 4))
            
            # Print all the values to csv file
            for metrics in all_metrics:
                metrics_str = ",".join([str(x) for x in metrics])
                out_fd.write(f"{metrics_str}\n")

rule copy_final_results_type1:
    input:
        "step_5/within_datasets_analysis.csv",
        "step_9/across_datasets_analysis.csv"
    output:
        "final_results_type1/within_datasets_analysis.csv",
        "final_results_type1/across_datasets_analysis.csv"
    run:
        shell("cp {input[0]} {output[0]}")
        shell("cp {input[1]} {output[1]}")

