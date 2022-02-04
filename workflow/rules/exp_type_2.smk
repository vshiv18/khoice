####################################################
# Name: exp_type_2.smk
# Description: Contains functions and rules for
#              the type of experiment type 2: this
#              type is where we generally look at the
#              kmer set in genomes w.r.t a pivot genome.
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

if exp_type == 2:
    # Make tmp directory for kmc
    if not os.path.isdir("tmp"):
        os.mkdir("tmp")
    
    # Parse out a pivot from the downloaded dataset
    if not os.path.isdir("input_data/"):
        for i in range(1, num_datasets+1):
            dir_prefix = f"data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"input_data/rest_of_set/dataset_{i}"):
                os.makedirs(f"input_data/rest_of_set/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"input_data/rest_of_set/dataset_{i}")
            
            # Extract the pivot and move to new directory
            dir_prefix = f"input_data/rest_of_set/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            pivot = random.choice(list_of_files)
            if not os.path.isdir(f"input_data/pivot/dataset_{i}"):
                os.makedirs(f"input_data/pivot/dataset_{i}")
            
            # Remove pivot from other set of genomes
            shutil.copy(pivot, f"input_data/pivot/dataset_{i}/pivot_{i}.fna.gz")
            os.remove(pivot)
    
    # Initialize complex ops directory to start writing files
    if not os.path.isdir("complex_ops"):
        os.makedirs("complex_ops/within_groups")
    
    for k in k_values:
        if not os.path.isdir(f"complex_ops/within_groups/k_{k}"):
            os.mkdir(f"complex_ops/within_groups/k_{k}")
        
        for num in range(1, num_datasets+1):
            kmc_input_files = []

            for data_file in os.listdir(f"input_data/rest_of_set/dataset_{num}"):
                if data_file.endswith(".fna.gz"):
                    base_name = data_file.split(".fna.gz")[0]
                    kmc_input_files.append(f"genome_sets/rest_of_set/k_{k}/dataset_{num}/{base_name}.transformed")
            
            if not os.path.isdir(f"complex_ops/within_groups/k_{k}/dataset_{num}/"):
                os.mkdir(f"complex_ops/within_groups/k_{k}/dataset_{num}/")
        
            with open(f"complex_ops/within_groups/k_{k}/dataset_{num}/within_dataset_{num}.txt", "w") as fd:
                fd.write("INPUT:\n")
                result_str = "("
                for i, path in enumerate(kmc_input_files):
                    fd.write(f"set{i+1} = {path}\n")
                    result_str += "set{} + ".format(i+1)
                result_str = result_str[:-2] + ")"
                fd.write("OUTPUT:\n")
                fd.write(f"within_databases/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined = {result_str}\n")
                fd.write("OUTPUT_PARAMS:\n-cs5000\n")

####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_all_genomes_in_dataset(wildcards):
    """ Returns a list of database in a certain dataset """
    input_files = []
    for data_file in os.listdir(f"input_data/rest_of_set/dataset_{wildcards.num}/"):
        if data_file.endswith(".fna.gz"):
            file_name = data_file.split(".fna.gz")[0]
            input_files.append(f"genome_sets/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{file_name}.transformed.kmc_pre")
    return input_files

def get_dataset_pivot(wildcards):
    """ Return a path to the pivot for certain dataset """
    pivot = os.listdir(f"input_data/pivot/dataset_{wildcards.num}/")
    assert len(pivot) == 1
    return [f"genome_sets/pivot/k_{wildcards.k}/dataset_{wildcards.num}/" + pivot[0].split(".fna.gz")[0] + ".transformed.kmc_pre"]

def get_within_group_histogram_files(wildcards):
    """ Returns all the histogram files for within group experiment in specific order """
    input_files = []
    for num in range(1, num_datasets+1):
        for k in k_values:
            for op in ['subtract', 'intersect']:
                input_files.append(f"within_dataset_results/k_{k}/dataset_{num}/{op}/dataset_{num}_pivot_{op}_group.hist.txt")
    return input_files

def summarize_histogram(sub_counts, inter_counts, num_genomes_in_dataset, across_group_analysis):
    """ 
        Takes in histogram of kmer occurrences, and returns the metrics for the 
        bar charts summarized below:

            %_1_occ - percentage of unique kmers that only occur in one genome
            %_25_or_less - percentage of unique kmers that occur in multiple genomes, in 25% of the genomes or less
            %_25_to_75 - percentage of unique kmers that occur in multiple genomes, in 25% to 75% of the genomes
            %_75_or_more - percentage of unique kmers that occur in multiple genomes, in 75% or more of the genomes
            unique_stat - weighted sum of kmer occurrents = SUM([occ * %_unique_occ for occ in range(255)]) 
    """
    # Perform a couple of assertions to start ...
    assert inter_counts[0] == 0, "intersection counts should have 0 unique kmers"
    assert sum(sub_counts[1:]) == 0, "all of kmers in sub_counts should be unique"
    print(inter_counts[:10])
    print(sub_counts[:10])

    # Start to calculate the metrics ...
    metrics = [0, 0, 0, 0, 0]
    total_unique_kmers = sum(sub_counts) + sum(inter_counts)

    boundaries = [0.25, 0.75]
    boundary_indices = [max(int(percent * num_genomes_in_dataset), 1) for percent in boundaries]
    print(boundary_indices)
    # Special cases where indices are customized ...
    if across_group_analysis:
        boundary_indices = [5, 20]

    metrics[0] = round(sub_counts[0]/total_unique_kmers, 3)
    metrics[1] = round(sum([inter_counts[i] for i in range(1, boundary_indices[0])])/total_unique_kmers, 3)
    metrics[2] = round(sum([inter_counts[i] for i in range(boundary_indices[0], boundary_indices[1])])/total_unique_kmers, 3)
    metrics[3] = round(sum([inter_counts[i] for i in range(boundary_indices[1], len(inter_counts))])/total_unique_kmers, 3)

    metrics[4] = (1 * sub_counts[0]/total_unique_kmers)
    metrics[4] += sum([((i+1) * (inter_counts[i]/total_unique_kmers)) for i in range(1, len(inter_counts))])
    metrics[4] = round(metrics[4], 4)

    rounding_error = abs(sum(metrics[0:4])-1) 
    assert rounding_error < 0.05, "Issue occurred with histogram summarization"
    return metrics


####################################################
# Section 3: Rules needed for this experiment type
####################################################

#   Section 3.1: Builds KMC databases, one rule
#   for pivot and one for all other genomes.

rule build_kmc_database_on_genome_exp_type_2:
    input:
        "input_data/rest_of_set/dataset_{num}/{genome}.fna.gz"
    output:
        temp("step_1/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_pre"),
        temp("step_1/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_suf")
    shell:
        "kmc -fm -m64 -k{wildcards.k} -ci1 {input} step_1/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome} tmp/"

rule build_kmc_database_on_pivot_exp_type_2:
    input:
        "input_data/pivot/dataset_{num}/pivot_{num}.fna.gz"
    output:
        temp("step_1/pivot/k_{k}/dataset_{num}/pivot_{num}.kmc_pre"),
        temp("step_1/pivot/k_{k}/dataset_{num}/pivot_{num}.kmc_suf")
    shell:
        "kmc -fm -m64 -k{wildcards.k} -ci1 {input} step_1/pivot/k_{wildcards.k}/dataset_{wildcards.num}/pivot_{wildcards.num} tmp/"


# Section 3.2: Converts each kmer database into 
# a set (multiplicity=1) for all databases.

rule transform_genome_to_set_exp_type2:
    input:
        "step_1/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_pre",
        "step_1/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_suf"
    output:
        "genome_sets/rest_of_set/k_{k}/dataset_{num}/{genome}.transformed.kmc_pre",
        "genome_sets/rest_of_set/k_{k}/dataset_{num}/{genome}.transformed.kmc_suf"
    shell:
        """
        kmc_tools transform step_1/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome} \
        set_counts 1 genome_sets/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome}.transformed
        """

rule transform_pivot_to_set_exp_type2:
    input:
        "step_1/pivot/k_{k}/dataset_{num}/pivot_{num}.kmc_pre",
        "step_1/pivot/k_{k}/dataset_{num}/pivot_{num}.kmc_suf"
    output:
        "genome_sets/pivot/k_{k}/dataset_{num}/pivot_{num}.transformed.kmc_pre",
        "genome_sets/pivot/k_{k}/dataset_{num}/pivot_{num}.transformed.kmc_suf"
    shell:
        """
        kmc_tools transform step_1/pivot/k_{wildcards.k}/dataset_{wildcards.num}/pivot_{wildcards.num} \
        set_counts 1 genome_sets/pivot/k_{wildcards.k}/dataset_{wildcards.num}/pivot_{wildcards.num}.transformed
        """

# Section 3.3: Takes all genomes within a dataset
# OTHER THAN the pivot, and then unions them together.

rule within_group_union_exp_type2:
    input:
        get_all_genomes_in_dataset
    output:
        "within_databases/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "within_databases/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    shell:
        "kmc_tools complex complex_ops/within_groups/k_{wildcards.k}/dataset_{wildcards.num}/within_dataset_{wildcards.num}.txt"


# Section 3.4: Performs the needed operations to generate
# the plots for the within_group plots.

rule pivot_intersect_within_group_exp_type2:
    input:
        "within_databases/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "genome_sets/pivot/k_{k}/dataset_{num}/pivot_{num}.transformed.kmc_pre"
    output:
        "within_dataset_results/k_{k}/dataset_{num}/intersect/dataset_{num}_pivot_intersect_group.kmc_pre",
        "within_dataset_results/k_{k}/dataset_{num}/intersect/dataset_{num}_pivot_intersect_group.kmc_suf"
    shell:
        """ 
        kmc_tools simple genome_sets/pivot/k_{wildcards.k}/dataset_{wildcards.num}/pivot_{wildcards.num}.transformed \
        within_databases/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined intersect  \
        within_dataset_results/k_{wildcards.k}/dataset_{wildcards.num}/intersect/dataset_{wildcards.num}_pivot_intersect_group -ocsum
        """

rule pivot_subtract_within_group_exp_type2:
    input:
        "within_databases/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "genome_sets/pivot/k_{k}/dataset_{num}/pivot_{num}.transformed.kmc_pre"
    output:
        "within_dataset_results/k_{k}/dataset_{num}/subtract/dataset_{num}_pivot_subtract_group.kmc_pre",
        "within_dataset_results/k_{k}/dataset_{num}/subtract/dataset_{num}_pivot_subtract_group.kmc_suf"
    shell:
        """ 
        kmc_tools simple genome_sets/pivot/k_{wildcards.k}/dataset_{wildcards.num}/pivot_{wildcards.num}.transformed \
        within_databases/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined kmers_subtract \
        within_dataset_results/k_{wildcards.k}/dataset_{wildcards.num}/subtract/dataset_{wildcards.num}_pivot_subtract_group 
        """

rule within_group_histogram_exp_type2:
    input:
        "within_dataset_results/k_{k}/dataset_{num}/{op}/dataset_{num}_pivot_{op}_group.kmc_pre",
        "within_dataset_results/k_{k}/dataset_{num}/{op}/dataset_{num}_pivot_{op}_group.kmc_suf"
    output:
        "within_dataset_results/k_{k}/dataset_{num}/{op}/dataset_{num}_pivot_{op}_group.hist.txt"
    shell:
        """
        kmc_tools transform \
        within_dataset_results/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.op}/dataset_{wildcards.num}_pivot_{wildcards.op}_group \
        histogram within_dataset_results/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.op}/dataset_{wildcards.num}_pivot_{wildcards.op}_group.hist.txt \
        """


# Section 3.5: Performs the within-group analysis ... 

rule within_group_analysis_exp_type2:
    input:
        get_within_group_histogram_files
    output:
        "within_dataset_analysis/within_dataset_analysis.csv"
    run:
        with open(output[0], "w") as output_fd:
            output_fd.write(f"group_num,k,percent_1_occ,percent_25_or_less,percent_25_to_75,percent_75_or_more,unique_stat\n")
            for i in range(0, len(input), 2):
                dataset_num = input[i].split("/")[2].split("_")[1]
                k_value = input[i].split("/")[1].split("_")[1]

                num_genomes_in_dataset = get_num_of_dataset_members(dataset_num)

                with open(input[i], "r") as sub_fd, open(input[i+1], "r") as inter_fd:
                    sub_results = sub_fd.readlines()
                    sub_counts = [int(record.split()[1]) for record in sub_results]

                    inter_results = inter_fd.readlines()
                    inter_counts = [int(record.split()[1]) for record in inter_results]
                
                metrics = summarize_histogram(sub_counts, inter_counts, num_genomes_in_dataset, False)
                metrics_str = ",".join([str(x) for x in metrics])

                output_fd.write(f"group_{dataset_num},{k_value},{metrics_str}\n")


