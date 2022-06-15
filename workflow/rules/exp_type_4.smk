####################################################
# Name: exp_type_4.smk
# Description: Contains functions and rules for
#              the type of experiment type 3: 
# Date: 6/14/22
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

if exp_type == 4:
    # Make tmp directory for kmc
    if not os.path.isdir("tmp"):
        os.mkdir("tmp")

    # Parse out a pivot from the downloaded dataset
    if not os.path.isdir("input_type4/"):
        for i in range(1, num_datasets+1):
            dir_prefix = f"data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"input_type4/rest_of_set/dataset_{i}"):
                os.makedirs(f"input_type4/rest_of_set/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"input_type4/rest_of_set/dataset_{i}")
            
            # Extract the pivot and move to new directory
            dir_prefix = f"input_type4/rest_of_set/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            pivot = random.choice(list_of_files)
            if not os.path.isdir(f"input_type4/pivot/dataset_{i}"):
                os.makedirs(f"input_type4/pivot/dataset_{i}")
            
            # Remove pivot from other set of genomes
            shutil.copy(pivot, f"input_type4/pivot/dataset_{i}/pivot_{i}.fna.gz")
            os.remove(pivot)
    
    # Initialize complex ops directory to start writing files
    if not os.path.isdir("complex_ops"):
        os.makedirs("complex_ops")

    # Build the complex ops files *
    for k in k_values:
        if not os.path.isdir(f"complex_ops/k_{k}"):
            os.mkdir(f"complex_ops/k_{k}")
        
        for num in range(1, num_datasets+1):
            kmc_input_files = []

            for data_file in os.listdir(f"input_type4/rest_of_set/dataset_{num}"):
                if data_file.endswith(".fna.gz"):
                    base_name = data_file.split(".fna.gz")[0]
                    kmc_input_files.append(f"genome_sets/rest_of_set/k_{k}/dataset_{num}/{base_name}.transformed")
            
            if not os.path.isdir(f"complex_ops/k_{k}/dataset_{num}/"):
                os.mkdir(f"complex_ops/k_{k}/dataset_{num}/")

            with open(f"complex_ops/k_{k}/dataset_{num}/ops_{num}.txt", "w") as fd:
                fd.write("INPUT:\n")
                result_str = "("
                for i, path in enumerate(kmc_input_files):
                    fd.write(f"set{i+1} = {path}\n")
                    result_str += "set{} + ".format(i+1)
                result_str = result_str[:-2] + ")"
                fd.write("OUTPUT:\n")
                fd.write(f"unions/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined = {result_str}\n")
                fd.write("OUTPUT_PARAMS:\n-cs5000\n")
        
            

####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_all_genomes_in_dataset(wildcards):
    """ Returns a list of database in a certain dataset """
    input_files = []
    for data_file in os.listdir(f"input_type4/rest_of_set/dataset_{wildcards.num}/"):
        if data_file.endswith(".fna.gz"):
            file_name = data_file.split(".fna.gz")[0]
            input_files.append(f"genome_sets/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{file_name}.transformed.kmc_pre")
    return input_files

####################################################
# Section 3: Rules needed for this experiment type
####################################################

#   Section 3.1: Builds KMC databases, one rule
#   for pivot and one for all other genomes.

rule build_kmc_database_on_genome_exp_type_4:
    input:
        "input_type4/rest_of_set/dataset_{num}/{genome}.fna.gz"
    output:
        "step_1_type4/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_pre",
        "step_1_type4/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_suf"
    shell:
        "kmc -fm -m64 -k{wildcards.k} -ci1 {input} step_1_type4/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome} tmp/"

rule build_kmc_database_on_pivot_exp_type_4:
    input:
        "input_type4/pivot/dataset_{num}/pivot_{num}.fna.gz"
    output:
        "step_1_type4/pivot/k_{k}/dataset_{num}/pivot_{num}.kmc_pre",
        "step_1_type4/pivot/k_{k}/dataset_{num}/pivot_{num}.kmc_suf"
    shell:
        "kmc -fm -m64 -k{wildcards.k} -ci1 {input} step_1_type4/pivot/k_{wildcards.k}/dataset_{wildcards.num}/pivot_{wildcards.num} tmp/"

# Section 3.2: Converts each kmer database into 
# a set (multiplicity=1) for all databases.

rule transform_genome_to_set_exp_type_4:
    input:
        "step_1_type4/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_pre",
        "step_1_type4/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_suf"
    output:
        "genome_sets/rest_of_set/k_{k}/dataset_{num}/{genome}.transformed.kmc_pre",
        "genome_sets/rest_of_set/k_{k}/dataset_{num}/{genome}.transformed.kmc_suf"
    shell:
        """
        kmc_tools transform step_1_type4/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome} \
        set_counts 1 genome_sets/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome}.transformed
        """

# Section 3.3: Takes all genomes within a dataset
# OTHER THAN the pivot, and then unions them together.
# Creates histogram

rule rest_of_set_union_exp_type_4:
    input:
        get_all_genomes_in_dataset
    output:
        "unions/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "unions/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    shell:
        "kmc_tools complex complex_ops/k_{wildcards.k}/dataset_{wildcards.num}/ops_{wildcards.num}.txt"

rule union_histogram_exp_type_4:
    input:
        "unions/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "unions/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    output:
        "unions/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.hist.txt"
    shell:
        """
        kmc_tools transform \
        unions/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined \
        histogram unions/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.hist.txt \
        """

rule transform_union_to_set_exp_type4:
    input: 
        "unions/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "unions/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    output:
        "genome_sets/unions/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.transformed.kmc_pre",
        "genome_sets/unions/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.transformed.kmc_suf"
    shell:
     """
        kmc_tools transform unions/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined \
        set_counts 1 genome_sets/unions/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined.transformed
     """

# Section 3.4: Performs intersection
rule pivot_intersect_exp_type4:
    input:
        "step_1_type4/pivot/k_{k}/dataset_{pivot_num}/pivot_{pivot_num}.kmc_pre",
        "genome_sets/unions/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.transformed.kmc_pre"
    output:
        "intersection_results/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_pre",
        "intersection_results/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_suf"
    shell:
        """
        kmc_tools simple genome_sets/unions/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined.transformed \
        step_1_type4/pivot/k_{wildcards.k}/dataset_{wildcards.pivot_num}/pivot_{wildcards.pivot_num} intersect \
        intersection_results/k_{wildcards.k}/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num} -ocsum
        """

rule intersection_histogram_exp_type4:
    input:
        "intersection_results/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_pre",
        "intersection_results/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_suf"
    output:
        "intersection_results/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.hist.txt"
    shell:
        """
        kmc_tools transform \
        intersection_results/k_{wildcards.k}/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num} \
        histogram intersection_results/k_{wildcards.k}/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num}.hist.txt\
        """

# Section 3.5 Generates text dump of pivot and intersection

rule pivot_text_dump_exp_type4:
    input:
        "step_1_type4/pivot/k_{k}/dataset_{num}/pivot_{num}.kmc_pre",
        "step_1_type4/pivot/k_{k}/dataset_{num}/pivot_{num}.kmc_suf"
    output:
        "text_dump/k_{k}/pivot/pivot_{num}.txt"
    shell: 
        """
        kmc_tools transform \
        step_1_type4/pivot/k_{wildcards.k}/dataset_{wildcards.num}/pivot_{wildcards.num} \
        dump -s text_dump/k_{wildcards.k}/pivot/pivot_{wildcards.num}.txt \
        """

rule intersection_text_dump_exp_type4:
    input:
        "intersection_results/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_pre",
        "intersection_results/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_suf"
    output:
        "text_dump/k_{k}/intersection/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.txt"
    shell:
        """
        kmc_tools transform \
        intersection_results/k_{wildcards.k}/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num} \
        dump -s text_dump/k_{wildcards.k}/intersection/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num}.txt
        """