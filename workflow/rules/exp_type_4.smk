####################################################
# Name: exp_type_4.smk
# Description: Contains functions and rules for
#              the type of experiment type 4: 
#
#              Selects an out-pivot genome from each species, 
#              determines what groups each k-mer occurs in 
#              and builds a confusion matrix.
#
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
            dir_prefix = f"{database_root}/trial_{curr_trial}/exp0_nonpivot_genomes/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"input_type4/rest_of_set/dataset_{i}"):
                os.makedirs(f"input_type4/rest_of_set/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"input_type4/rest_of_set/dataset_{i}")
            
            # Take the pivot and copy to needed folder
            pivot = f"{database_root}/trial_{curr_trial}/exp0_pivot_genomes/dataset_{i}/pivot_{i}.fna.gz"
            if not os.path.isdir(f"input_type4/pivot/"):
                os.makedirs(f"input_type4/pivot/")
            shutil.copy(pivot, f"input_type4/pivot/pivot_{i}.fna.gz")

            # Added pivot to rest_of_set if doing in-pivot
            if not out_pivot_exp4:
                shutil.copy(pivot, f"input_type4/rest_of_set/dataset_{i}")
    
    # Initialize complex ops directory to start writing files
    if not os.path.isdir("complex_ops_type_4"):
        os.makedirs("complex_ops_type_4")

    # Build the complex ops files
    for k in k_values:
        if not os.path.isdir(f"complex_ops_type_4/k_{k}"):
            os.mkdir(f"complex_ops_type_4/k_{k}")
        
        for num in range(1, num_datasets+1):
            kmc_input_files = []

            for data_file in os.listdir(f"input_type4/rest_of_set/dataset_{num}"):
                if data_file.endswith(".fna.gz"):
                    base_name = data_file.split(".fna.gz")[0]
                    kmc_input_files.append(f"genome_sets_type_4/rest_of_set/k_{k}/dataset_{num}/{base_name}.transformed")
            
            if not os.path.isdir(f"complex_ops_type_4/k_{k}/dataset_{num}/"):
                os.mkdir(f"complex_ops_type_4/k_{k}/dataset_{num}/")

            with open(f"complex_ops_type_4/k_{k}/dataset_{num}/ops_{num}.txt", "w") as fd:
                fd.write("INPUT:\n")
                result_str = "("
                for i, path in enumerate(kmc_input_files):
                    fd.write(f"set{i+1} = {path}\n")
                    result_str += "set{} + ".format(i+1)
                result_str = result_str[:-2] + ")"
                fd.write("OUTPUT:\n")
                fd.write(f"unions_type_4/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined = {result_str}\n")
                fd.write("OUTPUT_PARAMS:\n-cs5000\n")

    # Create filelists
    for k in k_values:
        if not os.path.isdir(f"filelists_type_4/k_{k}"):
            os.makedirs(f"filelists_type_4/k_{k}")
        intersect_files = []
        pivot_files = []
        for pivot in range(1, num_datasets+1):
            pivot_files.append("/text_dump_type_4/k_{k}/pivot/pivot_{pivot}.txt".format(k = k, pivot = pivot))
            for num in range(1, num_datasets+1):
                intersect_files.append("/text_dump_type_4/k_{k}/intersection/pivot_{pivot}/pivot_{pivot}_intersect_dataset_{num}.txt".format(k = k, pivot = pivot, num = num))
        with open(f"filelists_type_4/k_{k}/pivots_filelist.txt", "w") as fd:
            for f in pivot_files:
                fd.write(base_dir+f+"\n")
        with open(f"filelists_type_4/k_{k}/intersections_filelist.txt", "w") as fd:
            for f in intersect_files:
                fd.write(base_dir+f+"\n")
    

####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_all_genomes_in_dataset_exp_4(wildcards):
    """ Returns a list of database in a certain dataset """
    input_files = []
    if(exp_type ==4):
        for data_file in os.listdir(f"input_type4/rest_of_set/dataset_{wildcards.num}/"):
            if data_file.endswith(".fna.gz"):
                file_name = data_file.split(".fna.gz")[0]
                input_files.append(f"genome_sets_type_4/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{file_name}.transformed.kmc_pre")
    return input_files

def get_filelists_and_text_dumps_exp_4(wildcards):
    """ Generates text dumps needed for filelists, then returns filelists """
    needed_files = []
    needed_files.append("filelists_type_4/k_{k}/intersections_filelist.txt")
    needed_files.append("filelists_type_4/k_{k}/pivots_filelist.txt")
    for pivot in range(1, num_datasets+1):
            needed_files.append("text_dump_type_4/k_{k}/pivot/pivot_{pivot}.txt".format(k = wildcards.k, pivot = pivot))
            for num in range(1, num_datasets+1):
                needed_files.append("text_dump_type_4/k_{k}/intersection/pivot_{pivot}/pivot_{pivot}_intersect_dataset_{num}.txt".format(k = wildcards.k, pivot = pivot, num = num))
    return needed_files

####################################################
# Section 3: Rules needed for this experiment type
####################################################

#   Section 3.1: Builds KMC databases, one rule
#   for pivot and one for all other genomes.

rule build_kmc_database_on_genome_exp_type_4:
    input:
        "input_type4/rest_of_set/dataset_{num}/{genome}.fna.gz"
    output:
        "step_1_type_4/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_pre",
        "step_1_type_4/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_suf"
    shell:
        "kmc -fm -m64 -k{wildcards.k} -ci1 {input} step_1_type_4/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome} tmp/"

rule build_kmc_database_on_pivot_exp_type_4:
    input:
        "input_type4/pivot/pivot_{num}.fna.gz"
    output:
        "step_1_type_4/pivot/k_{k}/dataset_{num}/pivot_{num}.kmc_pre",
        "step_1_type_4/pivot/k_{k}/dataset_{num}/pivot_{num}.kmc_suf"
    shell:
        "kmc -fm -m64 -k{wildcards.k} -ci1 {input} step_1_type_4/pivot/k_{wildcards.k}/dataset_{wildcards.num}/pivot_{wildcards.num} tmp/"

#   Section 3.2: Converts each kmer database into 
#   a set (multiplicity=1) for all databases._

rule transform_genome_to_set_exp_type_4:
    input:
        "step_1_type_4/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_pre",
        "step_1_type_4/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_suf"
    output:
        "genome_sets_type_4/rest_of_set/k_{k}/dataset_{num}/{genome}.transformed.kmc_pre",
        "genome_sets_type_4/rest_of_set/k_{k}/dataset_{num}/{genome}.transformed.kmc_suf"
    shell:
        """
        kmc_tools transform step_1_type_4/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome} \
        set_counts 1 genome_sets_type_4/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome}.transformed

        # Remove step_1 rest of set
        rm {input[0]} {input[1]}
        """

#   Section 3.3: Takes all genomes within a dataset
#   OTHER THAN the pivot, and then unions them together.

rule rest_of_set_union_exp_type_4:
    input:
        get_all_genomes_in_dataset_exp_4
    output:
        "unions_type_4/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "unions_type_4/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    shell:
        "kmc_tools complex complex_ops_type_4/k_{wildcards.k}/dataset_{wildcards.num}/ops_{wildcards.num}.txt"

rule union_histogram_exp_type_4:
    input:
        "unions_type_4/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "unions_type_4/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    output:
        "unions_type_4/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.hist.txt"
    shell:
        """
        kmc_tools transform \
        unions_type_4/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined \
        histogram unions_type_4/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.hist.txt \
        """

rule transform_union_to_set_exp_type_4:
    input: 
        "unions_type_4/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "unions_type_4/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    output:
        "genome_sets_type_4/unions_type_4/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.transformed.kmc_pre",
        "genome_sets_type_4/unions_type_4/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.transformed.kmc_suf"
    shell:
     """
        kmc_tools transform unions_type_4/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined \
        set_counts 1 genome_sets_type_4/unions_type_4/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined.transformed

        # Remove unions rest of set
        rm {input[0]} {input[1]}
     """

#   Section 3.4: Performs intersection
#   between pivot and unioned databases.

rule pivot_intersect_exp_type_4:
    input:
        "step_1_type_4/pivot/k_{k}/dataset_{pivot_num}/pivot_{pivot_num}.kmc_pre",
        "genome_sets_type_4/unions_type_4/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.transformed.kmc_pre"
    output:
        "intersection_results_type_4/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_pre",
        "intersection_results_type_4/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_suf"
    shell:
        """
        kmc_tools simple genome_sets_type_4/unions_type_4/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined.transformed \
        step_1_type_4/pivot/k_{wildcards.k}/dataset_{wildcards.pivot_num}/pivot_{wildcards.pivot_num} intersect \
        intersection_results_type_4/k_{wildcards.k}/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num} -ocsum
        """

rule intersection_histogram_exp_type_4:
    input:
        "intersection_results_type_4/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_pre",
        "intersection_results_type_4/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_suf"
    output:
        "intersection_results_type_4/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.hist.txt"
    shell:
        """
        kmc_tools transform \
        intersection_results_type_4/k_{wildcards.k}/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num} \
        histogram intersection_results_type_4/k_{wildcards.k}/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num}.hist.txt\
        """

#   Section 3.5 Generates text dump of 
#   pivot and intersection kmc files.

rule pivot_text_dump_exp_type_4:
    input:
        "step_1_type_4/pivot/k_{k}/dataset_{num}/pivot_{num}.kmc_pre",
        "step_1_type_4/pivot/k_{k}/dataset_{num}/pivot_{num}.kmc_suf"
    output:
        "text_dump_type_4/k_{k}/pivot/pivot_{num}.txt"
    shell: 
        """
        kmc_tools transform \
        step_1_type_4/pivot/k_{wildcards.k}/dataset_{wildcards.num}/pivot_{wildcards.num} \
        dump -s text_dump_type_4/k_{wildcards.k}/pivot/pivot_{wildcards.num}.txt \
        """

rule intersection_text_dump_exp_type_4:
    input:
        "intersection_results_type_4/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_pre",
        "intersection_results_type_4/k_{k}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_suf"
    output:
        "text_dump_type_4/k_{k}/intersection/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.txt"
    shell:
        """
        kmc_tools transform \
        intersection_results_type_4/k_{wildcards.k}/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num} \
        dump -s text_dump_type_4/k_{wildcards.k}/intersection/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num}.txt
        """

#   Section 3.7 Runs python script to generate a 
#   confusion matrix and accuracy scores for one value of k.

rule run_merge_list_exp_type_4:
    input:
        get_filelists_and_text_dumps_exp_4
    output:
        "accuracies_type_4/values/k_{k}_accuracy_values.csv",
        "accuracies_type_4/confusion_matrix/k_{k}_confusion_matrix.txt"
    shell:
        """
        python3 {repo_dir}/src/merge_lists.py \
        -p {base_dir}/filelists_type_4/k_{wildcards.k}/pivots_filelist.txt \
        -i {base_dir}/filelists_type_4/k_{wildcards.k}/intersections_filelist.txt \
        -o {base_dir}/accuracies_type_4/ \
        -n {num_datasets} \
        -k {wildcards.k} \

        # Remove text dumps for this k
        rm text_dump_type_4/k_{wildcards.k}/intersection/pivot_*/pivot_*_intersect_dataset_*.txt
        rm text_dump_type_4/k_{wildcards.k}/pivot/pivot_*.txt
        """

#   Section 3.8 Concatenates accuracy score tables
#   for all values of k to final csv file.

rule concatenate_accuracies_exp_type4:
    input:
        expand("accuracies_type_4/values/k_{k}_accuracy_values.csv", k = k_values)
    output:
        "accuracies_type_4/accuracy_values.csv"
    shell:
        "cat accuracies_type_4/values/*.csv > accuracies_type_4/accuracy_values.csv"


#   Section 3.9: Overall rule for Experiment 4 to generate all the 
#   output files.

rule generate_exp4_output:
    input:
        "accuracies_type_4/accuracy_values.csv"