####################################################
# Name: exp_type_6.smk
# Description: Contains functions and rules for
#              the type of experiment type 6: 
#
#              Simulates reads from an out-pivot genome from
#              each species. Then, determines which species each 
#              kmer occurs in (across k) and builds a confusion matrix.
#
# Date: 6/30/22
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

if exp_type == 6:
    # Make tmp directory for kmc
    if not os.path.isdir("tmp"):
        os.mkdir("tmp")

    # Parse out a pivot from the downloaded dataset
    if not os.path.isdir("input_type_6/"):
        for i in range(1, num_datasets+1):
            dir_prefix = f"data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"input_type_6/rest_of_set/dataset_{i}"):
                os.makedirs(f"input_type_6/rest_of_set/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"input_type_6/rest_of_set/dataset_{i}")
            
            # Extract the pivot and move to new directory
            dir_prefix = f"input_type_6/rest_of_set/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            pivot = random.choice(list_of_files)

            # Print pivot name to list in trial folder
            with open('/home/ext-mcheng29/scratch16-blangme2/marie/trial_data/pivot_names.txt', 'a') as fd: # Remove hardcode
                fd.write(pivot+"\n")

            if not os.path.isdir(f"input_type_6/pivot/dataset_{i}"):
                os.makedirs(f"input_type_6/pivot/dataset_{i}")
            
            # Remove pivot from other set of genomes
            shutil.copy(pivot, f"input_type_6/pivot/dataset_{i}/pivot_{i}.fna.gz")
            os.remove(pivot)
    
    # Initialize complex ops directory to start writing files
    if not os.path.isdir("complex_ops_type_6"):
        os.makedirs("complex_ops_type_6")

    # Build the complex ops files
    for k in k_values:
        if not os.path.isdir(f"complex_ops_type_6/k_{k}"):
            os.mkdir(f"complex_ops_type_6/k_{k}")
        
        for num in range(1, num_datasets+1):
            kmc_input_files = []

            for data_file in os.listdir(f"input_type_6/rest_of_set/dataset_{num}"):
                if data_file.endswith(".fna.gz"):
                    base_name = data_file.split(".fna.gz")[0]
                    kmc_input_files.append(f"genome_sets_type_6/rest_of_set/k_{k}/dataset_{num}/{base_name}.transformed")
            
            if not os.path.isdir(f"complex_ops_type_6/k_{k}/dataset_{num}/"):
                os.mkdir(f"complex_ops_type_6/k_{k}/dataset_{num}/")

            with open(f"complex_ops_type_6/k_{k}/dataset_{num}/ops_{num}.txt", "w") as fd:
                fd.write("INPUT:\n")
                result_str = "("
                for i, path in enumerate(kmc_input_files):
                    fd.write(f"set{i+1} = {path}\n")
                    result_str += "set{} + ".format(i+1)
                result_str = result_str[:-2] + ")"
                fd.write("OUTPUT:\n")
                fd.write(f"unions_type_6/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined = {result_str}\n")
                fd.write("OUTPUT_PARAMS:\n-cs5000\n")

    # Create filelists

    read_types = ["illumina","ont"]
    
    for k in k_values:
        if not os.path.isdir(f"filelists_type_6/k_{k}"):
            os.makedirs(f"filelists_type_6/k_{k}")
        for read_type in read_types:
            if not os.path.isdir(f"filelists_type_6/k_{k}/{read_type}"):
                os.makedirs(f"filelists_type_6/k_{k}/{read_type}")
            intersect_files = []
            pivot_files = []
            for pivot in range(1, num_datasets+1):
                pivot_files.append("/text_dump_type_6/k_{k}/{read_type}/pivot/pivot_{pivot}.txt".format(k = k, pivot = pivot, read_type = read_type))
                for num in range(1, num_datasets+1):
                    intersect_files.append("/text_dump_type_6/k_{k}/{read_type}/intersection/pivot_{pivot}/pivot_{pivot}_intersect_dataset_{num}.txt".format(k = k, pivot = pivot, num = num, read_type = read_type))
            with open(f"filelists_type_6/k_{k}/{read_type}/pivots_filelist.txt", "w") as fd:
                for f in pivot_files:
                    fd.write(base_dir+f+"\n")
            with open(f"filelists_type_6/k_{k}/{read_type}/intersections_filelist.txt", "w") as fd:
                for f in intersect_files:
                    fd.write(base_dir+f+"\n")

    

####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_all_genomes_in_dataset_exp_6(wildcards):
    """ Returns a list of database in a certain dataset """
    input_files = []
    for data_file in os.listdir(f"input_type_6/rest_of_set/dataset_{wildcards.num}/"):
        if data_file.endswith(".fna.gz"):
            file_name = data_file.split(".fna.gz")[0]
            input_files.append(f"genome_sets_type_6/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{file_name}.transformed.kmc_pre")
    return input_files

def get_filelists_and_text_dumps_exp_6(wildcards):
    """ Generates text dumps needed for filelists, then returns filelists """
    needed_files = []
    needed_files.append("filelists_type_6/k_{k}/{read_type}/intersections_filelist.txt".format(k = wildcards.k, read_type = wildcards.read_type))
    needed_files.append("filelists_type_6/k_{k}/{read_type}/pivots_filelist.txt".format(k = wildcards.k, read_type = wildcards.read_type))
    for pivot in range(1, num_datasets+1):
            needed_files.append("text_dump_type_6/k_{k}/{read_type}/pivot/pivot_{pivot}.txt".format(k = wildcards.k, pivot = pivot, read_type = wildcards.read_type))
            for num in range(1, num_datasets+1):
                needed_files.append("text_dump_type_6/k_{k}/{read_type}/intersection/pivot_{pivot}/pivot_{pivot}_intersect_dataset_{num}.txt".format(k = wildcards.k, pivot = pivot, num = num, read_type = wildcards.read_type))
    return needed_files

####################################################
# Section 3: Rules needed for this experiment type
####################################################

#   Section 3.1: Simulate short, and long reads for 
#   every pivot genome

rule generate_raw_positive_short_reads_exp_6:
    input:
        "input_type_6/pivot/dataset_{num}/pivot_{num}.fna.gz"
    output:
        temp("input_type_6/pivot_raw_reads/illumina/dataset_{num}/pivot_{num}.fq")
    shell:
        """
        gzip -d {input} -k -f
        art_illumina -ss HS25 -i input_type_6/pivot/dataset_{wildcards.num}/pivot_{wildcards.num}.fna -na -l 150 -f 2.0 \
        -o input_type_6/pivot_raw_reads/illumina/dataset_{wildcards.num}/pivot_{wildcards.num}
        """

rule generate_raw_positive_long_reads_exp_6:
    input:
        "input_type_6/pivot/dataset_{num}/pivot_{num}.fna.gz"
    output:
        "input_type_6/pivot_raw_reads/ont/dataset_{num}/pivot_{num}.fastq"
    shell:
        """
        gzip -d {input} -k -f
        pbsim --depth 30.0 --prefix input_type_6/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num} \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 200 input_type_6/pivot/dataset_{wildcards.num}/pivot_{wildcards.num}.fna
        
        cat 'input_type_6/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}'*.fastq > {output}
        rm 'input_type_6/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}_'*.fastq
        ls  'input_type_6/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}'*.ref | xargs rm
        ls  'input_type_6/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}'*.maf | xargs rm
        """

rule convert_long_reads_to_fasta_and_subset_exp_6:
    input:
        "input_type_6/pivot_raw_reads/ont/dataset_{num}/pivot_{num}.fastq"
    output:
        "input_type_6/pivot_reads/ont/dataset_{num}/pivot_{num}.fa"
    shell:
        """
        num_lines=$(({num_reads_per_dataset} * 4))

        head -n $num_lines {input} > input_type_6/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}_subset.fastq
        seqtk seq -a input_type_6/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}_subset.fastq > {output}
        if [ $(grep -c '>' {output}) != {num_reads_per_dataset} ]; then echo 'number of reads assertion failed.'; exit 1; fi
        rm input_type_6/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}_subset.fastq
        """

rule convert_short_reads_to_fasta_and_subset_exp_6:
    input:
        "input_type_6/pivot_raw_reads/illumina/dataset_{num}/pivot_{num}.fq"
    output:
        "input_type_6/pivot_reads/illumina/dataset_{num}/pivot_{num}.fa"
    shell:
        """
        num_lines=$(({num_reads_per_dataset} * 4))

        head -n $num_lines {input} > input_type_6/pivot_raw_reads/illumina/dataset_{wildcards.num}/pivot_{wildcards.num}_subset.fq
        seqtk seq -a input_type_6/pivot_raw_reads/illumina/dataset_{wildcards.num}/pivot_{wildcards.num}_subset.fq > {output}
        if [ $(grep -c '>' {output}) != {num_reads_per_dataset} ]; then echo 'number of reads assertion failed.'; exit 1; fi
        rm input_type_6/pivot_raw_reads/illumina/dataset_{wildcards.num}/pivot_{wildcards.num}_subset.fq
        """

# Section 3.2: Builds KMC databases, one rule
# for pivot and one for all other genomes.

rule build_kmc_database_on_genome_exp_type_6:
    input:
        "input_type_6/rest_of_set/dataset_{num}/{genome}.fna.gz"
    output:
        "step_1_type_6/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_pre",
        "step_1_type_6/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_suf"
    shell:
        "kmc -fm -m64 -k{wildcards.k} -ci1 {input} step_1_type_6/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome} tmp/"

rule build_kmc_database_on_pivot_exp_type_6:
    input:
        "input_type_6/pivot_reads/{read_type}/dataset_{num}/pivot_{num}.fa"
    output:
        "step_1_type_6/pivot/k_{k}/{read_type}/pivot_{num}.kmc_pre",
        "step_1_type_6/pivot/k_{k}/{read_type}/pivot_{num}.kmc_suf"
    shell:
        "kmc -fm -m64 -k{wildcards.k} -ci1 {input} step_1_type_6/pivot/k_{wildcards.k}/{wildcards.read_type}/pivot_{wildcards.num} tmp/"

# Section 3.3: Converts each kmer database into 
# a set (multiplicity=1) for all databases.

rule transform_genome_to_set_exp_type_6:
    input:
        "step_1_type_6/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_pre",
        "step_1_type_6/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_suf"
    output:
        "genome_sets_type_6/rest_of_set/k_{k}/dataset_{num}/{genome}.transformed.kmc_pre",
        "genome_sets_type_6/rest_of_set/k_{k}/dataset_{num}/{genome}.transformed.kmc_suf"
    shell:
        """
        kmc_tools transform step_1_type_6/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome} \
        set_counts 1 genome_sets_type_6/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome}.transformed
        """

# Section 3.4: Takes all genomes within a dataset
# OTHER THAN the pivot, and then unions them together.

rule rest_of_set_union_exp_type_6:
    input:
        get_all_genomes_in_dataset_exp_6
    output:
        "unions_type_6/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "unions_type_6/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    shell:
        "kmc_tools complex complex_ops_type_6/k_{wildcards.k}/dataset_{wildcards.num}/ops_{wildcards.num}.txt"

# For sanity check

rule union_histogram_exp_type_6:
    input:
        "unions_type_6/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "unions_type_6/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    output:
        "unions_type_6/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.hist.txt"
    shell:
        """
        kmc_tools transform \
        unions_type_6/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined \
        histogram unions_type_6/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.hist.txt \
        """

rule transform_union_to_set_exp_type_6:
    input: 
        "unions_type_6/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "unions_type_6/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    output:
        "genome_sets_type_6/unions/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.transformed.kmc_pre",
        "genome_sets_type_6/unions/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.transformed.kmc_suf"
    shell:
     """
        kmc_tools transform unions_type_6/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined \
        set_counts 1 genome_sets_type_6/unions/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined.transformed
     """

# Section 3.5: Performs intersection
# between pivot and unioned databases.

rule pivot_intersect_exp_type_6:
    input:
        "step_1_type_6/pivot/k_{k}/{read_type}/pivot_{pivot_num}.kmc_pre",
        "genome_sets_type_6/unions/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.transformed.kmc_pre"
    output:
        "intersection_results_type_6/k_{k}/{read_type}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_pre",
        "intersection_results_type_6/k_{k}/{read_type}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_suf"
    shell:
        """
        kmc_tools simple genome_sets_type_6/unions/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined.transformed \
        step_1_type_6/pivot/k_{wildcards.k}/{wildcards.read_type}/pivot_{wildcards.pivot_num} intersect \
        intersection_results_type_6/k_{wildcards.k}/{wildcards.read_type}/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num} -ocsum
        """

rule intersection_histogram_exp_type_6:
    input:
        "intersection_results_type_6/k_{k}/{read_type}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_pre",
        "intersection_results_type_6/k_{k}/{read_type}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_suf"
    output:
        "intersection_results_type_6/k_{k}/{read_type}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.hist.txt"
    shell:
        """
        kmc_tools transform \
        intersection_results_type_6/k_{wildcards.k}/{wildcards.read_type}/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num} \
        histogram intersection_results_type_6/k_{wildcards.k}/{wildcards.read_type}/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num}.hist.txt\
        """

# Section 3.6 Generates text dump of 
# pivot and intersection kmc files.

rule pivot_text_dump_exp_type_6:
    input:
        "step_1_type_6/pivot/k_{k}/{read_type}/pivot_{num}.kmc_pre",
        "step_1_type_6/pivot/k_{k}/{read_type}/pivot_{num}.kmc_suf"
    output:
        "text_dump_type_6/k_{k}/{read_type}/pivot/pivot_{num}.txt"
    shell: 
        """
        kmc_tools transform \
        step_1_type_6/pivot/k_{wildcards.k}/{wildcards.read_type}/pivot_{wildcards.num} \
        dump -s text_dump_type_6/k_{wildcards.k}/{wildcards.read_type}/pivot/pivot_{wildcards.num}.txt \
        """

rule intersection_text_dump_exp_type_6:
    input:
        "intersection_results_type_6/k_{k}/{read_type}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_pre",
        "intersection_results_type_6/k_{k}/{read_type}/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.kmc_suf"
    output:
        "text_dump_type_6/k_{k}/{read_type}/intersection/pivot_{pivot_num}/pivot_{pivot_num}_intersect_dataset_{num}.txt"
    shell:
        """
        kmc_tools transform \
        intersection_results_type_6/k_{wildcards.k}/{wildcards.read_type}/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num} \
        dump -s text_dump_type_6/k_{wildcards.k}/{wildcards.read_type}/intersection/pivot_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}_intersect_dataset_{wildcards.num}.txt
        """

# Section 3.7 Runs python script to generate a 
# confusion matrix and accuracy scores for one value of k.

rule run_merge_list_exp_type_6:
    input:
        get_filelists_and_text_dumps_exp_6
    output:
        "accuracies_type_6/{read_type}/accuracy/k_{k}_accuracy.csv",
        "accuracies_type_6/{read_type}/confusion_matrix/k_{k}_confusion_matrix.csv"
    shell:
        """
        python3 {repo_dir}/src/merge_lists.py \
        -p {base_dir}/filelists_type_6/k_{wildcards.k}/{wildcards.read_type}/pivots_filelist.txt \
        -i {base_dir}/filelists_type_6/k_{wildcards.k}/{wildcards.read_type}/intersections_filelist.txt \
        -o {base_dir}/accuracies_type_6/{wildcards.read_type}/ \
        -n {num_datasets} \
        -k {wildcards.k} \
        """


# Section 3.8 Concatenates accuracy score tables
# for all values of k to final csv file.

rule concatenate_accuracies_exp_type_6:
    input:
        expand("accuracies_type_6/{read_type}/accuracy/k_{k}_accuracy.csv", k = k_values, read_type = {"illumina","ont"})
    output:
        "accuracies_type_6/short_accuracy_scores.csv",
        "accuracies_type_6/long_accuracy_scores.csv"
    shell:
        """
        cat accuracies_type_6/illumina/accuracy/k_*_accuracy.csv > accuracies_type_6/short_accuracy_scores.csv
        cat accuracies_type_6/ont/accuracy/k_*_accuracy.csv > accuracies_type_6/long_accuracy_scores.csv
        """
        
