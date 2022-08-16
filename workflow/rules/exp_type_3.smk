####################################################
# Name: exp_type_3.smk
# Description: Contains functions and rules for
#              the type of experiment type 3: 
#
#              Simulates reads from an out-pivot genome
#              and measures what percentage of each read's 
#              k-mers are found in databases of different species 
#              across values of k.
#
# Date: 2/10/22
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

if exp_type == 3:
    # Make tmp directory for kmc
    if not os.path.isdir("tmp"):
        os.mkdir("tmp")
    
    # Parse out a pivot from the downloaded dataset
    if not os.path.isdir("input_type3/"):
        for i in range(1, num_datasets+1):
            dir_prefix = f"data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"input_type3/rest_of_set/dataset_{i}"):
                os.makedirs(f"input_type3/rest_of_set/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"input_type3/rest_of_set/dataset_{i}")
            
            # Extract the pivot and move to new directory
            dir_prefix = f"input_type3/rest_of_set/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            pivot = random.choice(list_of_files)
            if not os.path.isdir(f"input_type3/pivot/dataset_{i}"):
                os.makedirs(f"input_type3/pivot/dataset_{i}")
            
            # Remove pivot from other set of genomes
            shutil.copy(pivot, f"input_type3/pivot/dataset_{i}/pivot_{i}.fna.gz")
            os.remove(pivot)

    # Initialize complex ops directory to start writing files
    if not os.path.isdir("complex_ops"):
        os.makedirs("complex_ops/within_groups")

    # Build the complex ops files for within groups
    for k in k_values:
        if not os.path.isdir(f"complex_ops/within_groups/k_{k}"):
            os.mkdir(f"complex_ops/within_groups/k_{k}")
        
        for num in range(1, num_datasets+1):
            kmc_input_files = []

            for data_file in os.listdir(f"input_type3/rest_of_set/dataset_{num}"):
                if data_file.endswith(".fna.gz"):
                    base_name = data_file.split(".fna.gz")[0]
                    kmc_input_files.append(f"genome_sets_type3/rest_of_set/k_{k}/dataset_{num}/{base_name}.transformed")
            
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
                fd.write(f"within_databases_type3/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined = {result_str}\n")
                fd.write("OUTPUT_PARAMS:\n-cs5000\n")

####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_all_genomes_in_dataset_type3(wildcards):
    """ Returns a list of database in a certain dataset """
    input_files = []
    for data_file in os.listdir(f"input_type3/rest_of_set/dataset_{wildcards.num}/"):
        if data_file.endswith(".fna.gz"):
            file_name = data_file.split(".fna.gz")[0]
            input_files.append(f"genome_sets_type3/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{file_name}.transformed.kmc_pre")
    return input_files

def get_all_histogram_files(wildcards):
    """ Returns a list of histogram files needed for final analysis """
    input_files = []
    for read_type in ["illumina", "ont"]:
        for pivot_num in range(1, num_datasets+1):
            for k in k_values:
                # Add the pivot histogram file
                input_files.append(f"genome_sets_type3/pivot/{read_type}/k_{k}/dataset_{pivot_num}/pivot_{pivot_num}.transformed.hist.txt")
                for num in range(1, num_datasets+1):
                    input_files.append(f"within_dataset_results_type3/{read_type}/pivot_{pivot_num}/k_{k}/dataset_{num}/intersect/dataset_{num}_pivot_intersect_group.hist.txt")
    return input_files

####################################################
# Section 3: Rules needed for this experiment type
####################################################

#   Section 3.1: Simulate short, and long reads for 
#   every pivot genome

rule generate_raw_positive_short_reads:
    input:
        "input_type3/pivot/dataset_{num}/pivot_{num}.fna.gz"
    output:
        temp("input_type3/pivot_raw_reads/illumina/dataset_{num}/pivot_{num}_illumina_reads.fq")
    shell:
        """
        gzip -d {input} -k -f
        art_illumina -ss HS25 -i input_type3/pivot/dataset_{wildcards.num}/pivot_{wildcards.num}.fna -na -l 150 -f 2.0 \
        -o input_type3/pivot_raw_reads/illumina/dataset_{wildcards.num}/pivot_{wildcards.num}_illumina_reads
        """

rule generate_raw_positive_long_reads:
    input:
        "input_type3/pivot/dataset_{num}/pivot_{num}.fna.gz"
    output:
        temp("input_type3/pivot_raw_reads/ont/dataset_{num}/pivot_{num}_ont_reads.fastq")
    run:
        shell("gzip -d {input} -k -f")
        shell("""
        pbsim --depth 30.0 --prefix input_type3/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}_ont_reads \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 200 input_type3/pivot/dataset_{wildcards.num}/pivot_{wildcards.num}.fna
        """)
        shell("cat 'input_type3/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}_ont_reads'*.fastq > {output}")
        shell("ls  'input_type3/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}_ont_reads_'*.fastq | xargs rm")
        shell("ls  'input_type3/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}_ont_reads_'*.ref | xargs rm")
        shell("ls  'input_type3/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}_ont_reads_'*.maf | xargs rm")

rule convert_long_reads_to_fasta_and_subset:
    input:
        "input_type3/pivot_raw_reads/ont/dataset_{num}/pivot_{num}_ont_reads.fastq"
    output:
        "input_type3/pivot_reads/ont/dataset_{num}/pivot_{num}_ont_reads.fa"
    run:
        num_lines = num_reads_per_dataset * 4
        shell("head -n{num_lines} {input} > input_type3/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}_ont_reads_subset.fastq")
        shell("seqtk seq -a input_type3/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}_ont_reads_subset.fastq > {output}")
        shell("if [ $(grep -c '>' {output}) != {num_reads_per_dataset} ]; then echo 'number of reads assertion failed.'; exit 1; fi")
        shell("rm input_type3/pivot_raw_reads/ont/dataset_{wildcards.num}/pivot_{wildcards.num}_ont_reads_subset.fastq")

rule convert_short_reads_to_fasta_and_subset:
    input:
        "input_type3/pivot_raw_reads/illumina/dataset_{num}/pivot_{num}_illumina_reads.fq"
    output:
        "input_type3/pivot_reads/illumina/dataset_{num}/pivot_{num}_illumina_reads.fa"
    run:
        num_lines = num_reads_per_dataset * 4
        shell("head -n{num_lines} {input} > input_type3/pivot_raw_reads/illumina/dataset_{wildcards.num}/pivot_{wildcards.num}_illumina_reads_subset.fq")
        shell("seqtk seq -a input_type3/pivot_raw_reads/illumina/dataset_{wildcards.num}/pivot_{wildcards.num}_illumina_reads_subset.fq > {output}")
        shell("if [ $(grep -c '>' {output}) != {num_reads_per_dataset} ]; then echo 'number of reads assertion failed.'; exit 1; fi")
        shell("rm input_type3/pivot_raw_reads/illumina/dataset_{wildcards.num}/pivot_{wildcards.num}_illumina_reads_subset.fq")

#   Section 3.2: Builds KMC databases, one rule
#   for pivot and one for all other genomes.

rule build_kmc_database_on_genome_exp_type_3:
    input:
        "input_type3/rest_of_set/dataset_{num}/{genome}.fna.gz"
    output:
        temp("step_1_type3/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_pre"),
        temp("step_1_type3/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_suf")
    shell:
        "kmc -fm -m64 -k{wildcards.k} -ci1 {input} step_1_type3/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome} tmp/"

rule build_kmc_database_on_pivot_reads_exp_type_3:
    input:
        "input_type3/pivot_reads/{read_type}/dataset_{num}/pivot_{num}_{read_type}_reads.fa"
    output:
        temp("step_1_type3/pivot/{read_type}/k_{k}/dataset_{num}/pivot_{num}.kmc_pre"),
        temp("step_1_type3/pivot/{read_type}/k_{k}/dataset_{num}/pivot_{num}.kmc_suf")
    shell:
        "kmc -fm -m64 -k{wildcards.k} -ci1 {input} step_1_type3/pivot/{wildcards.read_type}/k_{wildcards.k}/dataset_{wildcards.num}/pivot_{wildcards.num} tmp/"

# Section 3.3: Converts each kmer database into 
# a set (multiplicity=1) for all databases.

rule transform_genome_to_set_exp_type3:
    input:
        "step_1_type3/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_pre",
        "step_1_type3/rest_of_set/k_{k}/dataset_{num}/{genome}.kmc_suf"
    output:
        "genome_sets_type3/rest_of_set/k_{k}/dataset_{num}/{genome}.transformed.kmc_pre",
        "genome_sets_type3/rest_of_set/k_{k}/dataset_{num}/{genome}.transformed.kmc_suf"
    shell:
        """
        kmc_tools transform step_1_type3/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome} \
        set_counts 1 genome_sets_type3/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/{wildcards.genome}.transformed
        """

rule transform_pivot_reads_to_set_exp_type3:
    input:
        "step_1_type3/pivot/{read_type}/k_{k}/dataset_{num}/pivot_{num}.kmc_pre",
        "step_1_type3/pivot/{read_type}/k_{k}/dataset_{num}/pivot_{num}.kmc_suf"
    output:
        "genome_sets_type3/pivot/{read_type}/k_{k}/dataset_{num}/pivot_{num}.transformed.kmc_pre",
        "genome_sets_type3/pivot/{read_type}/k_{k}/dataset_{num}/pivot_{num}.transformed.kmc_suf"
    shell:
        """
        kmc_tools transform step_1_type3/pivot/{wildcards.read_type}/k_{wildcards.k}/dataset_{wildcards.num}/pivot_{wildcards.num} \
        set_counts 1 genome_sets_type3/pivot/{wildcards.read_type}/k_{wildcards.k}/dataset_{wildcards.num}/pivot_{wildcards.num}.transformed
        """

# Section 3.4: Takes all genomes within a dataset
# OTHER THAN the pivot, and then unions them together.

rule within_group_union_exp_type3:
    input:
        get_all_genomes_in_dataset_type3
    output:
        "within_databases_type3/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "within_databases_type3/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_suf"
    shell:
        "kmc_tools complex complex_ops/within_groups/k_{wildcards.k}/dataset_{wildcards.num}/within_dataset_{wildcards.num}.txt"

# Section 3.5: Performs the needed operations to figure out
# what % of kmers in pivot occur in other datasets, and 
# generate the histograms

rule pivot_intersect_within_group_exp_type3:
    input:
        "within_databases_type3/rest_of_set/k_{k}/dataset_{num}/dataset_{num}.transformed.combined.kmc_pre",
        "genome_sets_type3/pivot/{read_type}/k_{k}/dataset_{pivot_num}/pivot_{pivot_num}.transformed.kmc_pre"
    output:
        "within_dataset_results_type3/{read_type}/pivot_{pivot_num}/k_{k}/dataset_{num}/intersect/dataset_{num}_pivot_intersect_group.kmc_pre",
        "within_dataset_results_type3/{read_type}/pivot_{pivot_num}/k_{k}/dataset_{num}/intersect/dataset_{num}_pivot_intersect_group.kmc_suf"
    shell:
        """ 
        kmc_tools simple genome_sets_type3/pivot/{wildcards.read_type}/k_{wildcards.k}/dataset_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}.transformed \
        within_databases_type3/rest_of_set/k_{wildcards.k}/dataset_{wildcards.num}/dataset_{wildcards.num}.transformed.combined intersect  \
        within_dataset_results_type3/{wildcards.read_type}/pivot_{wildcards.pivot_num}/k_{wildcards.k}/dataset_{wildcards.num}/intersect/dataset_{wildcards.num}_pivot_intersect_group -ocsum
        """

rule intersection_histogram_exp_type3:
    input:
        "within_dataset_results_type3/{read_type}/pivot_{pivot_num}/k_{k}/dataset_{num}/intersect/dataset_{num}_pivot_intersect_group.kmc_pre",
        "within_dataset_results_type3/{read_type}/pivot_{pivot_num}/k_{k}/dataset_{num}/intersect/dataset_{num}_pivot_intersect_group.kmc_suf"
    output:
        "within_dataset_results_type3/{read_type}/pivot_{pivot_num}/k_{k}/dataset_{num}/intersect/dataset_{num}_pivot_intersect_group.hist.txt"
    shell:
        """
        kmc_tools transform \
        within_dataset_results_type3/{wildcards.read_type}/pivot_{wildcards.pivot_num}/k_{wildcards.k}/dataset_{wildcards.num}/intersect/dataset_{wildcards.num}_pivot_intersect_group \
        histogram within_dataset_results_type3/{wildcards.read_type}/pivot_{wildcards.pivot_num}/k_{wildcards.k}/dataset_{wildcards.num}/intersect/dataset_{wildcards.num}_pivot_intersect_group.hist.txt \
        """

rule pivot_histogram_exp_type3:
    input:
        "genome_sets_type3/pivot/{read_type}/k_{k}/dataset_{pivot_num}/pivot_{pivot_num}.transformed.kmc_pre",
        "genome_sets_type3/pivot/{read_type}/k_{k}/dataset_{pivot_num}/pivot_{pivot_num}.transformed.kmc_suf"
    output:
        "genome_sets_type3/pivot/{read_type}/k_{k}/dataset_{pivot_num}/pivot_{pivot_num}.transformed.hist.txt"
    shell:
        """
        kmc_tools transform \
        genome_sets_type3/pivot/{wildcards.read_type}/k_{wildcards.k}/dataset_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}.transformed \
        histogram genome_sets_type3/pivot/{wildcards.read_type}/k_{wildcards.k}/dataset_{wildcards.pivot_num}/pivot_{wildcards.pivot_num}.transformed.hist.txt \
        """

# Section 3.6: Performs the analysis

rule intersection_with_group_analysis_exp_type3:
    input:
        get_all_histogram_files
    output:
        "final_analysis_type3/final_analysis_type3.csv"
    run:
        with open(output[0], "w") as output_fd:
            output_fd.write(f"read_type,pivot_num,k,dataset_num,intersection_percent\n")
            pos = 0
            while pos < len(input):
                pivot_file = input[pos]; pos+=1

                # Gather the total number of kmers in 
                with open(pivot_file, "r") as input_fd:
                    pivot_results = input_fd.readlines()
                    pivot_counts = [int(record.split()[1]) for record in pivot_results]

                assert sum(pivot_counts[1:]) == 0, "issue with pivot histogram file"
                num_kmers_in_pivot = pivot_counts[0]

                for curr_pos in range(pos, pos+num_datasets):
                    # Gather current file identification
                    file_name = input[curr_pos]
                    read_type = file_name.split("/")[1]
                    pivot_num = file_name.split("/")[2].split("_")[1]
                    dataset_num = file_name.split("/")[4].split("_")[1]
                    k_value = file_name.split("/")[3].split("_")[1]

                    # Extract number of kmers that intersect with pivot
                    with open(input[curr_pos], "r") as input_fd:
                        data_lines = input_fd.readlines()
                        data_counts = [int(record.split()[1]) for record in data_lines]
                    
                    assert data_counts[0] == 0, "issue with histogram of intersection file"
                    num_kmers_intersect = sum(data_counts)

                    output_metrics = [read_type, pivot_num, k_value, dataset_num, round(num_kmers_intersect/num_kmers_in_pivot, 4)]
                    output_fd.write(",".join([str(x) for x in output_metrics]) + "\n")

                pos += num_datasets