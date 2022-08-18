####################################################
# Name: exp_type_7.smk
# Description: Contains functions and rules for
#              the type of experiment type 7:
# 
#              Simulates reads from an out-pivot genome from
#              each species. Then, determines which species each 
#              half-mem/mem occurs in and builds a confusion matrix.
#
# Date: 7/5/22
####################################################

####################################################
# Section 1: Code that is always run, when this 
#            experiment is called. It genereates the
#            needed files.
####################################################

if exp_type == 7:
    if not os.path.isdir("exp7_non_pivot_data/"):
        for i in range(1, num_datasets+1):
            dir_prefix = f"data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"exp7_non_pivot_data/dataset_{i}"):
                os.makedirs(f"exp7_non_pivot_data/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"exp7_non_pivot_data/dataset_{i}")
            
            # Extract the pivot and move to new directory
            dir_prefix = f"exp7_non_pivot_data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            pivot = random.choice(list_of_files)

            # Print pivot name to list in trial folder
            with open(trial_info_dir + "/exp7_pivot_names.txt", 'a') as fd:
                fd.write(pivot+"\n")

            if not os.path.isdir(f"exp7_pivot_data/pivot_ref/"):
                os.makedirs(f"exp7_pivot_data/pivot_ref/")
            
            # Remove pivot from other set of genomes
            shutil.copy(pivot, f"exp7_pivot_data/pivot_ref/pivot_{i}.fna")
            os.remove(pivot)

####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_dataset_non_pivot_genomes_exp7(wildcards):
    """ Returns list of non-pivot genomes for a given dataset """
    input_files = []
    #if exp_type == 7:
    for data_file in os.listdir(f"exp7_non_pivot_data/dataset_{wildcards.num}/"):
        if data_file.endswith(".fna"):
            file_name = data_file.split(".fna")[0]
            input_files.append(f"exp7_non_pivot_data/dataset_{wildcards.num}/{file_name}.fna")
    return input_files

def get_all_non_pivot_genomes_exp7(wildcards):
    """ Returns list of all non-pivot genomes """
    input_files = []
    #if exp_type == 7:  
    for i in range(1, num_datasets + 1):
        for data_file in os.listdir(f"exp7_non_pivot_data/dataset_{i}/"):
            if data_file.endswith(".fna"):
                file_name = data_file.split(".fna")[0]
                input_files.append(f"exp7_non_pivot_data/dataset_{i}/{file_name}.fna")
    return input_files

def get_combined_ref_dataset_exp7(wildcards):
    """ Returns list of forward combined refs for each dataset """
    input_files = []
    #if exp_type == 7:
    for i in range(1, num_datasets + 1):
        input_files.append(f"exp7_combined_refs/dataset_{i}/combined_ref_forward.fna")
    return input_files

def get_python_input_exp7(wildcards):
    """ Returns a list of all pivot SAM files for a read type """
    input_files = []
    #if(exp_type == 7):
    for i in range(1,num_datasets+1):
        #input_files.append(f"ref_lists_type_7/dataset_{i}_references.txt")
        for j in range(1, num_datasets+1):
            input_files.append(f"exp7_sam_files/{wildcards.mem_type}/{wildcards.read_type}/pivot_{i}_align_dataset_{j}.sam")
    return input_files

####################################################
# Section 3: Rules needed for this experiment type
####################################################

#   Section 3.1: Generating long and short reads from each 
#   out-pivot genome for each dataset

rule generate_raw_positive_short_reads_exp7:
    input:
        "exp7_pivot_data/pivot_ref/pivot_{num}.fna"
    output:
        "exp7_ms_data/illumina/pivot_{num}/pivot_{num}.fna"
    shell:
        """
        art_illumina -ss HS25 -i exp7_pivot_data/pivot_ref/pivot_{wildcards.num}.fna -na -l 150 -f 8.0 \
        -o exp7_ms_data/illumina/pivot_{wildcards.num}/pivot_{wildcards.num}

        seqtk seq -a exp7_ms_data/illumina/pivot_{wildcards.num}/pivot_{wildcards.num}.fq > {output}
        """

rule generate_raw_positive_long_reads_exp7:
    input:
        "exp7_pivot_data/pivot_ref/pivot_{num}.fna"
    output:
        "exp7_ms_data/ont/pivot_{num}/pivot_{num}.fna"
    shell:
        """
        pbsim --depth 8.0 --prefix exp7_ms_data/ont/pivot_{wildcards.num}/pivot_{wildcards.num} \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 200 exp7_pivot_data/pivot_ref/pivot_{wildcards.num}.fna
        
        cat 'exp7_ms_data/ont/pivot_{wildcards.num}/pivot_{wildcards.num}'*.fastq > exp7_ms_data/ont/pivot_{wildcards.num}/pivot_{wildcards.num}.fastq
        seqtk seq -a exp7_ms_data/ont/pivot_{wildcards.num}/pivot_{wildcards.num}.fastq > {output}
        rm 'exp7_ms_data/ont/pivot_{wildcards.num}/pivot_{wildcards.num}_'*.fastq
        ls  'exp7_ms_data/ont/pivot_{wildcards.num}/pivot_{wildcards.num}'*.ref | xargs rm
        ls  'exp7_ms_data/ont/pivot_{wildcards.num}/pivot_{wildcards.num}'*.maf | xargs rm
        """

#   Section 3.2: Contains rules for building a concatenated reference of all
#   the genomes in each database. Then, building reverse complements, and 
#   concatenating those as well. SPUMONI just uses the forward sequences, but
#   we want the forward + rev. comp. for the r-index

rule build_forward_ref_dataset_exp7:
    input:
        get_dataset_non_pivot_genomes_exp7
    output:
        "exp7_combined_refs/dataset_{num}/combined_ref_forward.fna"
    shell:
        "cat {input} > {output}"

rule build_rev_comp_ref_dataset_exp7:
    input:
        "exp7_combined_refs/dataset_{num}/combined_ref_forward.fna"
    output:
        "exp7_combined_refs/dataset_{num}/combined_ref_rcomp.fna"
    shell:
        "seqtk seq -r {input} | sed 's/^>/>revcomp_/' > {output}"

rule concat_forward_and_rev_comp_dataset_exp7:
    input:
        "exp7_combined_refs/dataset_{num}/combined_ref_forward.fna",
        "exp7_combined_refs/dataset_{num}/combined_ref_rcomp.fna"
    output:
        "exp7_combined_refs/dataset_{num}/combined_ref_all.fna"
    shell:
        "cat {input} > {output}"

rule concat_only_forward_refs_for_all_datasets_exp7:
    input:
        get_combined_ref_dataset_exp7
    output:
        "exp7_combined_refs/combined_ref_forward_only.fna"
    shell:
        "cat {input} > {output}"

#   Section 3.3: Generate FASTA index that can be used to determine the 
#   total size of the input databases in order to figure out what the 
#   random length matching statistic will be with respect to database

rule build_fasta_index_for_combined_forward_exp7:
    input:
        "exp7_combined_refs/combined_ref_forward_only.fna"
    output:
        "exp7_combined_ref_fai/combined_ref_forward_only.fna.fai"
    shell:
        """
        cp {input} exp7_combined_ref_fai/combined_ref_forward_only.fna
        samtools faidx exp7_combined_ref_fai/combined_ref_forward_only.fna
        """

#   Section 3.4: Generate matching statistics for pivot reads with 
#   respect to the full database of all the non-pivot genomes.

rule build_spumoni_index_exp7:
    input:
        "exp7_combined_refs/combined_ref_forward_only.fna"
    output:
        "exp7_indexes/spumoni_index/spumoni_full_ref.fa",
        "exp7_indexes/spumoni_index/spumoni_full_ref.fa.thrbv.ms"
    shell:
        # Only use forward references for spumoni
        """
        cp {input} exp7_indexes/spumoni_index/combined_ref_forward_only.fna
        spumoni build -r exp7_indexes/spumoni_index/combined_ref_forward_only.fna -M -n
        """

rule run_spumoni_to_generate_ms_exp7:
    input:
        "exp7_indexes/spumoni_index/spumoni_full_ref.fa",
        "exp7_ms_data/{read_type}/pivot_{pivot}/pivot_{pivot}.fna"
    output:
        "exp7_ms_data/{read_type}/pivot_{pivot}/pivot_{pivot}.fna.lengths"
    shell:
        """
        spumoni run -r {input[0]} -p {input[1]} -M -n 
        """

#   Section 3.5: Extract half-MEMs and MEMs using the pivot's MS lengths
#   and the pivot itself.

rule extract_mems_or_halfmems_based_on_ms_exp7:
    input:
        "exp7_ms_data/{read_type}/pivot_{pivot}/pivot_{pivot}.fna.lengths",
        "exp7_ms_data/{read_type}/pivot_{pivot}/pivot_{pivot}.fna"
    output:
        "exp7_{mem_type}_data/{read_type}/pivot_{pivot}.fastq"
    shell:
        """
        python3 {repo_dir}/src/extract_mems.py \
        -l {input[0]} \
        -p {input[1]} \
        --{wildcards.mem_type} \
        -t {thresh} \
        -o {output}
        """

rule subset_mems_or_halfmems_exp7:
    input:
        "exp7_{mem_type}_data/{read_type}/pivot_{pivot}.fastq",
        "exp7_combined_ref_fai/combined_ref_forward_only.fna.fai"
    output:
        "exp7_{mem_type}_data/{read_type}/pivot_{pivot}_subset.fastq"
    shell:
        """
        python3 {repo_dir}/src/subset_reads.py \
        -i {input[0]} \
        -o {output} \
        --{wildcards.mem_type} \
        -n {num_of_non_kmers} \
        -l {input[1]} # fasta index to figure out noise factor
        """

#   Section 3.6: Build r-index for each dataset individually over
#   the reference containing both the forward and reverse complement.

rule build_rindex_over_individual_databases_exp7:
    input:
        "exp7_combined_refs/dataset_{num}/combined_ref_all.fna"
    output:
        "exp7_combined_refs/dataset_{num}/combined_ref_all.fna.ri"
    shell:
        """
        cd {r_dir} 
        ri-buildfasta {base_dir}/{input}
        cd {base_dir}
        """

#   Section 3.7: Create SAM file using r-index to locate all the 
#   half-MEMs or MEMs in each of the individual databases.

rule align_pivot_against_rindex_of_database_exp7:
    input:
        "exp7_combined_refs/dataset_{num}/combined_ref_all.fna",
        "exp7_combined_refs/dataset_{num}/combined_ref_all.fna.ri",
        "exp7_{mem_type}_data/{read_type}/pivot_{pivot}_subset.fastq"
    output:
        "exp7_sam_files/{mem_type}/{read_type}/pivot_{pivot}_align_dataset_{num}.sam"
    shell:
        """
        cd {r_dir} 
        ri-align -m 1 locate {base_dir}/{input[0]} {base_dir}/{input[2]} \
        > {base_dir}/{output}
        """

#   Section 3.8: Analyze the SAM files with respect to each database
#   to determine how we should split each 

rule analyze_sam_for_each_type_exp7:
    input:
        get_python_input_exp7,
        "exp7_combined_ref_fai/combined_ref_forward_only.fna.fai"
    output:
        "exp7_output_data/{mem_type}/{read_type}/accuracy_values.csv",
        "exp7_output_data/{mem_type}/{read_type}/confusion_matrix.csv"
    shell:
        """
        python3 {repo_dir}/src/analyze_sam.py \
        -n {num_datasets} \
        -s {base_dir}/exp7_sam_files/{wildcards.mem_type}/{wildcards.read_type}/ \
        -o {base_dir}/exp7_output_data/{wildcards.mem_type}/{wildcards.read_type}/ \
        --{wildcards.mem_type} \
        -l exp7_combined_ref_fai/combined_ref_forward_only.fna.fai
        """

#   Section 3.9: Gather all the output files needed for full
#   analysis of experiment 7

rule gather_all_output_files_exp7:
    input:
        "exp7_output_data/half_mems/illumina/accuracy_values.csv",
        "exp7_output_data/half_mems/ont/accuracy_values.csv",
        "exp7_output_data/mems/illumina/accuracy_values.csv",
        "exp7_output_data/mems/ont/accuracy_values.csv"
    output:
        expand("exp7_final_output/{mem_type}_{read_type}.csv", mem_type=["half_mems", "mems"], read_type=["illumina", "ont"])
    shell:
        """
        cp {input[0]} exp7_final_output/half_mems_illumina.csv
        cp {input[1]} exp7_final_output/half_mems_ont.csv
        cp {input[2]} exp7_final_output/mems_illumina.csv
        cp {input[3]} exp7_final_output/mems_ont.csv
        """
