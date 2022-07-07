####################################################
# Name: exp_type_5.smk
# Description: Contains functions and rules for
#              the type of experiment type 5: 
#
#              Selects an out-pivot genome from each species, 
#              finds its half-mems/mems with respect to a combined 
#              database, and builds a confusion matrix.
#
# Date: 6/22/22
####################################################

####################################################
# Section 1: Code that is always run, when this 
#            experiment is called. It genereates the
#            needed files.
####################################################

if exp_type == 5:
    if not os.path.isdir("non_pivot_type_5/"):
        for i in range(1, num_datasets+1):
            dir_prefix = f"data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"non_pivot_type_5/dataset_{i}"):
                os.makedirs(f"non_pivot_type_5/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"non_pivot_type_5/dataset_{i}")
            
            # Extract the pivot and move to new directory
            dir_prefix = f"non_pivot_type_5/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            print(list_of_files)

            pivot = random.choice(list_of_files)
            if not os.path.isdir(f"pivot_ms_type_5/pivot_{i}"):
                os.makedirs(f"pivot_ms_type_5/pivot_{i}")
            
            # Remove pivot from other set of genomes
            shutil.copy(pivot, f"pivot_ms_type_5/pivot_{i}/pivot_{i}.fna")
            os.remove(pivot)
            
####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_all_non_pivot_genomes_exp_5():
    """ Returns list of all non-pivot genomes """
    input_files = []
    if exp_type == 5:
        for i in range(1, num_datasets + 1):
            for data_file in os.listdir(f"non_pivot_type_5/dataset_{i}/"):
                if data_file.endswith(".fna"):
                    file_name = data_file.split(".fna")[0]
                    input_files.append(f"non_pivot_type_5/dataset_{i}/{file_name}.fna")
    return input_files

def get_dataset_non_pivot_genomes_exp_5(wildcards):
    """ Returns list of non-pivot genomes for a given dataset """
    input_files = []
    if(exp_type == 5):
        for data_file in os.listdir(f"non_pivot_type_5/dataset_{wildcards.num}/"):
            if data_file.endswith(".fna"):
                file_name = data_file.split(".fna")[0]
                input_files.append(f"non_pivot_type_5/dataset_{wildcards}/{file_name}.fna")
    return input_files

def get_python_input_exp_5(wildcards):
    """ Returns list of text files with reference names for each dataset AND all pivot SAM files for mem type """
    input_files = []
    for i in range(1,num_datasets+1):
        input_files.append(f"ref_lists_type_5/dataset_{i}_references.txt")
        input_files.append(f"sam_type_5/{wildcards.mem_type}/pivot_{i}.sam")
    return input_files


####################################################
# Section 3: Rules needed for this experiment type
####################################################

# Section 3.1: 

rule concat_ref_exp_5:
    input:
        get_all_non_pivot_genomes_exp_5()
    output:
        "combined_type_5/combined_ref_forward.fna"
    shell:
        "cat {input} > {output}"

rule build_rev_comp_exp_5:
    input:
        "combined_type_5/combined_ref_forward.fna"
    output:
        "combined_type_5/combined_ref_rcomp.fna"
    shell:
        "seqtk seq -r {input} | sed 's/^>/>revcomp_/' > {output}"

rule concat_full_ref_exp_5:
    input:
        "combined_type_5/combined_ref_rcomp.fna",
        "combined_type_5/combined_ref_forward.fna"
    output:
        "combined_type_5/combined_ref_all.fna"
    shell:
        "cat {input} > {output}"


# Section 3.2: Generate matching statistics for pivot

rule spumoni_build_exp_5:
    input:
        "combined_type_5/combined_ref_all.fna"
    output:
        "combined_type_5/spumoni_full_ref.fa"
    shell:
        "spumoni build -r combined_type_5/combined_ref_all.fna -M -n"

rule spumoni_run_exp_5:
    input:
        "combined_type_5/spumoni_full_ref.fa",
        "pivot_ms_type_5/pivot_{num}/pivot_{num}.fna"
    output:
        "pivot_ms_type_5/pivot_{num}/pivot_{num}.fna.lengths"
    shell:
        """
        spumoni run -r combined_type_5/spumoni_full_ref.fa \
        -p pivot_ms_type_5/pivot_{wildcards.num}/pivot_{wildcards.num}.fna -M -n
        """

# Section 3.3: Extract Half-Mems and Mems from pivot MS

rule extract_mems_exp5:
    input:
        "pivot_ms_type_5/pivot_{num}/pivot_{num}.fna.lengths",
        "pivot_ms_type_5/pivot_{num}/pivot_{num}.fna"
    output:
        "{mem_type}_type_5/pivot_{num}.fastq"
    shell:
        """
        python3 {repo_dir}/src/extract_mems.py \
        -l {base_dir}/{input[0]} \
        -p {base_dir}/{input[1]} \
        --{wildcards.mem_type} \
        -t {thresh} \
        -o {base_dir}/{output}
        """

# Section 3.4: Build r-index

rule ref_ri_build_fasta_exp_5:
    input:
        "combined_type_5/combined_ref_all.fna"
    output:
        "combined_type_5/combined_ref_all.fna.ri"
    shell:
        """
        cd {r_dir} 
        ri-buildfasta {base_dir}/combined_type_5/combined_ref_all.fna 
        cd {base_dir}
        """

# Section 3.4: Create SAM file with located Half-Mems/Mems

rule ref_ri_align_mem_exp_5:
    input:
        "combined_type_5/combined_ref_all.fna",
        "{mem_type}_type_5/pivot_{num}.fastq",
        "combined_type_5/combined_ref_all.fna.ri"
        
    output:
        "sam_type_5/{mem_type}/pivot_{num}.sam"
    shell:
        """
        cd {r_dir} 
        ri-align locate {base_dir}/{input[0]} {base_dir}/{input[1]} \
        > {base_dir}/{output}
        """

# Section 3.5: Generate List of Reference Genomes

rule gen_ref_list_exp_5:
    input:
        get_dataset_non_pivot_genomes_exp_5
    output:
        "ref_lists_type_5/dataset_{num}_references.txt"
    shell:
        """
        for file in {input}; do
            grep '^>' $file | awk '{{print substr($1,2)}}' >> "ref_lists_type_5/dataset_{wildcards.num}_references.txt"
        done
        """

# Section 3.6 Outputs confusion matrix for either half-mem/mems

rule analyze_sam_exp_5:
    input:
        get_python_input_exp_5
    output:
        "output_type_5/{mem_type}/confusion_matrix.csv",
    shell:
        """
        python3 {repo_dir}/src/analyze_sam.py \
        -n {num_datasets} \
        -s {base_dir}/sam_type_5/{wildcards.mem_type}/ \
        -r {base_dir}/ref_lists_type_5/ \
        -o {base_dir}/output_type_5/{wildcards.mem_type}/ \
        --{wildcards.mem_type}
        """
