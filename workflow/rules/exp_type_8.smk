####################################################
# Name: exp_type_8.smk
# Description: Contains functions and rules for
#              the type of experiment type 7:
# 
#
# Date: 7/6/22
####################################################

####################################################
# Section 1: Code that is always run, when this 
#            experiment is called. It genereates the
#            needed files.
####################################################

# IMPORTANT: assumes input fasta unzipped (TODO: add unzip option)

if exp_type == 8:
    if not os.path.isdir("non_pivot_type_8/"):
        for i in range(1, num_datasets+1):
            dir_prefix = f"data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"non_pivot_type_8/dataset_{i}"):
                os.makedirs(f"non_pivot_type_8/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"non_pivot_type_8/dataset_{i}")
            
            # Extract the pivot and move to new directory
            dir_prefix = f"non_pivot_type_8/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            print(list_of_files)

            pivot = random.choice(list_of_files)
            if not os.path.isdir(f"pivot_type_8/pivot_ref/"):
                os.makedirs(f"pivot_type_8/pivot_ref/")
            
            # Remove pivot from other set of genomes
            shutil.copy(pivot, f"pivot_type_8/pivot_ref/pivot_{i}.fna")
            os.remove(pivot)
        
####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_all_non_pivot_genomes_exp_8():
    """ Returns list of all non-pivot genomes """
    input_files = []
    if exp_type == 8:
        for i in range(1, num_datasets + 1):
            for data_file in os.listdir(f"non_pivot_type_8/dataset_{i}/"):
                if data_file.endswith(".fna"):
                    file_name = data_file.split(".fna")[0]
                    input_files.append(f"non_pivot_type_8/dataset_{i}/{file_name}.fna")
    return input_files

def get_dataset_non_pivot_genomes_exp_8(wildcards):
    """ Returns list of non-pivot genomes for a given dataset """
    input_files = []
    for data_file in os.listdir(f"non_pivot_type_8/dataset_{wildcards.num}/"):
        if data_file.endswith(".fna"):
            file_name = data_file.split(".fna")[0]
            input_files.append(f"non_pivot_type_8/dataset_{wildcards}/{file_name}.fna")
    return input_files

def get_python_input_exp_8(wildcards):
    """ Returns list of text files with reference names for each dataset AND all pivot SAM files for a read type """
    input_files = []
    for i in range(1,num_datasets+1):
        input_files.append(f"ref_lists_type_8/dataset_{i}_references.txt")
        input_files.append(f"sam_type_8/{wildcards.mem_type}/{wildcards.read_type}/pivot_{i}.sam")
    print(input_files)
    return input_files

####################################################
# Section 3: Rules needed for this experiment type
####################################################

# Section 3.1: Generating long and short reads from reference pivot genome

rule generate_raw_positive_short_reads_exp_8:
    input:
        "pivot_type_8/pivot_ref/pivot_{num}.fna"
    output:
        "pivot_type_8/pivot_raw_reads/illumina/pivot_{num}.fq"
    shell:
        """
        art_illumina -ss HS25 -i pivot_type_8/pivot_ref/pivot_{wildcards.num}.fna -na -l 150 -f 2.0 \
        -o pivot_type_8/pivot_raw_reads/illumina/pivot_{wildcards.num}
        """

rule generate_raw_positive_long_reads_exp_8:
    input:
        "pivot_type_8/pivot_ref/pivot_{num}.fna"
    output:
        "pivot_type_8/pivot_raw_reads/ont/pivot_{num}.fastq"
    shell:
        """
        pbsim --depth 30.0 --prefix pivot_type_8/pivot_raw_reads/ont/pivot_{wildcards.num} \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 200 pivot_type_8/pivot_ref/pivot_{wildcards.num}.fna
        
        cat 'pivot_type_8/pivot_raw_reads/ont/pivot_{wildcards.num}'*.fastq > {output}
        rm 'pivot_type_8/pivot_raw_reads/ont/pivot_{wildcards.num}_'*.fastq
        ls  'pivot_type_8/pivot_raw_reads/ont/pivot_{wildcards.num}'*.ref | xargs rm
        ls  'pivot_type_8/pivot_raw_reads/ont/pivot_{wildcards.num}'*.maf | xargs rm
        """

rule convert_long_reads_to_fasta_and_subset_exp_8:
    input:
        "pivot_type_8/pivot_raw_reads/ont/pivot_{num}.fastq"
    output:
        "ms_type_8/ont/pivot_{num}/pivot_{num}.fna"
    run:
        num_lines = num_reads_per_dataset * 4
        shell("head -n{num_lines} {input} > pivot_type_8/pivot_raw_reads/ont/pivot_{wildcards.num}_subset.fastq")
        shell("seqtk seq -a pivot_type_8/pivot_raw_reads/ont/pivot_{wildcards.num}_subset.fastq > {output}")
        shell("if [ $(grep -c '>' {output}) != {num_reads_per_dataset} ]; then echo 'number of reads assertion failed.'; exit 1; fi")
        shell("rm pivot_type_8/pivot_raw_reads/ont/pivot_{wildcards.num}_subset.fastq")

rule convert_pivot_to_fasta_and_subset_exp_8:
    input:
        "pivot_type_8/pivot_raw_reads/illumina/pivot_{num}.fq"
    output:
        "ms_type_8/illumina/pivot_{num}/pivot_{num}.fna"
    run:
        num_lines = num_reads_per_dataset * 4
        shell("head -n{num_lines} {input} > pivot_type_8/pivot_raw_reads/illumina/pivot_{wildcards.num}_subset.fq")
        shell("seqtk seq -a pivot_type_8/pivot_raw_reads/illumina/pivot_{wildcards.num}_subset.fq > {output}")
        shell("if [ $(grep -c '>' {output}) != {num_reads_per_dataset} ]; then echo 'number of reads assertion failed.'; exit 1; fi")
        shell("rm pivot_type_8/pivot_raw_reads/illumina/pivot_{wildcards.num}_subset.fq")

# Section 3.2: Build reference combined reference database with reverse complement

rule concat_ref_exp_8:
    input:
        get_all_non_pivot_genomes_exp_8()
    output:
        "combined_type_8/combined_ref_forward.fna"
    shell:
        "cat {input} > {output}"

rule build_rev_comp_exp_8:
    input:
        "combined_type_8/combined_ref_forward.fna"
    output:
        "combined_type_8/combined_ref_rcomp.fna"
    shell:
        "seqtk seq -r {input} | sed 's/^>/>revcomp_/' > {output}"

rule concat_full_ref_exp_8:
    input:
        "combined_type_8/combined_ref_rcomp.fna",
        "combined_type_8/combined_ref_forward.fna"
    output:
        "combined_type_8/combined_ref_all.fna"
    shell:
        "cat {input} > {output}"

# Section 3.3: Generate matching statistics for pivot

rule spumoni_build_exp_8:
    input:
        "combined_type_8/combined_ref_all.fna"
    output:
        "combined_type_8/spumoni_full_ref.fa"
    shell:
        "spumoni build -r combined_type_8/combined_ref_all.fna -M -n"

rule spumoni_run_exp_8:
    input:
        "combined_type_8/spumoni_full_ref.fa",
        "ms_type_8/{read_type}/pivot_{num}/pivot_{num}.fna"
    output:
        "ms_type_8/{read_type}/pivot_{num}/pivot_{num}.fna.lengths"
    shell:
        """
        spumoni run -r {input[0]} \
        -p {input[1]} -M -n
        """

# Section 3.4: Extract Half-Mems and Mems from pivot MS with no threshold

rule extract_mems_exp_8:
    input:
        "ms_type_8/{read_type}/pivot_{num}/pivot_{num}.fna.lengths",
        "ms_type_8/{read_type}/pivot_{num}/pivot_{num}.fna"
    output:
        "{mem_type}_type_8/{read_type}/pivot_{num}.fastq"
    shell:
        """
        python3 {repo_dir}/src/extract_mems.py \
        -l {base_dir}/{input[0]} \
        -p {base_dir}/{input[1]} \
        --{wildcards.mem_type} \
        -t 2 \
        -o {base_dir}/{output}
        """

# Section 3.4: Build r-index

rule ref_ri_build_fasta_exp_8:
    input:
        "combined_type_8/combined_ref_all.fna"
    output:
        "combined_type_8/combined_ref_all.fna.ri"
    shell:
        """
        cd {r_dir} 
        ri-buildfasta {base_dir}/{input}
        cd {base_dir}
        """

# Section 3.4: Create SAM file with located Half-Mems across threshold values

rule ref_ri_align_mem_exp_8:
    input:
        "combined_type_8/combined_ref_all.fna",
        "{mem_type}_type_8/{read_type}/pivot_{num}.fastq",
        "combined_type_8/combined_ref_all.fna.ri"
    output:
        "sam_type_8/{mem_type}/{read_type}/pivot_{num}.sam"
    shell:
        """
        cd {r_dir} 
        ri-align locate {base_dir}/combined_type_8/combined_ref_all.fna {base_dir}/{wildcards.mem_type}_type_8/{wildcards.read_type}/pivot_{wildcards.num}.fastq \
        > {base_dir}/{output}
        """

# Section 3.5: Generate List of Reference Genomes

rule gen_ref_list_exp_8:
    input:
        get_dataset_non_pivot_genomes_exp_8
    output:
        "ref_lists_type_8/dataset_{num}_references.txt"
    shell:
        """
        for file in {input}; do
            grep '^>' $file | awk '{{print substr($1,2)}}' >> "ref_lists_type_8/dataset_{wildcards.num}_references.txt"
        done
        """

rule analyze_sam_exp_8:
    input:
        get_python_input_exp_8
    output:
        "output_type_8/{mem_type}/t_{t}/{read_type}/confusion_matrix.csv"
    shell:
        """
        python3 {repo_dir}/src/analyze_sam.py \
        -n {num_datasets} \
        -s {base_dir}/sam/{wildcards.mem_type}/{wildcards.read_type}/ \
        -r {base_dir}/ref_lists_type_8/ \
        -o {base_dir}/output_type_8/{wildcards.mem_type}/t_{wildcards.t}/{wildcards.read_type}/ \
        -t {wildcards.t} \
        --{wildcards.mem_type}
        """