####################################################
# Name: exp_type_7.smk
# Description: Contains functions and rules for
#              the type of experiment type 7: 
# Date: 7/5/22
####################################################

####################################################
# Section 1: Code that is always run, when this 
#            experiment is called. It genereates the
#            needed files.
####################################################

if exp_type == 7:
    if not os.path.isdir("non_pivot_type_7/"):
        for i in range(1, num_datasets+1):
            dir_prefix = f"data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"non_pivot_type_7/dataset_{i}"):
                os.makedirs(f"non_pivot_type_7/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"non_pivot_type_7/dataset_{i}")
            
            # Extract the pivot and move to new directory
            dir_prefix = f"non_pivot_type_7/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            print(list_of_files)

            pivot = random.choice(list_of_files)
            if not os.path.isdir(f"pivot_type_7/pivot_ref/"):
                os.makedirs(f"pivot_type_7/pivot_ref/")
            
            # Remove pivot from other set of genomes
            shutil.copy(pivot, f"pivot_type_7/pivot_ref/pivot_{i}.fna")
            os.remove(pivot)

####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_all_non_pivot_genomes_exp_7():
    """ Returns list of all non-pivot genomes """
    input_files = []
    for i in range(1, num_datasets + 1):
        for data_file in os.listdir(f"non_pivot_type_7/dataset_{i}/"):
            if data_file.endswith(".fna"):
                file_name = data_file.split(".fna")[0]
                input_files.append(f"non_pivot_type_7/dataset_{i}/{file_name}.fna")
    return input_files

####################################################
# Section 3: Rules needed for this experiment type
####################################################

# Section 3.1: Generating long and short reads from reference pivot genome

rule generate_raw_positive_short_reads_exp_7:
    input:
        "pivot_type_7/pivot_ref/pivot_{num}.fna"
    output:
        "pivot_type_7/pivot_raw_reads/illumina/pivot_{num}.fq"
    shell:
        """
        art_illumina -ss HS25 -i pivot_type_7/pivot_ref/pivot_{wildcards.num}.fna -na -l 150 -f 2.0 \
        -o pivot_type_7/pivot_raw_reads/illumina/pivot_{wildcards.num}
        """

rule generate_raw_positive_long_reads_exp_7:
    input:
        "pivot_type_7/pivot_ref/pivot_{num}.fna"
    output:
        "pivot_type_7/pivot_raw_reads/ont/pivot_{num}.fastq"
    shell:
        """
        pbsim --depth 30.0 --prefix pivot_type_7/pivot_raw_reads/ont/pivot_{wildcards.num} \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 200 pivot_type_7/pivot_ref/pivot_{wildcards.num}.fna
        
        cat 'pivot_type_7/pivot_raw_reads/ont/pivot_{wildcards.num}'*.fastq > {output}
        rm 'pivot_type_7/pivot_raw_reads/ont/pivot_{wildcards.num}_'*.fastq
        ls  'pivot_type_7/pivot_raw_reads/ont/pivot_{wildcards.num}'*.ref | xargs rm
        ls  'pivot_type_7/pivot_raw_reads/ont/pivot_{wildcards.num}'*.maf | xargs rm
        """

rule convert_long_reads_to_fasta_and_subset_exp_7:
    input:
        "pivot_type_7/pivot_raw_reads/ont/pivot_{num}.fastq"
    output:
        "ms_type_7/ont/pivot_{num}/pivot_{num}.fna"
    run:
        num_lines = num_reads_per_dataset * 4
        shell("head -n{num_lines} {input} > pivot_type_7/pivot_raw_reads/ont/pivot_{wildcards.num}_subset.fastq")
        shell("seqtk seq -a pivot_type_7/pivot_raw_reads/ont/pivot_{wildcards.num}_subset.fastq > {output}")
        shell("if [ $(grep -c '>' {output}) != {num_reads_per_dataset} ]; then echo 'number of reads assertion failed.'; exit 1; fi")
        shell("rm pivot_type_7/pivot_raw_reads/ont/pivot_{wildcards.num}_subset.fastq")

rule convert_pivot_to_fasta_and_subset_exp_7:
    input:
        "pivot_type_7/pivot_raw_reads/illumina/pivot_{num}.fq"
    output:
        "ms_type_7/illumina/pivot_{num}/pivot_{num}.fna"
    run:
        num_lines = num_reads_per_dataset * 4
        shell("head -n{num_lines} {input} > pivot_type_7/pivot_raw_reads/illumina/pivot_{wildcards.num}_subset.fq")
        shell("seqtk seq -a pivot_type_7/pivot_raw_reads/illumina/pivot_{wildcards.num}_subset.fq > {output}")
        shell("if [ $(grep -c '>' {output}) != {num_reads_per_dataset} ]; then echo 'number of reads assertion failed.'; exit 1; fi")
        shell("rm pivot_type_7/pivot_raw_reads/illumina/pivot_{wildcards.num}_subset.fq")


# Section 3.2: Concatonating ALL reference genomes and reverse complements

rule concat_ref_exp_7:
    input:
        get_all_non_pivot_genomes_exp_7()
    output:
        "combined_type_7/combined_ref_forward.fna"
    shell:
        "cat {input} > {output}"

rule build_rev_comp_exp_7:
    input:
        "combined_type_7/combined_ref_forward.fna"
    output:
        "combined_type_7/combined_ref_rcomp.fna"
    shell:
        "seqtk seq -r {input} | sed 's/^>/>revcomp_/' > {output}"

rule concat_full_ref_exp_7:
    input:
        "combined_type_7/combined_ref_rcomp.fna",
        "combined_type_7/combined_ref_forward.fna"
    output:
        "combined_type_7/combined_ref_all.fna"
    shell:
        "cat {input} > {output}"


# Section 3.3: Generate matching statistics for pivot reads

rule spumoni_build_exp_7:
    input:
        "combined_type_7/combined_ref_all.fna"
    output:
        "combined_type_7/spumoni_full_ref.fa"
    shell:
        "spumoni build -r combined_type_7/combined_ref_all.fna -M -n"

rule spumoni_run_exp_7:
    input:
        "combined_type_7/spumoni_full_ref.fa",
        "ms_type_7/{read_type}/pivot_{num}/pivot_{num}.fna"
    output:
        "ms_type_7/{read_type}/pivot_{num}/pivot_{num}.fna.lengths"
    shell:
        """
        spumoni run -r combined_type_7/spumoni_full_ref.fa \
        -p ms_type_7/{wildcards.read_type}/pivot_{wildcards.num}/pivot_{wildcards.num}.fna -M -n
        """
