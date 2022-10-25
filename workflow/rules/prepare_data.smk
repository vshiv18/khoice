####################################################
# Name: prepare_data.smk
# Description: Takes in a set of genomes of different
#              classes, generates the out-pivot 
#              dataset, and simulates both short
#              and long reads.
#
# Date: 10/25/22
####################################################

"""
    # Parse out a pivot from the downloaded dataset
    if not os.path.isdir("exp6_input/"):
        for i in range(1, num_datasets+1):
            dir_prefix = f"data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"exp6_input/rest_of_set/dataset_{i}"):
                os.makedirs(f"exp6_input/rest_of_set/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"exp6_input/rest_of_set/dataset_{i}")
            
            # Extract the pivot and move to new directory
            dir_prefix = f"exp6_input/rest_of_set/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            pivot = random.choice(list_of_files)

            # Print pivot name to list in trial folder
            with open(trial_info_dir + "/exp6_pivot_names.txt", 'a') as fd:
                fd.write(pivot+"\n")

            if not os.path.isdir(f"exp6_input/pivot/"):
                os.makedirs(f"exp6_input/pivot/")
            
            # Remove pivot from other set of genomes
            shutil.copy(pivot, f"exp6_input/pivot/pivot_{i}.fna.gz")
            os.remove(pivot)
"""

####################################################
# Section 1: Functions used by rules below
####################################################

def get_all_genomes_in_dataset_exp0(wildcards):
    """ Returns a list of database in a certain dataset """
    input_files = []
    for data_file in os.listdir(f"{database_root}/dataset_{wildcards.num}/"):
        if data_file.endswith(".fna.gz"):
            input_files.append(f"{database_root}/dataset_{wildcards.num}/{data_file}")
    return input_files

####################################################
# Section 2: Rules needed for this workflow
####################################################

#   Section 2.1: Take pivot genomes from dataset, and move
#   them to a separate folder

rule data_prep_choose_pivot_from_each_dataset:
    input:
        get_all_genomes_in_dataset_exp0
    output:
        "exp0_pivot_genomes/dataset_{num}/pivot_{num}.fna.gz",
        "exp0_pivot_genomes/dataset_{num}/pivot_name.txt",
        "exp0_nonpivot_genomes/dataset_{num}/nonpivot_names.txt"
    shell:
        """
        # Copy all genomes to non-pivot folder
        cp {input} exp0_nonpivot_genomes/dataset_{wildcards.num}/

        # Choose a pivot genome from the non-pivot folder
        pivot_genome_path=$(ls exp0_nonpivot_genomes/dataset_{wildcards.num}/*.fna.gz | shuf | head -n1)
        pivot_genome_name=$(echo $pivot_genome_path  | awk -F/ '{{print $NF}}')

        # Copy pivot to pivot folder, save the name, and remove it
        cp $pivot_genome_path {output[0]}
        echo $pivot_genome_name > {output[1]}
        rm $pivot_genome_path

        # Save all the names of the non-pivot genomes to file
        for non_pivot in $(ls exp0_nonpivot_genomes/dataset_{wildcards.num}/*.fna.gz); do
            echo $non_pivot | awk -F/ '{{print $NF}}' >> {output[2]}
        done
        """

#   Section 2.2: Simulate short and long reads from a pivot genome
#   from each dataset. And subset reads to reach a certain threshold
#   based on config parameters

rule generate_raw_positive_short_reads_exp0:
    input:
        "exp0_pivot_genomes/dataset_{num}/pivot_{num}.fna.gz"
    output:
        "exp0_pivot_reads/dataset_{num}/illumina/pivot_{num}.fa"
    shell:
        """
        gzip -d {input} -k -f
        art_illumina -ss HS25 -i exp0_pivot_genomes/dataset_{wildcards.num}/pivot_{wildcards.num}.fna -na -l 150 -f 50.0 \
        -o exp0_pivot_reads/dataset_{wildcards.num}/illumina/pivot_{wildcards.num}
        rm exp0_pivot_genomes/dataset_{wildcards.num}/pivot_{wildcards.num}.fna

        seqtk seq -a exp0_pivot_reads/dataset_{wildcards.num}/illumina/pivot_{wildcards.num}.fq > {output}
        rm exp0_pivot_reads/dataset_{wildcards.num}/illumina/pivot_{wildcards.num}.fq 
        """

rule generate_raw_positive_long_reads_exp0:
    input:
        "exp0_pivot_genomes/dataset_{num}/pivot_{num}.fna.gz"
    output:
        "exp0_pivot_reads/dataset_{num}/ont/pivot_{num}.fa"
    shell:
        """
        gzip -d {input} -k -f
        pbsim --depth 50.0 --prefix exp0_pivot_reads/dataset_{wildcards.num}/ont/pivot_{wildcards.num} \
        --hmm_model {pbsim_model} --accuracy-mean 0.95 --length-min 200 exp0_pivot_genomes/dataset_{wildcards.num}/pivot_{wildcards.num}.fna
        
        cat 'exp0_pivot_reads/dataset_{wildcards.num}/ont/pivot_{wildcards.num}'*.fastq > exp0_pivot_reads/dataset_{wildcards.num}/ont/pivot_{wildcards.num}.fastq
        seqtk seq -a exp0_pivot_reads/dataset_{wildcards.num}/ont/pivot_{wildcards.num}.fastq > exp0_pivot_reads/dataset_{wildcards.num}/ont/pivot_{wildcards.num}.fa
        rm 'exp0_pivot_reads/dataset_{wildcards.num}/ont/pivot_{wildcards.num}_'*.fastq
        ls  'exp0_pivot_reads/dataset_{wildcards.num}/ont/pivot_{wildcards.num}'*.ref | xargs rm
        ls  'exp0_pivot_reads/dataset_{wildcards.num}/ont/pivot_{wildcards.num}'*.maf | xargs rm
        rm exp0_pivot_reads/dataset_{wildcards.num}/ont/pivot_{wildcards.num}.fastq
        """

rule subset_raw_reads_exp0:
    input:
        "exp0_pivot_reads/dataset_{num}/{read_type}/pivot_{num}.fa"
    output:
        "exp0_pivot_reads/dataset_{num}/{read_type}/pivot_{num}_subset.fa"
    shell:
        """
        # Just use k=31 since it a middle-of-the-range value ...
        python3 {repo_dir}/src/subset_reads.py -i {input} -o {output} -n {num_kmers} -k 31 --kmers
        """

#   Section 2.2: Generate a summary file for the dataset that includes all the relevant
#   values of interest like pivot, non-pivot files and the number of reads of each.

rule generate_database_summary_file_exp0:
    input:
        expand("exp0_pivot_genomes/dataset_{num}/pivot_name.txt", num=range(1, num_datasets+1)),
        expand("exp0_nonpivot_genomes/dataset_{num}/nonpivot_names.txt", num=range(1, num_datasets+1)),
        expand("exp0_pivot_reads/dataset_{num}/{read_type}/pivot_{num}_subset.fa", num=range(1, num_datasets+1), read_type=['illumina', 'ont'])
    output:
        "exp0_summary/database_summary.txt"
    shell:
        """
        # Create the header line
        printf "Dataset #: | " > {output[0]}
        for i in $(seq 1 {num_datasets}); do
            printf "$i | " >> {output[0]}
        done
        printf "\n" >> {output[0]}

        # Grab the pivot genome names ..
        printf "Pivot Genome: | " >> {output[0]}
        for i in $(seq 1 {num_datasets}); do
            pivot_name=$(head -n1 "exp0_pivot_genomes/dataset_${{i}}/pivot_name.txt")
            printf "$pivot_name | " >> {output[0]}
        done
        printf "\n" >> {output[0]}

        # Print number of Illumina reads ...
        printf "# of Illumina Reads: | " >> {output[0]}
        for i in $(seq 1 {num_datasets}); do
            num_reads=$(grep -c "^>" "exp0_pivot_reads/dataset_${{i}}/illumina/pivot_${{i}}_subset.fa")
            printf "$num_reads | " >> {output[0]}
        done
        printf "\n" >> {output[0]}

        # Print number of ONT reads ...
        printf "# of ONT Reads: | " >> {output[0]}
        for i in $(seq 1 {num_datasets}); do
            num_reads=$(grep -c "^>" "exp0_pivot_reads/dataset_${{i}}/ont/pivot_${{i}}_subset.fa")
            printf "$num_reads | " >> {output[0]}
        done
        printf "\n" >> {output[0]}

        # Get the non-pivot file with the largest number of lines ...
        max_lines=0
        for i in $(seq 1 {num_datasets}); do
            num_lines=$(wc -l "exp0_nonpivot_genomes/dataset_${{i}}/nonpivot_names.txt" | awk '{{print $1}}')
            if (( $num_lines > $max_lines )); then
                max_lines=$((num_lines))
            fi 
        done

        # Print the non-pivot genomes to file ...
        printf "Non-Pivot genomes: \n" >> {output[0]}
        for line in $(seq 1 $max_lines); do
            printf " | " >> {output[0]}
            for i in $(seq 1 {num_datasets}); do
                genome_name=$(cat "exp0_nonpivot_genomes/dataset_${{i}}/nonpivot_names.txt" | awk -v var=$line 'NR==var')
                printf "$genome_name | " >> {output[0]}
            done
            printf "\n" >> {output[0]}
        done
        """
