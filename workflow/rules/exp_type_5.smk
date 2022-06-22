####################################################
# Name: exp_type_5.smk
# Description: Contains functions and rules for
#              the type of experiment type 5: 
# Date: 6/22/22
####################################################

####################################################
# Section 1: Code that is always run, when this 
#            experiment is called. It genereates the
#            needed files.
####################################################

if exp_type == 5:
    if not os.path.isdir("input_type_5/"):
        for i in range(1, num_datasets+1):
            dir_prefix = f"data/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            if not os.path.isdir(f"input_type_5/rest_of_set/dataset_{i}"):
                os.makedirs(f"input_type_5/rest_of_set/dataset_{i}")

            # Copy over all the files at first
            for data_file in list_of_files:
                shutil.copy(data_file, f"input_type_5/rest_of_set/dataset_{i}")
            
            # Extract the pivot and move to new directory
            dir_prefix = f"input_type_5/rest_of_set/dataset_{i}/"
            list_of_files = [dir_prefix + x for x in os.listdir(dir_prefix)]

            pivot = random.choice(list_of_files)
            if not os.path.isdir(f"input_type_5/pivot/dataset_{i}"):
                os.makedirs(f"input_type_5/pivot/dataset_{i}")
            
            # Remove pivot from other set of genomes
            shutil.copy(pivot, f"input_type_5/pivot/dataset_{i}/pivot_{i}.fna.gz")
            os.remove(pivot)
            
####################################################
# Section 2: Helper functions needed for these
#            experiment rules.
####################################################

def get_all_non_pivot_genomes(wildcards):
    """ Returns a list of database in a certain dataset """
    input_files = []
    for i in range(1, )
    for data_file in os.listdir(f"input_type_5/rest_of_set/dataset_{wildcards.num}/"):
        if data_file.endswith(".fna.gz"):
            file_name = data_file.split(".fna")[0]
            input_files.append(f"input_type_5/rest_of_set/dataset_{wildcards.num}/file_name.fna")
    return input_files

####################################################
# Section 3: Rules needed for this experiment type
####################################################

# Section 3.1: 

rule concat_ref_exp_5:
    input:
        get_all_non_pivot_genomes
    output:
        "combined/combined_ref.fna"
    shell:
        "cat input_type_5/rest_of_set/dataset_{wildcards.num}/*_genomic.fna > combined/dataset_{wildcards.num}/combined_ref.fna"

# Section 3.2: Generate matching statistics for pivot

rule spumoni_build_exp_5:
    input:
        "combined/dataset_{num}/combined_ref.fna"
    output:
        "combined/dataset_{num}/spumoni_full_ref.fa"
    shell:
        "spumoni build -r combined/dataset_{wildcards.num}/combined_ref.fna -M -n"

rule spumoni_run_exp_5:
    input:
        "combined/dataset_{num}/spumoni_full_ref.fa",
        "input_type_5/pivot/dataset_{pivot}/pivot_{pivot}.fna"
    output:
        "input_type_5/pivot/dataset_{num}/pivot_{pivot}.fna.lengths"
    shell:
        """
        spumoni run -r combined/dataset_{wildcards.num}/spumoni_full_ref.fa \
        -p input_type_5/pivot/dataset_{wildcards.pivot}/pivot_{wildcards.pivot}.fna -M -n
        """