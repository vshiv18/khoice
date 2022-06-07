# khoice - analysis of kmer length used for classification

This repository contains code to perform experiments focused on understanding the impact of kmer length (k) on sequence classification, such as metagenomic classification. Currently, the repository contains scripts to download genomes from the NCBI databases, and then looks at "discriminatory power" of each kmer by determining how many genomes the kmer occurs in within each group individually. Next, the workflow will look across groups to see how many groups each unique kmer occurs in to understand how strong it is for classification.

These experiments are inspired by a paper by [Nasko et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1554-6) where they showed the impact of RefSeq database growth on metagenomic classifiers that use kmers such as Kraken.

## Getting Started

To start out using `khoice`, you can download the repository and install the needed dependencies using the following commands:

```sh
git clone git@github.com:oma219/khoice.git
cd khoice
conda env create -f workflow/envs/khoice_exps.yaml
```
Now, you will have the `khoice_exps` environment that can be used to run the experiments. 

## Downloading Data

The first step before running any experiments is to download data from an NCBI database that you can use to create the plots. Currently, `khoice` contains a script called `src/download_genomes.py` that will allow you download genomes from different groups of species from the GenBank and RefSeq databases. **Importantly**, you need to activate the `khoice_exps` environment for this operation.

```sh
conda activate khoice_exps
python3 download_genomes.py -r -o <output_directory> -f <list_of_species> -n 5
```
The command above shows an example of how to run `src/download_genomes.py`, you first specify which database to download genomes from using either Refseq (`-r`) or GenBank (`-g`) option. Next, you need to specify an output directory where these genomes will be saved. Note, within that directory, a `data/` folder will be created to store all the genomes requested. 

Then, species can be requested either through file input where you specify a path to file with `-f` option, and it contains species name on each line. Or, you can use the console to type in species to download. Lastly, the `-n` option allows you to specify how many genomes you want to download from each species requested. The default is to download all the genomes available for that species in the requested database.

## Running Experiments

Now that we have data downloaded in our specified output directory, we can start to run the workflow to gather the data we need for our experiments. Currently, the workflow relies on the [KMC3: kmer-counter](https://github.com/refresh-bio/KMC) to perform the kmer set operations and manipulations.

```sh
conda activate khoice_exps
snakemake -np --config DATA_ROOT=<output_directory>
```
First, you can perform a dry-run of the workflow to make sure everything is as expected. It is important to use the same output directory path you used when downloading the data.

```sh
snakemake --cores N --config DATA_ROOT=<output_directory>
```

Next, you can actually run the workflow as shown above. After completing the workflow, you can obtain the `csv` files with the results in following folders: `step_5/` and `step_9` in your output directory.

