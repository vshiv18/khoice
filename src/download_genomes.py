#!/usr/bin/env python3

# Name: download_genomes.py 
# Description: This is a python script that will allow you to request certain
#              groups of genomes from NCBI databases like RefSeq and GenBank.
#
# Date: January 8th, 2022

import argparse
import os
import subprocess
import time

def run_command(command, loading_bar, sleep_time=1):
    """ 
    Takes in a command, and runs it. Presents a loading bar in case a command takes 
    some time to run like queries to databases or downloading genomes.
    """
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    
    if loading_bar:
        print("Loading ...", end='', flush="True")
        while process.poll() is None:
            time.sleep(sleep_time)
            print(".", end='', flush="True")
            (output, error) = process.communicate()
        (output, error) = process.communicate()
        print()
    else:
        (output, error) = process.communicate()

    output = output.decode("utf-8").strip() if output is not None else output
    error = error.decode("utf-8").strip() if error is not None else error
    return (output, error)

##################################################################
# NOTE: There is another way to find the FTP links for the 
# genomes that is slightly easier, it involves downloading the 
# assembly file, and then parsing it for what you are interested
# in. The approach is described in the following link:
# https://www.biostars.org/p/255212/
##################################################################

def find_refseq_ftp_paths(requested_species):
    """ Finds all the RefSeq FTP links for requested species """
    command = ("esearch -db assembly -query '(\"{}\"[Organism] OR {}[All Fields]) AND ((latest[filter] OR \"latest refseq\"[filter]) AND \"complete genome\"[filter] AND all[filter] NOT anomalous[filter])'"
               " | esummary | "
               " xtract -pattern DocumentSummary -block FtpPath -match '@type:RefSeq' -element FtpPath " 
               " | sed 's/$/\//' ")

    # Executes the command to figure out how many assemblies there are ...
    current_command = command.format(requested_species, requested_species) + " | wc -l"
    output, error = run_command(current_command, True)

    # Interprets output, and asks user if they want to download the genomes ...
    if output.isdigit() and int(output) > 0:
        response = input(f"\nWe found {output} complete genomes for {requested_species}. Does that\nsound right, would you like to download them (y/n)? ").upper().strip()
        while response != "Y" and response != "N":
            response = input("\nError in response, do you want to download or not (y/n)? ").upper().strip()
        if response == "N":
            return
    elif output.isdigit() and int(output) == 0:
        print("\nError: No genomes were found for this input.")
        exit(1)
    else:
        print(f"\nError occurred during database query: {output}")
        exit(1)

    current_command = command.format(requested_species, requested_species)
    output, error = run_command(current_command, False)
    output_lines = output.split("\n")

    genome_paths = []
    for ftp_path in output_lines:
        name = ftp_path.split("/")[-2]
        genome_paths.append(ftp_path + name + "_genomic.fna.gz")
    return genome_paths

def find_genbank_ftp_paths(requested_species):
    """ Finds all the GenBank FTP links for requested species """
    raise NotImplementedError("still working on this ...")

def download_genomes(ftp_paths, output_dir, num):
    """ Takes a list of paths, and downloads them to output directory"""
    with open(output_dir + "ftp_paths_for_curr_dataset.txt", "w") as fd:
        for paths in ftp_paths[:10]:
            fd.write(paths + "\n")
    
    filelist = output_dir + "ftp_paths_for_curr_dataset.txt"
    current_dataset_dir = output_dir + "dataset_" + str(num)

    # Executes a command to download genomes, remove filelist and
    # then decompress all the genomes
    
    print(f"\nStarting to download files to dataset_{num} folder ...")
    download_command = f"wget --input-file={filelist} --directory-prefix={current_dataset_dir} --no-parent --no-host-directories --cut-dirs=5"
    run_command(download_command, True, sleep_time=10)

    remove_command = f"rm {filelist} "
    run_command(remove_command, False)

    """
    Commented this out for now, KMC can read compressed files ....
    """
    #print(f"\nDecompressing those downloaded files now ...")
    #decompress_command = f"for i in $(find {current_dataset_dir} -name '*.fna.gz'); do gzip -d $i; done"
    #output, error = run_command(decompress_command, True)

def print_status(args):
    """ It will just state what parameters the user as specified prior to downloading any genomes """
    database = "RefSeq" if args.refseq else "GenBank"

    print("\nProvided Parameters:")
    print(f"\tGenome Database: {database}")
    print(f"\tOutput Directory: {args.output_dir}\n")

def main(args):
    """
    Main method:
        It will prompt the user continuously for a species name that they would like
        to download to the pre-specified directory, and it will continue to ask until
        user is done.
    """
    dataset_num = 1
    dataset_summaries = []

    downloaded_genomes = False
    print_status(args)
    while True:
        requested_species = input("Enter the name of species you want to download (type \"done\" to stop): ").strip()
        if requested_species.upper() == "DONE":
            break

        # Determine the paths to genomes for specified species
        if args.refseq:
            paths = find_refseq_ftp_paths(requested_species)
        elif args.genbank:
            paths = find_genbank_ftp_paths(requested_species)

        # Execute the download of those genomes ...
        download_genomes(paths, args.output_dir, dataset_num)
        downloaded_genomes = True

        dataset_summaries.append([dataset_num, requested_species, len(paths)])
        dataset_num += 1
        print()
    
    # Write out a README explaining what is in each folder ...
    if downloaded_genomes:
        with open(args.output_dir + "README_dataset_summary.txt", "w") as fd:
            fd.write("{:<20}{:<30}{:<20}\n".format("Dataset Number:", "Species:", "Number of Genomes:"))
            for summary in dataset_summaries:
                fd.write("{:<20d}{:<30}{:<20d}\n".format(summary[0], summary[1], summary[2]))

def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script allows users to download certain groups of genomes from the NCBI databases. It "
                                                  "will prompt users for species that you would like to download, and then those genomes "
                                                  "will be saved in the provided folder.")
    parser.add_argument("-r", "--refseq", action="store_true", default=False, help="download genomes from RefSeq")
    parser.add_argument("-g", "--genbank", action="store_true", default=False, help="download genomes from GenBank")
    parser.add_argument("-o", dest="output_dir", required=True, help="directory where all the genomes will be stored.")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Makes sure the arguments provided by user are valid, and make sense """
    if not args.refseq and not args.genbank:
        print("Error: Need to specify where to download genomes from either RefSeq (-r) or GenBank (-g)")
        exit(1)
    if args.refseq and args.genbank:
        print("Error: Can only choose one database to download from.")
        exit(1)
    if not os.path.isdir(args.output_dir):
        print("Error: The provided output directory is not valid.")
        exit(1)
    elif args.output_dir[-1] != "/":
        args.output_dir += "/"

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)