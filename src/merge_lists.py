#!/usr/bin/env python3

# Name: merge_lists.py
# Description: This is a python script is used by experiment 4 in order to merge
#              the kmer lists from kmc3 in order to build a confusion matrix.
#
# Date: June 14th, 2022

import argparse
import os

def main(args):
    """ main method for the script """

    # Read in all the paths into these two lists
    pivot_files = []
    intersect_files = [] 
    with open(args.pivot_filelist, "r") as input_fd:
        pivot_files = [x.strip() for x in input_fd.readlines()]
    with open(args.intersect_list, "r") as input_fd:
        intersect_files = [x.strip() for x in input_fd.readlines()]

    # Verify each of the files in the file_lists is valid
    for file_path in (pivot_files + intersect_files):
        if not os.path.isfile(file_path):
            print(f"Error: At least one of the file paths in the file lists is not valid ({file_path})")

    print("[log] all file paths are valid, now will start to merge the lists.")

    # Define the confusion matrix
    num_datasets = args.num_datasets
    confusion_matrix = [[0 for i in range(num_datasets)] for j in range(num_datasets)]

    # Fill in the confusion matrix from the upper-left to bottom right
    curr_pivot_num = 0
    for curr_pivot_file in pivot_files:

        ################################################################################
        # TODO - Build a dictionary holding the kmers and counts for the current pivot
        # You can insert a method here to keep the main method nice and consice
        ################################################################################

        curr_intersect_num = 0
        for curr_intersect_file in intersect_files:

            ##############################################################################
            # TODO - Go though all the kmers in the intersect, and update the dictionary
            # This can be a method as well if you want.
            ##############################################################################

            curr_intersect_num += 1 # might to be last line of inner for-loop

        # Fill in the current row of the confusion matrix
        
        ############################################################
        # TODO: Take the kmer counts, and datasets they occur in
        #       and fill in the confusion matrix.
        ############################################################

        curr_pivot_num += 1 # might to be last line of for-loop
    
    # Print out your confusion matrix to a csv file

def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script helps to analyze the data in experiment 4 in order to"
                                                 "merge the kmer-lists to form a confusion matrix.")
    parser.add_argument("-n", "--num", dest="num_datasets", required=True, help="number of datasets in this experiment", type=int)
    parser.add_argument("-p", "--pivot_list", dest="pivot_filelist", required=True, help="path to text file with a list of pivot kmer lists")
    parser.add_argument("-i", "--intersect_list", dest="intersect_list", required=True, help="path to text file with a list of kmer lists representing the intersections")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Checks the arguments to make sure they are valid """
    if not os.path.isfile(args.pivot_filelist) or not os.path.isfile(args.intersect_list):
        print("Error: One of the provided files is not valid.")
        exit(1)
    if args.num_datasets <= 0:
        print("Error: The number of datasets needs to be positive integer.")
        exit(1)

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)
