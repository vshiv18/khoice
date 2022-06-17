#!/usr/bin/env python3

# Name: merge_lists.py
# Description: This is a python script is used by experiment 4 in order to merge
#              the kmer lists from kmc3 in order to build a confusion matrix.
#
# Date: June 14th, 2022

import argparse
from calendar import c
import os
import csv
from tkinter.tix import COLUMN

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
            exit(1)

    print("[log] all file paths are valid, now will start to merge the lists.")

    # Define the confusion matrix
    # Last column counts kmers only in pivot
    num_datasets = args.num_datasets
    confusion_matrix = [[0 for i in range(num_datasets + 1)] for j in range(num_datasets)]

    # Fill in the confusion matrix from the upper-left to bottom right
    curr_pivot_num = 0
    curr_intersect_num = 0
    for curr_pivot_file in pivot_files:
        print("[log] processing pivot {pivot} with k = {k}".format(pivot = curr_pivot_num + 1, k = args.k))
        # Build dictionary of kmers with counts
        curr_pivot_dict = build_dictionary(curr_pivot_file)
        for i in range(num_datasets):
            # Update dictionary with intersect num
            curr_pivot_dict = update_dictionary(curr_pivot_dict,intersect_files[curr_intersect_num], curr_intersect_num, num_datasets)
            curr_intersect_num += 1 # might to be last line of inner for-loop
            
        unique_pivot_count = 0
        for kmer in curr_pivot_dict:
            # Sum counts of kmers only in pivot
            if(len(curr_pivot_dict[kmer]) == 1):
                unique_pivot_count += curr_pivot_dict[kmer][0]
        
        # Fill in the current row of the confusion matrix
        for kmer in curr_pivot_dict:
            count = curr_pivot_dict[kmer][0]
            matches = curr_pivot_dict[kmer][1:]
            for match in matches:
                confusion_matrix[curr_pivot_num][match] += 1 / len(matches) * count
        # Last column is pivot only
        confusion_matrix[curr_pivot_num][num_datasets] = unique_pivot_count
        curr_pivot_num += 1 

    # Print out confusion matrix to a csv file
    output_matrix = args.output_path + "confusion_matrix/k_"+ args.k +"_confusion_matrix.csv"
    with open(output_matrix,"w+") as csvfile:
        writer = csv.writer(csvfile)
        for row in confusion_matrix:
            writer.writerow(row)
    
    # Calculate accuracy and print to csv file (k, pivot_1, ... pivot_n)
    accuracies = [args.k]
    for score in calculate_accuracy(confusion_matrix,num_datasets):
        accuracies.append(score)
    output_acc = args.output_path + "accuracy/k_"+ args.k +"_accuracy.csv"
    with open(output_acc, "w+") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(accuracies)



def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script helps to analyze the data in experiment 4 in order to"
                                                 "merge the kmer-lists to form a confusion matrix.")
    parser.add_argument("-n", "--num", dest="num_datasets", required=True, help="number of datasets in this experiment", type=int)
    parser.add_argument("-p", "--pivot_list", dest="pivot_filelist", required=True, help="path to text file with a list of pivot kmer lists")
    parser.add_argument("-i", "--intersect_list", dest="intersect_list", required=True, help="path to text file with a list of kmer lists representing the intersections")
    parser.add_argument("-o", "--output_path", dest = "output_path", required=True, help = "path to output directory") # Match snakemake format
    parser.add_argument("-k", "--k_value", dest = "k", required = True, help = "k value for this experiment")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Checks the arguments to make sure they are valid """
    check_path(args.pivot_filelist)
    check_path(args.intersect_list)
    if args.num_datasets <= 0:
        print("Error: The number of datasets needs to be positive integer.")
        exit(1)

def check_path(file):
    if(not os.path.isfile(file)):
        print("Error: One of the provided files is not valid: "+file)
        exit(1)

def build_dictionary(pivot_path):
    """ Builds dictionary of pivot kmers with dict[0] as count """
    dict = {}
    with open(pivot_path, 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            kmer = line.split()[0]
            count = line.split()[1]
            dict[kmer] = [int(count)]
    return dict

def update_dictionary(pivot_dict, intersect_path, intersect_num, num_datasets):
    """ Updates dictionary with intersect kmers """
    with open(intersect_path, 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            kmer = line.split()[0]
            if(kmer in pivot_dict):
                pivot_dict[kmer].append(intersect_num % num_datasets)
    return pivot_dict

def calculate_accuracy(confusion_matrix, num_datasets):
    """ Calculates accuracy score given confusion matrix """
    accuracies = []
    for pivot in range(num_datasets):
        tp = confusion_matrix[pivot][pivot]
        fp = fn = tn = 0
        for row in range(num_datasets):
            for column in range(num_datasets + 1):
                curr = confusion_matrix[row][column]
                if(column == pivot and row != pivot):
                    fp += curr
                elif(row == pivot and column != pivot):
                    fn += curr
                elif(row != pivot):
                    tn += curr
        accuracies.append((tp + tn)/(tp + tn + fp + fn))
    return accuracies

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)
