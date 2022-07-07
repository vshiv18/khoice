# Name: merge_lists.py
# Description: This is a python script is used by experiment 5 in order to
#              read and analyze the SAM file generated by ri-index locate
#              to generate a confusion matrix
#
# Date: June 14th, 2022

from distutils.command.config import config
import os
import argparse
import csv
import pysam

def main(args):
    # Set up confusion matrix and reference list
    confusion_matrix = [[0 for i in range(args.num_datasets)] for j in range(args.num_datasets)]
    refs_list = build_refs_list(args.ref_dir, args.num_datasets)
    for i in range(args.num_datasets):
        
        # Selects individual sam file to build map of read matches
        sam_file = pysam.AlignmentFile(args.sam_dir+"pivot_{n}.sam".format(n = i + 1), "r")
        
        print("\n[log] building a dictionary of the read alignments for pivot {n}".format(n = i + 1))
        read_mappings = {}

        for read in sam_file.fetch():
            dataset = find_class_of_reference_name(refs_list, read.reference_name)
            query_length = int(read.query_name.split("_")[3])
            # Only uses reads above threshold
            if query_length >= args.t:
                if read.query_name not in read_mappings:   
                    read_mappings[read.query_name] = [query_length]
                read_mappings[read.query_name].append(dataset)
        sam_file.close()
        print(read_mappings)
        for key in read_mappings:
            mem_len = read_mappings[key][0]
            curr_set = set(read_mappings[key][1:])
            for dataset in curr_set:
                confusion_matrix[i][dataset -1] += 1/len(curr_set) * mem_len 
        
    
    # Create output file
    output_matrix = args.output_dir + "confusion_matrix.csv"
    #output_precision = args.output_dir + "precision.csv"
    
    # Write confusion matrix to output csv
    with open(output_matrix,"w+") as csvfile:
        writer = csv.writer(csvfile)
        for row in confusion_matrix:
            writer.writerow(row)


def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script helps to analyze the SAM file from experiment 5"
                                                 "in order to form a confusion matrix.")
    parser.add_argument("-n", "--num", dest="num_datasets", required=True, help="number of datasets in this experiment", type=int)
    parser.add_argument("-s", "--sam_file", dest="sam_dir", required=True, help="path to directory with SAM files to be analyzed")
    parser.add_argument("-r", "--ref_lists", dest="ref_dir", required=True, help="path to directory with dataset reference lists")
    parser.add_argument("-o", "--output_path", dest = "output_dir", required=True, help="path to directory for output matrix and accuracies")
    parser.add_argument("--half_mems", action="store_true", default=False, dest="half_mems", help="sam corresponds to half-mems")
    parser.add_argument("--mems", action="store_true", default=False, dest="mems", help="sam corresponds to mems (it can either be mems or half-mems, not both)")
    parser.add_argument("-t", "--threshold", dest = "t", required=False, default=0, help="optional threshold value for experiment 8")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Checks for invalid arguments """
    if(not os.path.isdir(args.sam_dir)):
        print("Error: sam directory does not exist.")
        exit(1)
    if(not os.path.isdir(args.ref_dir)):
        print("Error: reference list directory does not exist.")
        exit(1)
    if(not os.path.isdir(args.output_dir)):
        print("Error: output directory does not exist.")
        exit(1)
    if(args.num_datasets < 1):
        print("Error: number of datasets must be a positive integer")
    # Verify only one of the options are chosen: half-mems or mems
    if (args.half_mems and args.mems) or (not args.half_mems and not args.mems):
        print("Error: exactly one type needs to be chosen (--half-mems or --mems)")
        exit(1)

def build_refs_list(ref_lists_dir, num_datasets):
    """ Creates a list of sets containing reference names for each dataset """
    ref_list = []
    for i in range(1,num_datasets + 1):
        curr_file = ref_lists_dir+"dataset_{num}_references.txt".format(num = i)
        curr_set = set()
        with open(curr_file, "r") as input_fd:
            all_lines = [x.strip() for x in input_fd.readlines()]
            for line in all_lines:
                curr_set.add(line)
                curr_set.add("revcomp_"+line)
        ref_list.append(curr_set)
    return ref_list

def find_class_of_reference_name(ref_list, ref_name):
    """ Finds the class that a particular reference name occurs in """
    datasets = []
    for i, name_set in enumerate(ref_list):
        if ref_name in name_set:
            datasets.append(i+1)
    if len(datasets) != 1:
        print(f"Error: this read hits {len(datasets)} which is not expected.")
        exit(1)
    return datasets[0]

#def calculate_accuracy(confusion_matrix, num_datasets):
    """ Calculates accuracy score given confusion matrix """
    accuracies = []
    for pivot in range(num_datasets):
        tp = confusion_matrix[pivot][pivot]
        fp = fn = tn = 0
        for row in range(num_datasets):
            for column in range(num_datasets):
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