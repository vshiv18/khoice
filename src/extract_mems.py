#!/usr/bin/env python3

# Name: extract_mems.py
# Description: This script is meant to take a FASTA file for which
#              matching statistics have been computed for and then
#              extract either the half-mems or mems and output
#              them to stdout.
# Date: June 21, 2022
# Adapted by Vikram Shivakumar for Sigmoni

import argparse
import os
from Bio import SeqIO
import numpy as np

def print_out_half_mems(curr_seq, lengths_array, curr_seq_id, threshold, out_file, read_num):
    """ Given a pattern and lengths array (matching statistics), print out half mems in FASTQ file """
    for i in range(len(curr_seq)):
        curr_half_mem_length = lengths_array[i]

        # Filter out half-mems below threshold
        if curr_half_mem_length >= threshold:
            write_mem_length = curr_half_mem_length

            # Limit the length that is written to save space
            if curr_half_mem_length >= 1000:
                write_mem_length = 1000

            random_quality = "#" * write_mem_length
            curr_half_mem = curr_seq[i:i+write_mem_length]

            out_file.write(f">read_{read_num}_halfmem_{curr_seq_id}_length_{curr_half_mem_length}\n{curr_half_mem}\n+\n{random_quality}\n")
            curr_seq_id += 1
    return curr_seq_id

def print_out_mems(curr_seq, lengths_array, curr_seq_id, threshold, out_file, read_num):
    """ Given a pattern and lengths array (matching statistics), print out MEMs in FASTQ file """

    # Print out MEM at position 0 (if it is above threshold)
    curr_mem_length = lengths_array[0]
    if curr_mem_length >= threshold:
        write_mem_length = curr_mem_length
        
        # Added to reduce size of output FASTQ (mostly for out-genome case)
        if curr_mem_length >= 1000:
            write_mem_length = 1000

        random_quality = "#" * write_mem_length
        curr_mem = curr_seq[0:write_mem_length]
        out_file.write(f">read_{read_num}_mem_{curr_seq_id}_length_{curr_mem_length}\n{curr_mem}\n+\n{random_quality}\n")
        curr_seq_id += 1

    # Scan MS for peaks which correspond to MEMs
    for i in range(1, len(curr_seq)):

        prev_mem_length = lengths_array[i - 1]
        curr_mem_length = lengths_array[i]

        if curr_mem_length >= threshold:
            # Take only MEMs
            if curr_mem_length >= prev_mem_length:
                write_mem_length = curr_mem_length

                # Limit mem length for memory
                if curr_mem_length >= 1000:
                    write_mem_length = 1000
                    
                random_quality = "#" * write_mem_length
                curr_mem = curr_seq[i:i+write_mem_length]
                out_file.write(f">read_{read_num}_mem_{curr_seq_id}_length_{curr_mem_length}\n{curr_mem}\n+\n{random_quality}\n")
                curr_seq_id += 1
    return curr_seq_id

def parse_ms(fname, names=None):
    with open(fname,'r') as f:
        if names:
            return {n : np.fromstring(x, dtype=int, sep=' ').tolist() for n, (_, x) in zip(open(names, 'r').read().splitlines(), 
                                                                    SeqIO.FastaIO.FastaTwoLineParser(f))}
        return {i : np.fromstring(x, dtype=int, sep=' ').tolist() for i, x in SeqIO.FastaIO.FastaTwoLineParser(f)}

def parse_fasta(fname):
    with open(fname,'r') as f:
        return {i : x for i, x in SeqIO.FastaIO.FastaTwoLineParser(f)}

def main(args):
    """ Extracts either the half-mems/mems from a set of patterns matched to database """

    # Open the reads, the MS lengths, and the output FASTQ
    # pattern_fd = open(args.pattern_file, "r")
    # length_fd = open(args.length_file, "r")
    patterns = parse_fasta(args.pattern_file)
    lengths = parse_ms(args.length_file)
    out_fd = open(args.output_file, "w")

    # Verify the number of sequences in each file is the same
    assert len(patterns) == len(lengths), "assertion failed: mismatch in the number of sequences in each file."
    
    curr_id = 0
    seq_id = 0
    for r in patterns.keys():
        lengths_array = lengths[r]
        curr_seq = patterns[r]
        assert len(lengths_array) == len(curr_seq), "assertion failed: mis-match in terms of the lengths of lengths array"
        if args.half_mems: # deal with half-MEMs
            curr_id = print_out_half_mems(curr_seq, lengths_array, curr_id, args.threshold, out_fd, seq_id)
        else: # deal with MEMs
            curr_id = print_out_mems(curr_seq, lengths_array, curr_id, args.threshold, out_fd, seq_id)
        seq_id += 1
    # Start to extract the requested output (half-mems or mems) and output to stdout

    # Close input/output
    out_fd.close()

def parse_arguments():
    """ Parses the command-line arguments, and returns them """
    parser = argparse.ArgumentParser(description="Takes in matching statistic lengths, and extracts either half-mems or mems, prints to stdout.")
    parser.add_argument("-l", dest="length_file", help="path to matching statistic length file", required=True)
    parser.add_argument("-p", dest="pattern_file", help="path to pattern file that correspond to the matching statistics.", required=True)
    parser.add_argument("--half_mems", action="store_true", default=False, dest="half_mems", help="extract half-mems and print to stdout")
    parser.add_argument("--mems", action="store_true", default=False, dest="mems", help="extract mems and print to stdout (it can either be mems or half-mems, not both)")
    parser.add_argument("-t", dest="threshold", default=0, type=int, help="threshold for half-mem length to use")
    parser.add_argument("-o", dest="output_file", help="path to output file", required=True)
    args = parser.parse_args()
    return args

def check_args(args):
    """ Checks the command-line arguments and makes sure they are valid """

    # Verify the files exist
    if not os.path.isfile(args.length_file):
        print("Error: this matching statistic length file does not exist.")
        exit(1)
    if not os.path.isfile(args.pattern_file):
        print("Error: this pattern file does not exist.")
        exit(1)
    
    # Verify the files are the right type
    #if not args.length_file.endswith(".lengths"):
    #    print("Error: the lengths file has the incorrect file extension.")
    #    exit(1)
    if not args.pattern_file.endswith(".fa") and not args.pattern_file.endswith(".fna"):
        print("Error: the pattern file has the incorrect file extension.")
        exit(1)
    
    # Verify only one of the options are chosen: half-mems or mems
    if (args.half_mems and args.mems) or (not args.half_mems and not args.mems):
        print("Error: exactly one type needs to be chosen (--half-mems or --mems)")
        exit(1)
    
    # Verify if using half-mems, that we specified a threshold
    if args.half_mems and args.threshold == 0:
        print("Error: need to specify a threshold if using half-mems")
        exit(1)
    if args.half_mems and args.threshold < 0:
        print("Error: need to specify a positive value as the threshold")
        exit(1)

if __name__ == "__main__":
    args = parse_arguments()
    check_args(args)
    main(args)
