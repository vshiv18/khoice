
import argparse
import os
import random

def main(args):
    input = open(args.input, "r")
    output = open(args.output,"w+")
    input_lines = [x.strip() for x in input.readlines()]
    line_num = 0
    entry_sum = 0



    if(args.kmers): 
        # Write kmers until over requested num
        # Build dictionary of headers and reads
        input_dict = {}
        for i in range(0,len(input_lines),2):
            if(">" in input_lines[i]):
                input_dict[input_lines[i]] = input_lines[i+1]
                
        # Chose random read until reach kmer ceiling
        while(entry_sum < args.num_included and line_num < len(input_lines)):
            rand_entry = random.choice(list(input_dict.keys()))
            rand_seq = input_dict[rand_entry]
            output.write(rand_entry+"\n"+rand_seq+"\n")
            entry_sum += len(rand_seq) - args.k + 1
    elif(args.half_mems):
        print("Not implemented")
    elif(args.mems):
        print("Not implemented")

    print("Returned {entry_sum} kmers with requested {requested}".format(entry_sum = entry_sum, requested = args.num_included))
    # If entires requested greater than input
    if(line_num == len(input_lines)):
        print("Error: Number of entries requested >= input entries ({entry_sum})".format(entry_sum=entry_sum))
        exit(1)

def parse_arguments():
    """ Defines the command-line argument parser, and return arguments """
    parser = argparse.ArgumentParser(description="This script subsets a certain amount of entries (kmers, mems, half-mems) given a file of reads")
    parser.add_argument("-i", "--input", dest="input", required=True, help="path to fq/fastq file with all reads")
    parser.add_argument("-o", "--output", dest="output", required=True, help="path to output fasta/fastq file")
    parser.add_argument("-n", "--num_included", dest="num_included", required=True, help="number of kmers/half-mems/mems included", type = int)
    parser.add_argument("-k", "--k_value", dest="k", required=False, help="length of kmers", type = int)
    parser.add_argument("--kmers", action="store_true", default=False, dest="kmers", help="subset reads for kmer analysis")
    parser.add_argument("--half_mems", action="store_true", default=False, dest="half_mems", help="subset reads for mem analysis")
    parser.add_argument("--mems", action="store_true", default=False, dest="mems", help="subset reads for mem analysis")
    args = parser.parse_args()
    return args

def check_arguments(args):
    """ Checks for valid arguments """
    if not os.path.isfile(args.input):
        print("Error: this input fasta/fastq file does not exist.")
        exit(1)
    if args.num_included < 1:
        print("Error: number of requested entries must be greater than 0")
        exit(1)
    if args.k < 1:
        print("Error: k value must be greater than 0")
        exit(1)
    # Can only chose one option
    chosen = 0
    for option in [args.kmers, args.mems, args.half_mems]:
        if option: chosen +=1
    if(chosen != 1):
        print("Error: exactly one of the following options must be chosen: --kmers, --half_mems, --mems")
        exit(1)
    

if __name__ == "__main__":
    args = parse_arguments()
    check_arguments(args)
    main(args)