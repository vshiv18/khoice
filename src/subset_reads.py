
import argparse
import os
import random

def main(args):
    input = open(args.input, "r")
    output = open(args.output,"w+")
    input_lines = [x.strip() for x in input.readlines()]
    entry_sum = 0
    entry_type = ""

    if args.kmers: 
        entry_type="kmers"

        # Build a dictionary of reads to sample from
        input_dict = {}
        for i in range(0,len(input_lines),2):
            if ">" in input_lines[i]:
                input_dict[input_lines[i]] = input_lines[i+1]

        # Chose and remove random read until reach kmer ceiling
        batch_size = 10000
        while entry_sum < args.num_included and len(input_dict) > 0:
            # Adjust batch size if it is larger then remaining number of reads
            if len(input_dict) < batch_size:
                batch_size = len(input_dict)
            
            # Sample a batch of reads from dictionary
            rand_entries = random.sample(list(input_dict.keys()), batch_size)
            curr_batch_id = 0

            # Extract reads one by one, write to output file, and keep track of how many kmers we have
            while entry_sum < args.num_included and len(input_dict) > 0 and curr_batch_id < batch_size:
                rand_seq = input_dict.pop(rand_entries[curr_batch_id],None)
                output.write(rand_entries[curr_batch_id]+"\n"+rand_seq+"\n")
                entry_sum += len(rand_seq) - args.k + 1
                curr_batch_id += 1

    elif args.half_mems:
        entry_type="half mems"
        input_dict = {}
        for i in range(0,len(input_lines),4):
            if(">" in input_lines[i]):
                # Lazy way
                input_dict[input_lines[i]] = [input_lines[i+1], input_lines[i+2], input_lines[i+3]]
        batch_size = args.num_included
        print("Dictionary Built...")
        if(len(input_dict) < args.num_included):
            batch_size = len(input_dict)
        rand_entries = random.sample (list(input_dict.keys()), batch_size)
        for entry in rand_entries: 
            entry_sum+=1
            output.write(entry+"\n")
            for i in range(3):
                output.write(input_dict[entry][i]+"\n")
    elif(args.mems):
        entry_type="bp of mems"
        input_dict = {}
        for i in range(0,len(input_lines),4):
            if(">" in input_lines[i]):
                # Lazy way
                input_dict[input_lines[i]] = [input_lines[i+1], input_lines[i+2], input_lines[i+3]]
        batch_size = args.num_included
        print("Dictionary Built...")

        batch_size = 10000
        while entry_sum < args.num_included and len(input_dict) > 0:
            if len(input_dict) < batch_size:
                batch_size = len(input_dict)
            rand_entries = random.sample (list(input_dict.keys()), batch_size)
            curr_batch_id = 0
            while entry_sum < args.num_included and curr_batch_id < batch_size:
                rand_entry = rand_entries[curr_batch_id]
                rand_seqs = input_dict.pop(rand_entry,None)
                output.write(rand_entry+"\n")
                for i in range(3):
                    output.write(rand_seqs[i]+"\n")
                
                entry_sum += len(rand_seqs[0])
                curr_batch_id += 1

    print("[log] Returned {entry_sum} {type} with requested {requested}".format(entry_sum = entry_sum, type = entry_type, requested = args.num_included))
    
    # If entries requested greater than input
    if len(input_dict) == 0 or entry_sum == len(input_dict):
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
    if not args.k == None and args.k < 1:
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