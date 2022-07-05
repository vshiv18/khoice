from distutils.command.config import config
import os
import argparse
import csv

k_values = [str(x) for x in range(7, 23, 1)] + [str(x) for x in range(23, 37, 2)] + [str(x) for x in range(38,53,3)]

def main(args):
    accuracies = [[0 for i in range(args.num_datasets +1)] for j in range(28)]
    k_index = 0
    for k in k_values:
        curr_matrix_path = "{matrix_dir}/k_{k}_confusion_matrix.csv".format(matrix_dir = args.matrix_dir, k = k)
        with open(curr_matrix_path, 'r') as fp:
            lines = fp.readlines()
            curr_row = 0
            
            for line in lines:
                
                entries = [float(x) for x in line.strip().split(",")]
                row_sum = sum(entries)
                true = entries[curr_row]
                accuracies[k_index][0] = k
                accuracies[k_index][curr_row + 1] = float(true)/row_sum
                curr_row += 1
        k_index+=1
    if(args.short):
        out_path = args.output_dir+"short_accuracies.csv"
    else:
        out_path = args.output_dir+"long_accuracies.csv"
    with open(out_path,"w+") as csvfile:
        writer = csv.writer(csvfile)
        for row in accuracies:
            writer.writerow(row)

def parse_arguments():
    parser = argparse.ArgumentParser(description="This script calculates accuracy from a set of confusion matricies")
    parser.add_argument("-n", "--num", dest="num_datasets", required=True, help="number of datasets in this experiment", type=int)
    parser.add_argument("-m", "--matrices", dest = "matrix_dir", required=True, help="directory containing confusion matrices")
    parser.add_argument("-o", "--output", dest = "output_dir", required=False, help="path to output directory")
    parser.add_argument("--long", action="store_true", default=False, dest="long", help="confusion matrices correspond to long reads")
    parser.add_argument("--short", action="store_true", default=False, dest="short", help="confusion matrices correspond to long reads")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()
    main(args)