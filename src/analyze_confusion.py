from distutils.command.config import config
import os
import argparse
import csv

k_values = [str(x) for x in range(7, 23, 1)] + [str(x) for x in range(23, 37, 2)] + [str(x) for x in range(38,53,3)]

def main(args):
    accuracies = [[0 for i in range(args.num_datasets +1)] for j in range(28)]
    all_values = []
    k_index = 0
    num_datasets = args.num_datasets
    
    for k in k_values:
        curr_matrix_path = "{matrix_dir}/k_{k}_confusion_matrix.csv".format(matrix_dir = args.matrix_dir, k = k)
        curr_matrix = []
        with open(curr_matrix_path, 'r') as fp:
            lines = fp.readlines()
            for line in lines:
                entries = [float(x) for x in line.strip().split(",")]
                curr_matrix.append(entries)
        for pivot in range(num_datasets):
            tp = curr_matrix[pivot][pivot]
            fp = fn = tn = 0
            for row in range(num_datasets):
                for column in range(num_datasets + 1):
                    curr = curr_matrix[row][column]
                    if(column == pivot and row != pivot):
                        fp += curr
                    elif(row == pivot and column != pivot):
                        fn += curr
                    elif(row != pivot):
                        tn += curr
            all_values.append([k,pivot,tp,tn,fp,fn])
        k_index+=1
    if(args.short):
        output_acc = args.output_dir + "short_accuracy_values.csv"
    elif(args.long):
        output_acc = args.output_dir + "long_accuracy_values.csv"
    
    with open(output_acc, "w+") as csvfile:
        writer = csv.writer(csvfile)
        for row in all_values:
            writer.writerow(row)

def parse_arguments():
    parser = argparse.ArgumentParser(description="This script calculates precision and recall from a set of confusion matricies")
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