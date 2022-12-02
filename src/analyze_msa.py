
import sys
import os
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt

def get_entropy_scores(line_set):
    aln_results = []
    for seq in line_set:
        if "*" not in seq:
            aln_results.append(seq.split()[1])
    assert all(len(x) == len(aln_results[0]) for x in aln_results), "all alignments should be equal length"

    length = len(aln_results[0])
    results = []
    for i in range(length):
        # Get all the values in a column of MSA
        entries = [x[i] for x in aln_results]

        # Get the probability of each element, and compute entropy
        values, value_counts = np.unique(entries, return_counts=True)
        probs = value_counts/len(entries)

        entropy = 0.0
        for i in probs:
            entropy -= i * math.log(i)

        # Store the entroy of each column
        results.append(entropy)
    return results

### Main Method Code
msa_file = sys.argv[1]
with open(msa_file, "r") as in_fd:
    all_msa_lines = in_fd.readlines()

curr_line_set = []
final_results = []
for i, line in enumerate(all_msa_lines[3:]):
    if len(line) > 1 and len(line.split()) == 2: # alignment line
        curr_line_set.append(line.strip())
    elif len(curr_line_set) > 0:
        final_results += get_entropy_scores(curr_line_set)
        curr_line_set = []

### Plotting
x = [x for x in range(1, len(final_results)+1)]
plt.bar(x, height=final_results, width=1.0)
plt.xlabel("Base Position in Protein Sequence")
plt.ylabel("Shannon Entropy")
plt.savefig(msa_file + ".png", dpi=500, bbox_inches="tight")

plt.figure()
rolling_avg = np.convolve(final_results, np.ones(250)/250, mode='valid')
x = [x for x in range(1, len(rolling_avg)+1)]
plt.bar(x, height=rolling_avg, width=1.0)
plt.xlabel("Base Position in Protein Sequence")
plt.ylabel("Avg Shannon Entropy (over 250 aa windows)")
plt.savefig(msa_file + ".rolling.png", dpi=500, bbox_inches="tight")




