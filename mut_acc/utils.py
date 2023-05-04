import matplotlib.pyplot as plt
import numpy as np
import os



def read_seq(path):
    f = open(path)
    lines = f.readlines()
    lines = [line[:-1] for line in lines]
    seq = ''.join(lines)
    return seq



def read_sequences(path):
    f = open(path)
    lines = f.readlines()
    names = []
    seq = []
    len_ = []
    for i in range(0, len(lines), 2):
        names.append(lines[i][1:])
        seq.append(lines[i+1])
        len_.append(len(lines[i+1]))
    return [names, seq, len_]



def write_trans_table(path, table, dtype, caption, label, idx2aa):
    with open(path, "w") as file:
        file.write("\\begin{table*}[]\n\t\centering\n\t\\begin{tabular}")
        file.write("{c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c}\n")

        # Write the headers
        file.write("\t\t")
        for i in range(len(table)):
            if idx2aa[i] == "_":
                file.write(" & \\_")
            else:
                file.write(f" & {idx2aa[i]}")
        file.write(" \\\\\n")

        for i, row in enumerate(table):
            file.write("\t\t\\hline\n")
            if idx2aa[i] == "_":
                file.write("\t\t\\_")
            else:
                file.write(f"\t\t{idx2aa[i]}")
            for val in row:
                if dtype == int:
                    file.write(f" & {int(val):2d}")
                elif dtype == float:
                    file.write(f" & {val:.2f}")
            file.write(" \\\\\n")

        file.write("\t\\end{tabular}\n")
        file.write("\t\\caption{" + caption + "}\n")
        file.write("\t\\label{tab:" + label + "}\n")
        file.write("\\end{table*}\n")
        


def plot_dist(idx2aa, dist, title, xlabel, path, rotate=False):
    plt.rcParams.update({'font.size': 26})
    _, axs = plt.subplots(figsize=(12, 10))
    axs.bar(idx2aa, dist, color="white", edgecolor="black")
    axs.set_xlabel(xlabel)
    axs.set_ylabel("Percentage found in genome")
    # axs.set_title(title)
    if rotate:
        plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight", transparent=True)
    
    
    
def plot_cmp(idx2aa, entropy, probs_aa, trans, path):
    _, ax = plt.subplots(figsize=(13, 10))
    ax.bar(np.arange(len(idx2aa)) - 0.2, entropy * probs_aa, 0.4, label='NextStrain')
    ax.bar(np.arange(len(idx2aa)) + 0.2, trans.diagonal(), 0.4, label='Our method')
    ax.set_xticks(np.arange(len(idx2aa)), idx2aa)
    ax.set_xlabel('Amino acid')
    ax.set_ylabel('Ratio of synonymous mutation')
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, transparent=True)
    
    

def write_synnonsyn_table(idx2aa, path, table, dtype, caption, label):
    with open(path, "w") as file:
        file.write("\\begin{table}[]\n\t\centering\n\t\\begin{tabular}")
        file.write("{c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c}\n")

        # Write the headers
        file.write("\t\t& Synonymous & Non-synonymous \\\\\n")

        for i, row in enumerate(table):
            file.write("\t\t\\hline\n")
            if idx2aa[i] == "_":
                file.write("\t\t\\_")
            else:
                file.write(f"\t\t{idx2aa[i]}")
            for val in row:
                if dtype == int:
                    file.write(f" & {int(val):2d}")
                elif dtype == float:
                    file.write(f" & {val:.2f}")
            file.write(" \\\\\n")

        file.write("\t\\end{tabular}\n")
        file.write("\t\\caption{" + caption + "}\n")
        file.write("\t\\label{tab:" + label + "}\n")
        file.write("\\end{table}\n")