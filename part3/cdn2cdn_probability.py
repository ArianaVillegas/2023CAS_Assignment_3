from neutral_network import table, BasicNeutralNetwork
import numpy as np
import pandas as pd

sim_count = 0
total_cdns = 0

#For basic neutral network codons:
base_pairs = BasicNeutralNetwork(1)
base_pairs._build()

bps_list = list(base_pairs.mutations.keys())
print(bps_list)
print("\n")

size = len(bps_list)
total_cdns = size

#cdn2cdn_table = [[0 for i in range(size)] for j in range(size)]
cdn2cdn_table = np.zeros((total_cdns, total_cdns))

firstpos = None
secondpos = None

for key in base_pairs.mutations.keys():
    if firstpos == None:
        firstpos = 0
    else:
        firstpos += 1

    for key2 in base_pairs.mutations.keys():
        if secondpos == None:
            secondpos = 0
        elif secondpos == (total_cdns-1):
            secondpos = 0
        else:
            secondpos += 1

        if key[0] == key2[0]:
            sim_count += 1
        if key[1] == key2[1]:
            sim_count += 1
        if key[2] == key2[2]:
            sim_count += 1

        if sim_count >= 2:
            cdn2cdn_table[firstpos][secondpos] = (1/total_cdns)
        else:
            cdn2cdn_table[firstpos][secondpos] = 0

        sim_count = 0

for i, row in enumerate(cdn2cdn_table):
    print(f"{bps_list[i]} ", end="")
    for val in row:
        print(f"{val:.6f} ", end="")
    print()
print("\n")

#For specifc genome codons:
inpath = "../part2/aamut_fitness_all.csv"
gene = "S"
df = pd.read_csv(inpath)
df["keep"] = df.apply(lambda row: row["gene"] == gene, axis=1)
df = df[df["keep"]]

'''
...not sure what comes next after this inpath.
...will we also be using the sequence txt file?
'''
