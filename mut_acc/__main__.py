# Adapted from Project 2
from mut_acc.neutral_network import table, get_mut_table
from mut_acc.neutral_network import make_trans_table, make_probs_table
from mut_acc.utils import read_seq, read_sequences, plot_dist, plot_cmp
from scipy.stats import entropy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

plt.style.use('dark_background')

# Set global values
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 6

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "../data")
TABLE_DIR = os.path.join(BASE_DIR, "tables")
FIGURE_DIR = os.path.join(BASE_DIR, "figures")
STRICT = False
if not os.path.exists(TABLE_DIR):
   os.makedirs(TABLE_DIR)
if not os.path.exists(FIGURE_DIR):
   os.makedirs(FIGURE_DIR)
   

# Read data
genome_influenza = read_seq(os.path.join(DATA_DIR, "genome_influenza.txt"))
meta_influenza = pd.read_csv(os.path.join(DATA_DIR, "meta_influenza.csv"))
nextstrain_influenza = pd.read_csv(os.path.join(DATA_DIR, "nextstrain_influenza.tsv"), sep='\t')
print('len influenza', len(genome_influenza))

genome_covid = read_seq(os.path.join(DATA_DIR, "genome_covid.txt"))
meta_covid = pd.read_csv(os.path.join(DATA_DIR, "meta_covid.csv"))
nextstrain_covid = pd.read_csv(os.path.join(DATA_DIR, "nextstrain_covid.tsv"), sep='\t')
print('len covid', len(genome_covid))


# Set index of aminoacids and triplets
idx2aa = list(set(table.values()))
idx2aa.sort()
aa2idx = {aa: i for i, aa in enumerate(idx2aa)}

idx2tri = list(set(table.keys()))
idx2tri.sort()
tri2idx = {aa: i for i, aa in enumerate(idx2tri)}


# Build probs table of aminoacids
probs_influenza_aa, probs_influenza_tri = make_probs_table(genome_influenza, meta_influenza)
probs_influenza_aa = [round(probs_influenza_aa[key]/sum(probs_influenza_aa.values()), 4) for key in idx2aa]
probs_influenza_tri = [round(probs_influenza_tri[key]/sum(probs_influenza_tri.values()), 4) for key in idx2tri]
plot_dist(idx2aa, probs_influenza_aa, "Aminoacid distribution - Influenza", "Aminoacid", 
          os.path.join(FIGURE_DIR, "influenza_aa_dist.png"), False)
plot_dist(idx2tri, probs_influenza_tri, "Triplet distribution - Influenza", "Triplet", 
          os.path.join(FIGURE_DIR, "influenza_tri_dist.png"), True)

probs_covid_aa, probs_covid_tri = make_probs_table(genome_covid, meta_covid)
probs_covid_aa = [round(probs_covid_aa[key]/sum(probs_covid_aa.values()), 4) for key in idx2aa]
probs_covid_tri = [round(probs_covid_tri[key]/sum(probs_covid_tri.values()), 4) for key in idx2tri]
plot_dist(idx2aa, probs_covid_aa, "Aminoacid distribution - SARS-CoV-2", "Aminoacid", 
          os.path.join(FIGURE_DIR, "covid_aa_dist.png"), False)
plot_dist(idx2tri, probs_covid_tri, "Triplet distribution - SARS-CoV-2", "Triplet", 
          os.path.join(FIGURE_DIR, "covid_tri_dist.png"), True)


# Build single point mutation table
trans_table = make_trans_table(1, aa2idx).T
entropy_table = np.asarray([entropy(row[row != 0], base=2) for row in trans_table])


# Mutations per aminoacid our method
trans_influenza = (trans_table * probs_influenza_aa).T
trans_influenza = np.round(trans_influenza / np.sum(trans_influenza, axis=0), 3)
print('Influenza')
print(trans_influenza.diagonal())
print(np.mean(trans_influenza.diagonal()), (np.sum(trans_influenza) - np.sum(trans_influenza.diagonal()))/len(trans_influenza))

trans_influenza_entropy = (entropy_table * probs_influenza_aa)

trans_covid = (trans_table * probs_covid_aa).T
trans_covid = np.round(trans_covid / np.sum(trans_covid, axis=0), 3)
print('Covid')
print(trans_covid.diagonal())
print(np.mean(trans_covid.diagonal()), (np.sum(trans_covid) - np.sum(trans_covid.diagonal()))/len(trans_covid))

trans_covid_entropy = (entropy_table * probs_covid_aa)


print('NEXT STRAIN')


# Mutations per aminoacid nextstrain
muts_influenza, entropy_influenza = get_mut_table(genome_influenza, meta_influenza, nextstrain_influenza)
muts_influenza = dict(sorted(muts_influenza.items()))
entropy_influenza = dict(sorted(entropy_influenza.items()))
fig, ax = plt.subplots(figsize=(13, 10))
ax.boxplot(muts_influenza.values())
ax.set_xticklabels(muts_influenza.keys())
# ax.set_yscale('log')
ax.set_xlabel('Amino acid')
ax.set_ylabel('Number of mutations')
plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, "influenza_next.png"), transparent=True)

entropy_influenza = np.array(list(entropy_influenza.values())) * 60
print('entropy_influenza\n', np.round(entropy_influenza * probs_influenza_aa, 3))
print(np.round(np.mean(entropy_influenza * probs_influenza_aa), 3))

muts_covid, entropy_covid = get_mut_table(genome_covid, meta_covid, nextstrain_covid)
muts_covid = dict(sorted(muts_covid.items()))
entropy_covid = dict(sorted(entropy_covid.items()))
fig, ax = plt.subplots(figsize=(13, 10))
ax.boxplot(muts_covid.values())
ax.set_xticklabels(muts_covid.keys())
# ax.set_yscale('log')
ax.set_xlabel('Amino acid')
ax.set_ylabel('Number of mutations')
plt.tight_layout()
plt.savefig(os.path.join(FIGURE_DIR, "covid_next.png"), transparent=True)

entropy_covid = np.array(list(entropy_covid.values())) * 120
print('entropy_covid\n', np.round(entropy_covid * probs_covid_aa, 3))
print(np.round(np.mean(entropy_covid * probs_covid_aa), 3))


# Save compare plots
plot_cmp(idx2aa, entropy_influenza, probs_influenza_aa, trans_influenza, os.path.join(FIGURE_DIR, "influenza_entropy.png"))
plot_cmp(idx2aa, entropy_covid, probs_covid_aa, trans_covid, os.path.join(FIGURE_DIR, "covid_entropy.png"))