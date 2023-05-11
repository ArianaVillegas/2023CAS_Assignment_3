# Adapted from Project 2
import numpy as np


table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }




class BasicNeutralNetwork:
    def __init__(self, n_mutations) -> None:
        self.n_mutations = n_mutations
        self.mutations = {}
        self._build()
        
    def _build(self):
        for codon in table:
            syn = []
            nonsyn = []
            for cmp_codon in table:
                if sum([(codon[i] != cmp_codon[i]) for i in range(3)]) == self.n_mutations:
                    if table[codon] == table[cmp_codon]:
                        syn.append(cmp_codon)
                    else:
                        nonsyn.append(cmp_codon)
            self.mutations[codon] = (syn, nonsyn)
            
    def get_mutation(self, codon):
        return self.mutations[codon]
            
    def get_mutations(self):
        return self.mutations



def make_trans_table(num_hops, aa2idx):
    # Change number from 1, 2, or 3, for respective n-point strict mutation; leave as 1 for non-strict n-point mutation
    base_muts = BasicNeutralNetwork(num_hops)
    count_table = np.zeros((len(aa2idx), len(aa2idx)))

    for codonA, muts in base_muts.get_mutations().items():
        amino = table[codonA]
        row = aa2idx[amino]
        count_table[row][row] += len(muts[0])

        for codonB in muts[1]:
            amino = table[codonB]
            col = aa2idx[amino]
            count_table[row][col] += 1

    trans_table = count_table/np.sum(count_table, axis=0)
    return trans_table



def make_probs_table(seq, meta):
    cnt = {key: 0 for key in list(set(table.keys()))}
    for _, row in meta.iterrows():
        for i in range(row['start'], row['end'], 3):
            cnt[seq[i:i+3]] += 1
    aa = {key: 0 for key in list(set(table.values()))}
    for tri in cnt:
        aa[table[tri]] += cnt[tri]
    return aa, cnt



def get_mut_table(seq, meta, nextstrain):
    aa_muts = {key: [] for key in list(set(table.values()))}
    aa_cnt = {key: 0 for key in list(set(table.values()))}
    aa_entropy = {key: 0 for key in list(set(table.values()))}
    meta.index = meta['gene']
    for _, row in nextstrain.iterrows():
        gene = meta.loc[row['gene']]
        start = gene['start'] + row['position']*3
        amino = seq[start:start+3]
        aa_cnt[table[amino]] += 1
        aa_muts[table[amino]].append(row['events'])
        aa_entropy[table[amino]] += row['entropy']
    aa_entropy = {key: round(aa_entropy[key] / aa_cnt[key], 3) for key in list(set(table.values()))}
    return aa_muts, aa_entropy