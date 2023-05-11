from part3.neutral_network import BasicNeutralNetwork, table
import numpy as np
import scipy
import matplotlib.pyplot as plt
import time

# index to amino acid and vice versa
idx2aa = list(set(table.values()))
idx2aa.sort()
aa2idx = {aa: i for i, aa in enumerate(idx2aa)}
dtype = np.double

synnonsyn_count_table = np.zeros((len(idx2aa), 2), dtype=dtype)

base_muts = BasicNeutralNetwork(1)
base_muts._build()
for condon, muts in base_muts.mutations.items():
    amino = table[condon]
    row = aa2idx[amino]
    synnonsyn_count_table[row][0] += len(muts[0])
    synnonsyn_count_table[row][1] += len(muts[1])

synonsyn_trans_tablenp = (synnonsyn_count_table.T/np.sum(synnonsyn_count_table.T, axis=0, dtype=dtype)).T

def make_trans_table(num_hops):
    base_muts = BasicNeutralNetwork(num_hops)
    base_muts._build()

    count_table = np.zeros((len(idx2aa), len(idx2aa)), dtype=dtype)

    for condonA, muts in base_muts.mutations.items():
        amino = table[condonA]
        row = aa2idx[amino]
        count_table[row][row] += len(muts[0])

        for codonB in muts[1]:
            amino = table[codonB]
            col = aa2idx[amino]
            count_table[row][col] += 1

    trans_table = count_table/np.sum(count_table, axis=0, dtype=dtype)
    return trans_table

def print_trans_table(trans_table):
    for row in trans_table:
        for val in row:
            print(f"{val:6.4f} ", end="")
        print()

def make_trans_prop(prop):
    r"""
    Makes a transition table where the proportion of mutating genes to total
    genes is given by `prop`.
    """
    trans_table = make_trans_table(1)
    n, _ = trans_table.shape
    res = prop * trans_table + (1 - prop) * np.eye(n, dtype=dtype)
    for i in range(n):
        res[i] = res[i] / np.sum(res[i])
    return res


def experiment(prop, n_steps, state=None):
    trans_table = make_trans_prop(prop)
    n, _ = trans_table.shape

    if state is None:
        state = np.random.rand(n)
        state /= np.sum(state)

    dist = np.zeros((n_steps, n), dtype=dtype)
    nonsyn = np.zeros(n_steps, dtype=dtype)
    dist[0] = state
    nonsyn[0] = (state @ synonsyn_trans_tablenp)[1] # [syn, nonsyn]
    for i in range(1, n_steps):
        state = state @ trans_table
        state = state / np.sum(state) # Re-normalize because of numerical error
        dist[i] = state
        nonsyn[i] = (state @ synonsyn_trans_tablenp)[1] # [syn, nonsyn]
    return dist, nonsyn, trans_table

def savefig(ax, name, extension="png", do_save=True):
    if not do_save:
        return

    savefig.dir_name = "figures"
    if not os.path.exists(savefig.dir_name):
        os.path.makedir(savefig.dir_name)

    path = os.path.join(savefig.dir_name, name + "." + extension)
    plt.savefig(path, bbox_inches="tight", dpi=500, transparent=True)

def distance(vec1, vec2):
    return np.sqrt((vec1 - vec2) @ (vec1 - vec2))

def l2(vec):
    return distance(np.zeros(vec.shape), vec)

def plot_trans_table(trans_table, **params):
    figsize = params["figsize"]
    do_save = params["do_save"]
    extension = params["extension"]

    fig, ax = plt.subplots(figsize=figsize)
    img = ax.imshow(np.log(trans_table + 1e-7), aspect="equal", interpolation="none", cmap="plasma")
    fig.colorbar(img, ax=ax)
    ax.set_yticks(np.arange(len(idx2aa)), idx2aa)
    ax.set_xticks(np.arange(len(idx2aa)), idx2aa)
    ax.set_title("Transition table (log probability)")
    ax.set_xlabel("Amino acid (to)")
    ax.set_ylabel("Amino acid (from)")
    savefig(ax, "part2_table", extension, do_save)

def stdize(vecs):
    n, d = vecs.shape
    stdized = np.zeros((n, d))
    means = np.mean(vecs, axis=1)
    stds = np.std(vecs, axis=1)
    for feature in range(d):
        for sample in range(n):
            stdized[sample, feature] = \
                (vecs[sample, feature] - means[feature]) / stds[feature]
    return stdized

def pca(vecs, **params):
    n, d = vecs.shape
    stdized = stdize(vecs)
    cov = np.cov(stdized.T)
    assert cov.shape == (d, d)
    eigenvals, eigenvecs = np.linalg.eig(cov)

    sortedvals = np.zeros((n, 1))
    sortedvecs = np.zeros((n, d))

    print("Percent of variance explained:")
    percents = 100.0 * eigenvals / np.sum(eigenvals)
    for percent in percents:
        print(f"\t{percent:4.2f}")
    return eigenvecs

def eigen(trans_table, **params):
    r"""
    Gets the eigenvectors of the transition table, normalize them, and then
    generate the principal components

    Args:
        trans_table: the transition table
        **params: the dictionary that tracks our parameters

    Returns:
        eigenvectors, principal_components
    """
    eigenvals, eigenvecs = scipy.linalg.eig(trans_table, left=True, right=False)

    for eigenvec in eigenvecs:
        eigenvec = eigenvec / l2(eigenvec)

    pcs = pca(eigenvecs, **params)
    return eigenvecs, pcs

def plot_eigenvectors(eigenvecs, **params):
    figsize = params["figsize"]
    extension = params["extension"]
    do_save = params["do_save"]

    lim = np.max([np.abs(np.min(eigenvecs)), np.abs(np.max(eigenvecs))])
    fig, ax = plt.subplots(figsize=figsize)
    img = ax.imshow(eigenvecs.T, vmin=-lim, vmax=lim, cmap="seismic")
    fig.colorbar(img, ax=ax)
    ax.set_yticks(np.arange(len(idx2aa)), idx2aa)
    ax.set_xticks(np.arange(len(idx2aa)), np.arange(len(idx2aa)))
    savefig(ax, "part2_eigenvectors", extension, do_save)

if __name__=="__main__":
    import os
    import argparse

    parser = argparse.ArgumentParser(description="Runs the transition table experiments")
    parser.add_argument("--n-trials", type=int, default=100,
                        help="The number of repeated trials to run")
    parser.add_argument("--n-steps", type=int, default=10000,
                        help="How many time steps to simulate")
    parser.add_argument("--proportion", type=float, default=1e-3,
                        help="The proportion of amino acids which receive the mutation")
    parser.add_argument("--save", type=bool, default=False,
                        help="Whether to save the figures or display them")
    parser.add_argument("--extension", type=str, default="pdf",
                        help="Which extension to use to save the figures")
    parser.add_argument("--seed", type=int, default=None,
                        help="The seed used for the pseudo-random number generator")
    parser.add_argument("--theme", type=str, default="light",
                        help="Whether to use a \"dark\" or \"light\" theme for the plots")
    args = parser.parse_args()

    # A parameter dictionary for these experiments
    params = {}

    # Set Matplotlib parameters
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.size"] = 24
    if args.theme == "dark":
        plt.style.use("dark_background")
    factor = 8
    params["figsize"] = (factor * 1.618, factor)
    params["do_save"] = args.save
    params["extension"] = args.extension

    # Load command-line args into parameter dictionary
    params["n_trials"] = args.n_trials
    params["n_timesteps"] = args.n_steps
    params["prop"] = args.proportion

    # Seed the PRNG
    if args.seed is None:
        np.random.seed(int(time.time() * 1000) % 10000)
    else:
        np.random.seed(args.seed)

    # Create and plot the transition table
    trans_table = make_trans_prop(params["prop"])
    plot_trans_table(trans_table, **params)

    eigenvecs, pcs = eigen(trans_table, **params)
    plot_eigenvectors(eigenvecs, **params)

    for i, row in enumerate(eigenvecs.T):
        print(f"        {idx2aa[i]}:", end="")
        for val in row:
            print(f"\t{val:7.4f}", end="")
        print()

#     exp_from_state(state=None, eigenvecs, n_timesteps, fig, ax)
    dist, _, _ = experiment(params["prop"], params["n_timesteps"])
    fig, ax = plt.subplots(figsize=params["figsize"])
    img = ax.imshow(dist.T, aspect="auto", interpolation="none", cmap="plasma")
    fig.colorbar(img, ax=ax)
    ax.set_yticks(np.arange(len(idx2aa)), idx2aa)
    ax.set_title("Change in distribution of amino acids")
    ax.set_xlabel("Iterations")
    ax.set_ylabel("Amino acid")
    savefig(ax, "part2_distribution", params["extension"], params["do_save"])

    n, d = eigenvecs.shape
    eigendists = np.zeros((n, params["n_timesteps"]))
    for i in range(n):
        for j in range(params["n_timesteps"]):
            eigendists[i, j] = np.dot(eigenvecs[i], dist[j] / l2(dist[j]))

    fig, ax = plt.subplots(figsize=params["figsize"])
#     ax.imshow(eigendists, aspect="auto", cmap="plasma", interpolation="none")
#     fig.colorbar(img, ax=ax)
#     ax.set_yticks(np.arange(len(idx2aa)), np.arange(len(idx2aa)))
    for i, eigendist in enumerate(eigendists):
        ax.plot(eigendist, color="tab:blue")
        ax.text(params["n_timesteps"], eigendist[-1], f"eigenvector {i:>2d}", horizontalalignment="right", verticalalignment="bottom", fontsize="xx-small")

    ax.set_title("Change in distance from eigenvectors")
    ax.set_xlabel("Iterations")
    ax.set_ylabel("Distance from eigenvectors")
#     ax.set_ylim([0.0, 1.3])
    savefig(ax, "part2_distance_eig", params["extension"], params["do_save"])

    # Probability progression
    nonsyns = np.zeros((params["n_trials"], params["n_timesteps"]), dtype=dtype)
    dists = np.zeros((params["n_trials"], params["n_timesteps"]), dtype=dtype)
    _, n = dist.shape
    even = np.ones(n, dtype=dtype) / n
    for trial in range(params["n_trials"]):
        dist, nonsyn, _ = experiment(params["prop"], params["n_timesteps"])

        for j, tmp in enumerate(dist):
            dists[trial, j] = distance(tmp, even)
        nonsyns[trial] = nonsyn

    mins = np.min(nonsyns, axis=0)
    maxs = np.max(nonsyns, axis=0)
    means = np.mean(nonsyns, axis=0)

    fig, ax = plt.subplots(figsize=params["figsize"])
    ax.fill_between(np.arange(params["n_timesteps"]), mins, maxs, alpha=0.2)
    ax.plot(means)
#     ax.hlines(0.5, 0, params["n_timesteps"], color="r", linestyle="--", label="min")
#     ax.hlines(1.0, 0, params["n_timesteps"], color="r", linestyle="--", label="max")
    ax.set_title("Change in probability of nonsynonymous mutation")
    ax.set_xlabel("Iterations")
    ax.set_ylabel("Probability of nonsynonymouse mutation")
    savefig(ax, "part2_nonsyn", params["extension"], params["do_save"])

    mins = np.min(dists, axis=0)
    maxs = np.max(dists, axis=0)
    means = np.mean(dists, axis=0)
    fig, ax = plt.subplots(figsize=params["figsize"])
    ax.fill_between(np.arange(params["n_timesteps"]), mins, maxs, alpha=0.2)
    ax.plot(means)
    ax.set_title("Change in distance from uniform distribution")
    ax.set_xlabel("Iterations")
    ax.set_ylabel("Distance from uniform distribution")
    savefig(ax, "part2_distance", params["extension"], params["do_save"])

    if not params["do_save"]:
        plt.show()
