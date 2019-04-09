import os

import matplotlib.pyplot as plt
import numba
import numpy as np
import seaborn as sns
#from sklearn.manifold import TSNE
import umap
import hdbscan
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs


def chemical_space(fname):
    """
    from text file with smiles data, create a chemical space representation
    :param fname:
    :return:
    """
    ligands = []
    X = []

    with open(fname, "r") as f:
        entries = f.read().splitlines()

        for e in entries:
            smiles = e.split(",")[2]
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name", str(e.split(",")[0] + "/" + e.split(",")[1]))
            ligands.append(mol)

        for l in ligands:
            AllChem.Compute2DCoords(l)
            arr = np.zeros((0,))
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
            DataStructs.ConvertToNumpyArray(fp, arr)
            X.append(arr)

    #return TSNE(n_components=3, metric=tanimoto_dist).fit_transform(X)
    return umap.UMAP(n_neighbors=5, min_dist=0.2, metric=tanimoto_dist).fit_transform(X)


@numba.njit()
def tanimoto_dist(a, b):
    """
    calculate the tanimoto distance between two fingerprint arrays
    :param a:
    :param b:
    :return:
    """
    dotprod = np.dot(a, b)
    tc = dotprod / (np.sum(a) + np.sum(b) - dotprod)
    return 1.0 - tc


def main():
    sns.set_context('poster')
    sns.set_style('white')
    sns.set_color_codes()
    plot_kwds = {'alpha': 0.5, 's': 80, 'linewidth': 0}

    targets = ["CDK2", "DHFR", "Thrombin", "HIVRT", "A2Ar"]

    out = "Z:/patel_set"

    for t in targets:
        umap_X = chemical_space(fname=os.path.join(out, t, "all_ligands.dat"))
        # plt.scatter(all_X.T[0], all_X.T[1], color='g', **plot_kwds)
        # plt.show()

        cluster_umap = hdbscan.HDBSCAN(min_cluster_size=2, gen_min_span_tree=True)
        cluster_umap.fit(umap_X)
        print("target: {}\n".format(t), len(umap_X), "to", len(set(cluster_umap.labels_)))
        # ax = cluster_umap.minimum_spanning_tree_.plot(edge_cmap='viridis',
        #                                               edge_alpha=0.6,
        #                                               node_size=20,
        #                                               edge_linewidth=1)
        #
        # # ax = cluster_tsne.single_linkage_tree_.plot(cmap='viridis',
        # #                                             colorbar=True,)
        # ax.figure.show()

        x = [umap_X.T[0][i] for i, l in enumerate(cluster_umap.labels_) if l != -1]
        y = [umap_X.T[1][i] for i, l in enumerate(cluster_umap.labels_) if l != -1]
        hue = [l for i, l in enumerate(cluster_umap.labels_) if l != -1]

        plt.scatter(x, y, c=hue, cmap='plasma')
        plt.title("{} clusters".format(t))
        plt.show()
        plt.close()

if __name__ == "__main__":
    main()
