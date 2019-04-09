from __future__ import print_function, division
import tempfile
import luigi
import os
import pickle
from collections import OrderedDict
import numpy as np
from numpy import ma
import seaborn as sns
import pandas as pd
from scipy.spatial import distance

from ccdc.io import MoleculeWriter
from ccdc.protein import Protein

from hotspots.pdb_python_api import PDBResult
from hotspots.hs_io import HotspotWriter, HotspotReader
from hotspots.calculation import Runner
from hotspots.result import Extractor, Results
from hotspots.hs_pharmacophore import PharmacophoreModel


class GetHotspots(luigi.Task):
    """
    Creates "main" Result
    """
    base = luigi.Parameter()
    target = luigi.Parameter()
    pdb = luigi.Parameter()
    chain = luigi.Parameter()

    def requires(self): pass

    def output(self):
        return luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "out.zip"))

    def run(self):
        h = Runner()

        # hotspot calculation settings
        s = h.Settings()
        s.apolar_translation_threshold = 15
        s.polar_translation_threshold = 15
        s.polar_contributions = False
        s.nrotations = 3000

        hr = h.from_pdb(pdb_code=self.pdb, charged_probes=False, buriedness_method='ghecom', nprocesses=1,
                        settings=s, cavities=None)

        out_settings = HotspotWriter.Settings()
        out_settings.charged = False

        with HotspotWriter(os.path.dirname(self.output().path), grid_extension=".grd",zip_results=True,
                           settings=out_settings) as w:
            w.write(hr)


class GetBestVolume(luigi.Task):
    """
    Creates "best volume" Result
    """
    base = luigi.Parameter()
    target = luigi.Parameter()
    pdb = luigi.Parameter()
    chain = luigi.Parameter()

    def requires(self): return GetHotspots(self.base, self.target, self.pdb, self.chain)

    def output(self): return luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "bestvol", "out.zip"))

    def run(self):
        hs = HotspotReader(self.input().path).read()

        settings = Extractor.Settings()
        settings.pharmacophore = False
        settings.cutoff = 12
        settings.mvon = True

        extractor = Extractor(hs)
        best = extractor.extract_best_volume(volume=250)[0]

        out_settings = HotspotWriter.Settings()
        out_settings.charged = False

        with HotspotWriter(os.path.dirname(self.output().path), grid_extension=".grd", zip_results=True,
                           settings=out_settings) as w:
            w.write(best)


class GetHotspotPharmacophore(luigi.Task):
    """
    Creates "Hotspot Pharmacophore"
    """
    base = luigi.Parameter()
    target = luigi.Parameter()
    pdb = luigi.Parameter()
    chain = luigi.Parameter()

    def requires(self):
        return GetBestVolume(self.base, self.target, self.pdb, self.chain)

    def output(self):
        return {"pymol": luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "hs_pharmacophore", "hotspot_pharmacophore.py")),
                "points": luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "hs_pharmacophore", "hotspot_points.pkl"))}

    def run(self):
        bestvol = HotspotReader(self.input().path).read()
        pharmacophore = bestvol.get_pharmacophore_model()

        pharmacophore.write(self.output()['pymol'].path)

        points = pharmacophore._comparision_dict()
        with open(self.output()['points'].path, 'wb') as w:
            pickle.dump(points, w)


class GetRepresentitives(luigi.Task):
    """
    creates (best volume) Results
    """
    base = luigi.Parameter()
    target = luigi.Parameter()
    pdb = luigi.Parameter()
    chain = luigi.Parameter()

    def requires(self): pass

    def output(self): return luigi.LocalTarget(os.path.join(self.base, self.target, "representatives.dat"))

    def run(self):
        # generate representatives
        ref = PharmacophoreModel.from_pdb(pdb_code=self.pdb, chain=self.chain, identifier=self.pdb)
        if len(ref.representatives) > 2:
            ligands = ref.representatives
        else:
            ligands = ref.all_ligands

        with open(self.output().path, "w") as f:
            for l in ligands:
                f.write("{},{},{}\n".format(l.structure_id, l.chemical_id, l.smiles))


class GetReferencePharmacophore(luigi.Task):
    """
    creates (best volume) Results
    """
    base = luigi.Parameter()
    target = luigi.Parameter()
    pdb = luigi.Parameter()
    chain = luigi.Parameter()

    def requires(self): return GetRepresentitives(self.base, self.target, self.pdb, self.chain)

    def output(self):
        return {"pymol": luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "reference", "reference_pharmacophore.py")),
                "grids": luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "reference", "out.zip")),
                "aligned_mols": luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "reference", "aligned.mol2")),
                "points": luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "reference", "reference_points.pkl")),
                }

    def run(self):
        # create pharmacophore
        ref = PharmacophoreModel.from_pdb(pdb_code=self.pdb, chain=self.chain,
                                          representatives=self.input().path, identifier=self.pdb)
        ref.rank_features(max_features=6, feature_threshold=5)

        # write pymol file
        ref.write(self.output()["pymol"].path)

        # write Results file
        temp = tempfile.mkdtemp()
        PDBResult(self.pdb).download(temp)
        result = Results(protein=Protein.from_file(os.path.join(temp, "{}.pdb".format(self.pdb))),
                         super_grids=ref.dic)

        out_settings = HotspotWriter.Settings()
        out_settings.charged = False
        with HotspotWriter(os.path.dirname(self.output()["grids"].path), grid_extension=".grd", zip_results=True,
                           settings=out_settings) as w:
            w.write(result)

        # write aligned molecules
        with MoleculeWriter(self.output()['aligned_mols'].path) as w:
            for l in ref.aligned_ligands:
                w.write(l)

        # points
        points = ref._comparision_dict()
        with open(self.output()['points'].path, 'wb') as w:
            pickle.dump(points, w)


class Compare(luigi.Task):
    """
    - Pharmacophoric points labelled with "type_score_rank"
    - Order Hotspot and Referenec Pharmacophore by rank
    - Calculate all by all distance
    - < 2A = pass
    """
    base = luigi.Parameter()
    target = luigi.Parameter()
    pdb = luigi.Parameter()
    chain = luigi.Parameter()

    def requires(self): return {"hotspot": GetHotspotPharmacophore(self.base, self.target, self.pdb, self.chain),
                                "ref": GetReferencePharmacophore(self.base, self.target, self.pdb, self.chain)}

    def output(self): return {"df": luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "dataframe.csv")),
                              "apolar": luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "reference", "apolar.csv")),
                              "donor": luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "reference", "donor.csv")),
                              "acceptor": luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "reference", "acceptor.csv")),
                              }

    def run(self):
        def get_labels(pharmacophore_dic):
            """
            generate label dictionary
            :param pharmacophore_dic:
            :return:
            """
            dic = {}
            labels = ["{1:.2f}\n{0}".format(v[0], k) for k, v in pharmacophore_dic.iteritems()]
            r = range(0, len(labels))
            for i in r:
                dic.update({i: labels[i]})
            return dic

        # unpickle and order
        with open(self.input()['hotspot']['points'].path, 'rb') as f:
            hotspot = pickle.load(f)

        with open(self.input()['ref']['points'].path, 'rb') as f:
            ref = pickle.load(f)

        # generate ordered data
        a = OrderedDict(sorted(hotspot.items(), reverse=True))
        b = OrderedDict(sorted(ref.items(), reverse=True))

        # create list of probe types for mask generation
        a_probes = [i[0] for i in a.values()]
        b_probes = [i[0] for i in b.values()]

        # create mask on iteraction type: show = 0, hide = 1
        blank = np.ones((len(a_probes), len(b_probes)))
        masks = {"donor": blank.copy(), "acceptor": blank.copy(), "apolar": blank.copy()}
        for ident in masks.keys():
            for i, a_val in enumerate(a_probes):
                for j, b_val in enumerate(b_probes):
                    if a_val == ident and b_val == ident:
                        masks[ident][i, j] = 0

        # get tick labels
        a_dic = get_labels(a)
        b_dic = get_labels(b)

        # calculate all by all distances
        all_distances = distance.cdist(np.array([i[1] for i in a.values()]),
                                       np.array([j[1] for j in b.values()]),
                                       'euclidean')

        # mask obvious fails
        m2 = ma.masked_where(all_distances > 6, all_distances).mask.astype(int)
        for k, m1 in masks.items():
            masks[k] = ma.mask_or(m1.astype(int), m2)

        for probe, array in masks.items():
            mdf = pd.DataFrame(array.astype(int))
            mdf.to_csv(self.output()[probe].path)

        # create dataframe
        df = pd.DataFrame(all_distances)
        df = df.rename(a_dic, axis='rows')
        df = df.rename(b_dic, axis='columns')

        df.to_csv(self.output()['df'].path)


class Plot(luigi.Task):
    """
    Plot heatmap
    """
    base = luigi.Parameter()
    target = luigi.Parameter()
    pdb = luigi.Parameter()
    chain = luigi.Parameter()

    def requires(self): return Compare(self.base, self.target, self.pdb, self.chain)

    def output(self): return luigi.LocalTarget(os.path.join(self.base, self.target, self.pdb, "comparision.png"))

    def run(self):

        masks = {"apolar": 0, "acceptor": 0, "donor": 0}
        for probe in masks.keys():
            masks[probe] = pd.read_csv(self.input()[probe].path, index_col=0).values

        df = pd.read_csv(self.input()['df'].path, index_col=0)
        sns.set_style("white")

        print(df.values.shape)
        print(masks["apolar"].shape)

        ax = sns.heatmap(df, annot=True, mask=masks["donor"], cmap="Blues_r", annot_kws={"size": 8}, fmt='.1f',
                         cbar=False, vmin=0, vmax=3)
        ax = sns.heatmap(df, annot=True, mask=masks["acceptor"], cmap="Reds_r", annot_kws={"size": 8}, fmt='.1f',
                         cbar=False, vmin=0, vmax=3)
        ax = sns.heatmap(df, annot=True, mask=masks["apolar"], cmap="Greys_r", annot_kws={"size": 8}, fmt='.1f',
                         cbar=False, vmin=0, vmax=3)

        ax.set_xlabel("Reference Pharmacophore", fontsize=8, fontweight='bold')
        ax.set_ylabel("Hotspot Pharmacophore", fontsize=8, fontweight='bold')
        ax.tick_params(labelsize=8)

        # output
        ax.figure.savefig(self.output().path)


class LotsOTasks(luigi.WrapperTask):
    """
    The organising class,
    """

    def requires(self):

        base = "/home/pcurran/patel"

        targets = {"CDK2": ["1hcl", "1aq1"],
                   "DHFR": ["1drf", "2w9t"],
                   "Thrombin": ["1c4v", "1vr1"],
                   "HIVRT": ["1tvr", "1dlo"],
                   "A2Ar": ["2ydo"]}

        chains = {"1hcl": "A",
                  "1aq1": "A",
                  "1drf": "A",
                  "2w9t": "A",
                  "1c4v": "2",
                  "1vr1": "H",
                  "1tvr": "A",
                  "1dlo": "A",
                  "2ydo": "A",
                  }

        for target, pdbs in targets.items():
            for pdb in pdbs:
                chain = chains[pdb]
                yield Plot(base=base, target=target, pdb=pdb, chain=chain)


if __name__ == '__main__':
    luigi.build([LotsOTasks()], workers=7)
    #luigi.build([LotsOTasks(), AggregateData()], workers=7)