from fragment_hotspot_maps.fragment_hotspot_maps import Hotspots
from ccdc.protein import Protein
from ccdc import io
from ccdc.io import MoleculeReader
from ccdc.utilities import Grid
import sys
from ccdc.molecule import Coordinates
from os.path import join, dirname, exists, abspath
from os import mkdir
import argparse

from ccdc.io import MoleculeReader, MoleculeWriter


class Runner(argparse.ArgumentParser):
    """ class to initiate a HotspotResults instance"""
    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)
        self.add_argument(
            'prot_fname',
            help='pdb file path'
        )
        self.add_argument(
            '-g', '--ghecom_executable',
            default=None,
            help='path to ghecom executable, if None, Ligsite will be used'
        )
        self.args = self.parse_args()

        self.in_dir = dirname(abspath(self.args.prot_fname))
        self.out_dir = join(self.in_dir, "out")

    def prepare_protein(self):
        """default protein preparation settings on the protein"""
        self.prot = Protein.from_file(self.args.prot_fname)
        self.prot.remove_all_waters()
        for lig in self.prot.ligands:
            self.prot.remove_ligand(lig.identifier)
        self.prot.remove_all_metals()

    def run(self, prepare=False):
        """from fragment hotspot calc from protein"""
        h = Hotspots()
        if prepare:
            self.prepare_protein()
        else:
            self.prot = Protein.from_file(self.args.prot_fname)

        h.out_dir = self.out_dir
        result = h.from_protein(prot=self.prot,
                                charged_probes=True,
                                ghecom_executable=self.args.ghecom_executable,
                                output_directory=self.out_dir)

        result.out_dir = self.out_dir
        result.output_pymol_file()
        for p, mols in result.sampled_probes.items():
            if not exists(self.out_dir):
                mkdir(self.out_dir)
            with MoleculeWriter(join(self.out_dir, "{}_probes.mol2".format(p))) as w:
                for mol in mols:
                    w.write(mol)

        objs = result.extract_hotspots(out_dir=dirname(self.out_dir), pharmacophores=True)
        for i, hr in enumerate(objs):
            out = join(dirname(self.out_dir), "hotspot_boundaries", "{}".format(i))

            if not exists(join(dirname(self.out_dir), "hotspot_boundaries")):
                mkdir((join(dirname(self.out_dir), "hotspot_boundaries")))

            if not exists(out):
                mkdir(out)
            hr.out_dir = out
            hr.output_pymol_file()
            hr.pharmacophore.out_dir = out
            hr.pharmacophore.write(join(out, "pharmacophore.json"))

def main():
    r = Runner()
    r.run()

if __name__ == "__main__":
    main()

