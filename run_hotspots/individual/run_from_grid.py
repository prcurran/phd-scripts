from fragment_hotspot_maps.fragment_hotspot_maps import Hotspots
from ccdc.protein import Protein
from ccdc import io
from ccdc.io import MoleculeReader
from ccdc.utilities import Grid
import sys
from ccdc.molecule import Coordinates
from os.path import join, dirname, exists
from os import mkdir
import argparse

from ccdc.io import MoleculeReader


class Runner(argparse.ArgumentParser):
    """ class to initiate a HotspotResults instance"""
    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)
        self.add_argument(
            'out_dir',
            help='Directory containing the Hotspot Grids'
        )
        self.args = self.parse_args()

    def run(self):
        """create HotspotResult, runs boundary calculation, generates pharmacophores"""
        h = Hotspots()
        # read
        result = h._from_grid_dic(super_grids={"apolar": Grid.from_file(join(self.args.out_dir, "apolar.grd")),
                                               "donor": Grid.from_file(join(out_dir, "donor.grd")),
                                               "acceptor": Grid.from_file(join(out_dir, "acceptor.grd")),
                                               "negative": Grid.from_file(join(out_dir, "negative.grd")),
                                               "positive": Grid.from_file(join(out_dir, "positive.grd"))
                                               },
                                  prot=Protein.from_file(join(out_dir,"protein.pdb")),
                                  sampled_probes={"apolar": MoleculeReader(join(out_dir, "apolar_probes.mol2")),
                                                  "donor": MoleculeReader(join(out_dir, "donor_probes.mol2")),
                                                  "acceptor": MoleculeReader(join(out_dir, "acceptor_probes.mol2")),
                                                  "negative": MoleculeReader(join(out_dir, "negative_probes.mol2")),
                                                  "positive": MoleculeReader(join(out_dir, "positive_probes.mol2"))
                                                  }
                                    )
        # boundary
        objs = result.extract_hotspots(out_dir=dirname(out_dir), pharmacophores=True)

        for i, hr in enumerate(objs):
            out = join(dirname(out_dir), "hotspot_boundaries", "{}".format(i))

            if not exists(join(dirname(out_dir), "hotspot_boundaries")):
                mkdir((join(dirname(out_dir), "hotspot_boundaries")))

            if not exists(out):
                mkdir(out)
            hr.out_dir = out
            hr.output_pymol_file()
            hr.pharmacophore.out_dir = out
            hr.pharmacophore.write(join(out,"pharmacophore.json"))

if __name__ == "__main__":
    r = Runner()
    r.run()
