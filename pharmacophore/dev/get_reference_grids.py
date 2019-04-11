import tempfile
from os.path import join

from ccdc.io import MoleculeReader, MoleculeWriter
from ccdc.protein import Protein
from hotspots.calculation import Results
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.hs_pharmacophore import PharmacophoreModel
from hotspots.hs_utilities import Helper


def to_grid(target, pdb):
    out_dir = "Z:/patel_set/{}/{}".format(target, pdb)
    mols = MoleculeReader(join(out_dir, "reference_pharmacophore", "aligned_mols.mol2"))
    p = PharmacophoreModel.from_ligands(ligands=mols,
                                        identifier="test")
    result = Results(super_grids=p.dic,
                     protein=Protein.from_file(join(out_dir, "hs", "{}.pdb".format(pdb)))
                     )

    out = Helper.get_out_dir(join(out_dir, "reference_pharmacophore", "grids"))

    settings = HotspotWriter.Settings()
    settings.isosurface_threshold = [2, 5, 10]

    with HotspotWriter(path=out, zip_results=True, settings=settings) as w:
        w.write(result)


def main():
    targets = {"CDK2": ["1hcl", "1aq1"],
               "DHFR": ["1drf", "2w9t"],
               "Thrombin": ["1c4v", "1vr1"],
               "HIVRT": ["1tvr", "1dlo"],
               "A2Ar": ["2ydo"]}

    for t, values in targets.items():
        for v in values:
            print t, v
            to_grid(t, v)


if __name__ == "__main__()":
    main()
