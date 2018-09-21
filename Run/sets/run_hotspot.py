from fragment_hotspot_maps.fragment_hotspot_maps import Hotspots
from ccdc.protein import Protein
from ccdc import io
from ccdc.io import MoleculeReader
from ccdc.utilities import Grid
import sys

from os.path import join, dirname, exists
from os import mkdir


h = Hotspots()

#batch
s = sys.argv[1]
fname = "/local/pcurran/Patel_set/{}/protein.pdb".format(s)
out_dir = "/local/pcurran/Patel_set/{}/out".format(s)

# prep
prot = Protein.from_file(fname)
prot.remove_all_waters()
for lig in prot.ligands:
    prot.remove_ligand(lig.identifier)
prot.remove_all_metals()

#run
result = h.from_protein(prot=prot,
                        fname=fname,
                        ghecom_executable="/home/pcurran/src",
                        charged_probes=True
                        )

result.out_dir = out_dir

for p, mols in result.sampled_probes.items():
    with io.MoleculeWriter(join(out_dir, "{}_probes.mol2".format(p))) as w:
        for mol in mols:
            w.write(mol)

# extract
objs = result.extract_hotspots(out_dir=dirname(out_dir), pharmacophores=True)

for i, hr in enumerate(objs):
    out = join(dirname(out_dir), "hotspot_boundaries", "{}".format(i))

    if not exists(join(dirname(out_dir), "hotspot_boundaries")):
        mkdir((join(dirname(out_dir), "hotspot_boundaries")))

    if not exists(out):
        mkdir(out)

    hr.output_pymol_file(out)
    hr.pharmacophore.write(join(out, "pharmacophore.json"))

