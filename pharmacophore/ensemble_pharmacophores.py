import tempfile
import os

from pdb_python_api import PDBResult
from hotspots.hs_pharmacophore import PharmacophoreModel
from ccdc.protein import Protein

from hotspots.calculation import Runner
from hotspots.hs_io import HotspotWriter

tmp = tempfile.mkdtemp()
out_dir = "/home/pcurran/patel/CDK2/1aq1_ensemble"
centre_pdb = PDBResult("1aq1")
centre_pdb.download(out_dir)

ensemble_members = [("3QTU","X44"), ("2VTL","LZ5"), ("1OIT","HDT"), ("3R6X","X84")]

pdbs = [PDBResult(x[0]) for x in ensemble_members]
for i, pdb in enumerate(pdbs):
    pdb.download(tmp)
    pdb.clustered_ligand = ensemble_members[i][1]

#prots = [Protein.from_file(os.path.join(tmp, member[0] + ".pdb")) for member in ensemble_members]

proteins, ligands = PharmacophoreModel._align_proteins(centre_pdb, "A", pdbs)

# create ensemble
for p in proteins:
    print p.identifier
    p.remove_all_metals()
    for l in p.ligands:
        p.remove_ligand(l.identifier)
    p.add_hydrogens()

    h = Runner()

    # hotspot calculation settings
    s = h.Settings()
    s.apolar_translation_threshold = 15
    s.polar_translation_threshold = 15
    s.polar_contributions = False
    s.nrotations = 3000
    s.sphere_maps = True

    hr = h.from_protein(protein=p, charged_probes=False, buriedness_method='ghecom', nprocesses=3,
                    settings=s, cavities=None)

    out_settings = HotspotWriter.Settings()
    out_settings.charged = False

    out = os.path.join(out_dir, p.identifier)
    if not os.path.exists(out):
        os.mkdir(out)

    with HotspotWriter(out, grid_extension=".grd", zip_results=True,
                       settings=out_settings) as w:
        w.write(hr)

# read

# Mean map

