import os

from hotspots.hs_pharmacophore_extension import PharmacophoreModel
from pdb_python_api import PDBResult

# targets = {"CDK2": ["1hcl", "1aq1"],
#            "DHFR": ["1drf", "2w9t"],
#            "Thrombin": ["1c4v", "1vr1"],
#            "HIVRT": ["1tvr", "1dlo"],
#            "A2Ar": ["2ydo"]}

targets = {"A2Ar": ["2ydo"]}

chain_dic = {"1hcl": "A",
             "1aq1": "A",
             "1drf": "A",
             "2w9t": "A",
             "1c4v": "2",
             "1vr1": "H",
             "1tvr": "A",
             "1dlo": "A",
             "2ydo": "A",
             }

out_dir = "Z:/patel_set"

target = []
pdb = []
all_ligands = []
clustered_ligands = []
silhouette = []

for target, pdbs in targets.items():
    out = os.path.join(out_dir, target)
    reps = os.path.join(out, "representatives.dat")

    if not os.path.exists(out):
        os.mkdir(out)

    for pdb in pdbs:
        print "Running PDB: {} ...".format(pdb)
        out = os.path.join(out_dir, target)
        out = os.path.join(out, pdb)
        if not os.path.exists(out):
            os.mkdir(out)

        out = os.path.join(out, "reference_pharmacophore")
        if not os.path.exists(out):
            os.mkdir(out)

        if os.path.exists(reps):
            representatives = reps
        else:
            representatives = None

        f = PDBResult(identifier=pdb).download(out_dir=out)
        p = PharmacophoreModel.from_pdb(pdb_code=pdb,
                                        chain=chain_dic[pdb],
                                        out_dir=out,
                                        representatives=representatives
                                        )
        p.rank_features(max_features=5)

        fname = "reference_pharmacophore"
        p.write(os.path.join(out, "{}.py".format(fname)))
        p.write(os.path.join(out, "{}.cm".format(fname)))