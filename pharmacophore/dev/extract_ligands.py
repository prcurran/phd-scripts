from ccdc.io import MoleculeReader, MoleculeWriter
from ccdc.protein import Protein
from os.path import join, dirname, splitext
from os import listdir
from ccdc.cavity import Cavity

from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.ML.Cluster import Butina
from pprint import pprint



excluded_hetids = ['3CO', '3CN', 'REO', 'BA', 'CUZ', 'SF4', 'EOH', 'ZNO',
                   'NAO', 'EOM', 'CUO', 'CCN', 'NAW', 'BR', 'NI2', 'NI1',
                   'O2', 'RU', 'TAS', 'TL', 'OF2', 'OF3', 'OF1', 'RB', 'BF4',
                   'IOD', 'C2O', 'PB', '6MO', 'GD', 'PER', 'GA', 'EHN', 'U1',
                   'TE', 'C2C', 'OH', 'NA5', 'ZN2', 'NA6', 'OX', 'OS', 'AR',
                   'CD1', 'CMO', 'OCL', 'EMC', 'OCN', 'MO3', 'NGN', 'BEF',
                   'HDZ', 'MO5', 'MO4', 'VO4', 'HO', 'PCL', 'CB5', 'HG',
                   'NO2', 'NO3', 'PR', 'SCN', 'IUM', 'PT', 'FEL', 'YT3',
                   'TCN', 'ZO3', 'PD', 'PI', 'MH3', 'AF3', 'ZN', 'OXY', 'MN3',
                   'MN5', 'MN6', 'S', 'O', 'W', 'NMO', 'EU', 'NO', 'PT4',
                   'MW2', 'MG', '543', 'IRI', 'MN', 'PC4', 'MOH', 'RE', 'OC3',
                   'OC2', 'I42', 'OC4', 'OC7', 'OC6', 'GD3', '3OF', 'MM4',
                   'IN', 'DOD', 'CD5', 'PEO', 'FE', 'PBM', 'CD3', 'ETI', 'NI',
                   'AST', 'NA', 'AU3', 'FE2', 'SEK', 'MO1', 'EU3', 'AUC',
                   'CLO', 'CO', 'DUM', 'CL', 'MO7', 'CA', 'MW1', 'CE', 'XE',
                   'CUA', 'V', 'CS', 'CR', 'CU', 'SR', 'H2S', 'ND4', 'BRO',
                   'KR', 'NME', 'SM', 'SB', 'SE', 'AZI', 'FLO', 'YB', 'MSM',
                   'SE4', 'BF2', '2NO', '2MO', 'DIS', 'LA', 'LI', '2BM',
                   'MOO', 'PS5', 'TB', 'CYN', 'MBR', '202', 'AG', 'IR', 'HOH',
                   '3NI', 'NH4', 'SO2', 'AU', 'MTO', 'NH3', 'HGI', 'CNN',
                   'HGC', 'ARS', 'ART', 'TFH', '6WO', 'E1H', 'ARF', 'Y1',
                   'CFT', 'IDO', 'KO4', 'NRU', '4MO', 'ACT', 'CME', 'HEM', 'GOL',
                   'PO4'

                   ]


def ClusterFps(fps, cutoff=0.2):
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs


def extracted_ligands(fnames, pdb_codes, cavities):
    """get ligands"""
    ligands_by_cavtiy = {}
    for ident in range(0, len(cavities)):
        ligands_by_cavtiy.update({ident: []})

    for i, fname in enumerate(fnames):
        print i, "/", len(fnames)
        prot = Protein.from_file(fname)
        prot.remove_all_waters()
        prot.remove_all_metals()
        prot.detect_ligand_bonds()
        for ligand in prot.ligands:
            centroid = ligand.centre_of_geometry()
            for j, cavity in enumerate(cavities):
                if contains(cavity, centroid):
                    try:
                        het = ligand.identifier.split(":")[1][0:3]
                    except IndexError:
                        het = ligand.identifier
                    if len(ligand.atoms) > 5 and het not in excluded_hetids:
                        print het
                        ligand.add_hydrogens()
                        try:
                            ligand.identifier = "{0} ({1})".format(pdb_codes[i].upper(), het)
                            print ligand.identifier
                        except IndexError:
                            print ligand.identifier
                        ligands_by_cavtiy[j].append(ligand)
    return ligands_by_cavtiy


def get_fnames(base):
    """get fnames"""
    return [join(base, f) for f in listdir(base)], [f.split(".")[0] for f in listdir(base)]


def contains(cav, point, tolerance=0):
    mini = cav.bounding_box[0]
    maxi = cav.bounding_box[1]
    return all([mini.x - tolerance < point[0] < maxi.x + tolerance,
                mini.y - tolerance < point[1] < maxi.y + tolerance,
                mini.z - tolerance < point[2] < maxi.z + tolerance])


def deduplicate(ligand_by_cavity):
    new_dict = {}
    for i, key in enumerate(ligand_by_cavity.keys()):
        seen = []
        new_dict.update({i:[]})
        for mol in ligand_by_cavity[key]:
            if mol.identifier.split(" ")[1] not in seen:
                seen.append(mol.identifier.split(" ")[1])
                new_dict[i].append(mol)
            else:
                print "seen"
    return new_dict
######################################################################################################################
target = "CDK2"
pdb = "1aq1"
in_dir = "Z:/fragment-hotspot-results/patel_set/{}/{}/gold_standard".format(target, pdb)
prot = join(dirname(in_dir), "{}.pdb".format(pdb))
######################################################################################################################


fnames, pdb_codes = get_fnames(join(in_dir, "output"))
cavities = Cavity.from_pdb_file(prot)
ligand_by_cavity_mix = extracted_ligands(fnames, pdb_codes, cavities)
#ligand_by_cavity = deduplicate(ligand_by_cavity_mix)
ligand_by_cavity = ligand_by_cavity_mix

for c, mols in ligand_by_cavity.items():
    with MoleculeWriter(join(in_dir, "cavity_mol_{}.sdf".format(c))) as w:
        for m in mols:
            w.write(m)

for k in range(len(cavities)):
    if len(ligand_by_cavity[k]) > 0:
        in_mols = join(in_dir, "cavity_mol_{}.sdf".format(k))
        ms = [x for x in Chem.ForwardSDMolSupplier(in_mols) if x is not None]

        #print len(ms)
        #fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in ms]
        fps = [MACCSkeys.GenMACCSKeys(x) for x in ms]

        print len(fps), len(ms)
        clusters = ClusterFps(fps, cutoff=0.3)

        member = [x[0] for x in clusters if len(x) > 1]
        # m = [x[1] for x in clusters if len(x) > 1]
        # member.extend(m)

        print clusters

        mols = [ms[int(m)-1] for m in member]
        print "in", len(ms)
        print "out", len(mols)

        writer = Chem.SDWriter(join(in_dir, "cluster_mol_{}.sdf".format(k)))
        for j, m in enumerate(mols):
            writer.write(m)
