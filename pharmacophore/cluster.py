from rdkit import Chem
from rdkit.Chem import AllChem
import gzip
from ccdc.protein import Protein
from os import listdir
from os.path import join, dirname, splitext
from ccdc.io import MoleculeWriter, MoleculeReader

def ClusterFps(fps,cutoff=0.2):
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs

def GetMolecules(base):
    mols = []
    files = [join(base, f) for f in listdir(base)]

    with MoleculeWriter(join(dirname(base), "out.sdf")) as w:
        for i, file in enumerate(files):
            print i
            prot = Protein.from_file(file)
            ligands = [l for l in prot.ligands]
            mols.extend(ligands)
            for ligand in ligands:
                w.write(ligand)
    return mols

base = "C:/Users/pcurran/Desktop/paper1/CDK2/output"
#mols = GetMolecules(base)

in_mols = join(dirname(base), "out.sdf")

ms = [x for x in Chem.ForwardSDMolSupplier(in_mols) if x is not None]
print len(ms)
fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in ms]
clusters = ClusterFps(fps,cutoff=0.7)
member = [x[0] for x in clusters]
print len(member)
mols = [ms[int(m)-1] for m in member]

writer = Chem.SDWriter(join(dirname(base), "cluster.sdf"))
for j, m in enumerate(mols):
    #fname = join(dirname(base), "cluster", "cluster_id{}.sdf".format(j))
    #writer = Chem.SDWriter(join(dirname(base), "cluster", "cluster_id{}.sdf".format(j)))
    writer.write(m)


