from ccdc.io import MoleculeReader
from pprint import pprint

mols = MoleculeReader("Z:/fragment-hotspot-results/patel_set/2/gold_standard/cluster_mol_1.sdf")

ids = [str(mol.identifier) for mol in mols]

pprint(ids)