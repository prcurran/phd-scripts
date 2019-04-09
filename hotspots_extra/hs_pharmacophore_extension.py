# from __future__ import print_function, division
#
# import os
# import tempfile
#
# import numpy as np
# from ccdc.protein import Protein
# from pdb_python_api import Query, PDB, PDBResult
# from rdkit import Chem, DataStructs
# # from rdkit.Chem import MACCSkey, AllChem, Draw
# from rdkit.Chem import MACCSkeys, AllChem
#
# from sklearn.cluster import KMeans
# from sklearn.metrics import silhouette_score
# from ccdc import io
# from tqdm import tqdm
#
# from hotspots import hs_pharmacophore
#
#
# class PharmacophoreModel(hs_pharmacophore.PharmacophoreModel):
#     """
#     a class to handle the the generation of diverse set of ligands
#     """
#
#     # @staticmethod
#     # def _img(ligands, out_dir, fname="ligands.png"):
#     #     out = []
#     #     for l in ligands:
#     #         lig = l.rdmol
#     #         lig.SetProp("_Name", str(l.structure_id + "/" + l.chemical_id))
#     #         out.append(lig)
#     #
#     #     img = Draw.MolsToGridImage(mols=out,
#     #                                molsPerRow=10,
#     #                                subImgSize=(200, 200),
#     #                                legends=[o.GetProp("_Name") for o in out])
#     #
#     #     img.save(os.path.join(out_dir, fname))
#
#     @staticmethod
#     def _to_file(ligands, out_dir, fname="ligands.dat"):
#         with open(os.path.join(out_dir, fname), "w") as f:
#             for l in ligands:
#                 f.write("{},{},{}\n".format(l.structure_id, l.chemical_id, l.smiles))
#
#     @staticmethod
#     def from_pdb(pdb_code, chain, out_dir=None, representatives=None, identifier="gold_standard"):
#         """
#
#         :return:
#         """
#         temp = tempfile.mkdtemp()
#         ref = PDBResult(pdb_code)
#         ref.download(out_dir=temp, compressed=False)
#
#         if representatives:
#             print("Reading representative PDB codes ...")
#             reps = []
#             f = open(representatives, "r")
#             entries = f.read().splitlines()
#
#             for entry in entries:
#                 pdb_code, hetid, smiles = entry.split(",")
#                 reps.append((pdb_code, hetid))
#
#         else:
#             accession_id = PDBResult(pdb_code).protein.sub_units[0].accession_id
#             results = PharmacophoreModel._run_query(accession_id)
#             ligands = PharmacophoreModel._get_ligands(results)
#
#             top = os.path.dirname(os.path.dirname(out_dir))
#             PharmacophoreModel._to_file(ligands, top, fname="all_ligands.dat")
#
#             k = int(round(len(ligands) / 5))
#             if k < 2:
#                 k = 2
#             cluster_dict, s = PharmacophoreModel._cluster_ligands(n=k, ligands=ligands)
#             reps = [l[0] for l in cluster_dict.values() if len(l) != 0]
#
#         targets = []
#
#         for rep in reps:
#             try:
#                 r = PDBResult(identifier=rep.structure_id)
#                 r.clustered_ligand = rep.chemical_id
#
#             except:
#                 try:
#                     r = PDBResult(identifier=rep[0])
#                     r.clustered_ligand = rep[1]
#                 except:
#                     raise AttributeError
#
#             r.download(out_dir=temp, compressed=False)
#             targets.append(r)
#
#         prots, ligands = PharmacophoreModel._align_proteins(reference=ref,
#                                                             reference_chain=chain,
#                                                             targets=targets)
#
#         try:
#             with open(os.path.join(top, "silhouette.dat"), "w") as sfile:
#                 sfile.write("{}".format(s))
#             PharmacophoreModel._to_file(reps, top, fname="representatives.dat")
#
#         except:
#             pass
#
#         if out_dir:
#             with io.MoleculeWriter(os.path.join(out_dir, "aligned_mols.mol2")) as w:
#                 for l in ligands:
#                     w.write(l)
#
#         return PharmacophoreModel.from_ligands(ligands=ligands,
#                                                identifier=identifier)
#
#     @staticmethod
#     def _run_query(accession_id):
#         """
#
#         :return:
#         """
#         # create query
#         q = Query()
#         q.add_term(query_type="UpAccessionIdQuery",
#                    query_parameters={'accessionIdList': '{}'.format(accession_id)})
#
#         q.add_term(query_type="NoLigandQuery",
#                    query_parameters={'haveLigands': 'yes'},
#                    conjunction='and')
#
#         q.add_term(query_type="ResolutionQuery",
#                    query_parameters={'refine.ls_d_res_high.comparator': 'between',
#                                      'refine.ls_d_res_high.min': '0',
#                                      'refine.ls_d_res_high.max': '2.5'},
#                    conjunction='and')
#
#         return PDB.search(q.query)
#
#     @staticmethod
#     def _get_ligands(results):
#         """
#
#         :return:
#         """
#         ligs = []
#         uniques = []
#         for entry in results:
#             for l in entry.filtered_ligands:
#                 try:
#                     l.rdmol = Chem.MolFromSmiles(l.smiles)
#                     l.rdmol.SetProp("_Name", str(entry.identifier + "/" + l.chemical_id))
#                     l.fingerprint = MACCSkeys.GenMACCSKeys(l.rdmol)
#                     if l.chemical_id not in uniques:
#                         ligs.append(l)
#                         uniques.append(l.chemical_id)
#                 except AttributeError:
#                     continue
#         return ligs
#
#     @staticmethod
#     def _cluster_ligands(n, ligands):
#         """
#         generate an all by all similarity matrix
#         :return:
#         """
#         cluster_dict = {a: [] for a in range(int(n))}
#         num = len(ligands)
#         sim = np.zeros((num, num))
#         for i in range(num):
#             for j in range(num):
#                 sim[i, j] = DataStructs.FingerprintSimilarity(ligands[i].fingerprint,
#                                                               ligands[j].fingerprint)
#
#         kmeans_model = KMeans(n_clusters=n, random_state=1).fit(sim)
#         labels = kmeans_model.labels_
#         s = silhouette_score(sim, labels, metric='euclidean')
#
#         for i, ident in enumerate(kmeans_model.labels_):
#             ligands[i].cluster_id = ident
#             cluster_dict[int(ident)].append(ligands[i])
#
#         return cluster_dict, s
#
#     @staticmethod
#     def _align_proteins(reference, reference_chain, targets):
#         """
#
#         :return:
#         """
#         print("Aligning proteins to {}, chain {}...".format(reference.identifier, reference_chain))
#         aligned_prots = []
#         aligned_ligands = []
#
#         reference = Protein.from_file(reference.fname)
#         reference.add_hydrogens()
#
#         for t in tqdm(targets):
#             prot = Protein.from_file(t.fname)
#             prot.detect_ligand_bonds()
#             prot.add_hydrogens()
#             for l in prot.ligands:
#                 if str(t.clustered_ligand) == str(l.identifier.split(":")[1][0:3]):
#                     try:
#                         bs = Protein.BindingSiteFromMolecule(protein=prot,
#                                                              molecule=l,
#                                                              distance=6)
#                         chain = bs.residues[0].identifier.split(":")[0]
#                     except:
#                         break
#                     break
#
#                 else:
#                     continue
#             if not chain:
#                 print("\n        {} failed! No chain detected".format(t.identifier))
#                 break
#             try:
#                 binding_site_superposition = Protein.ChainSuperposition()
#                 (bs_rmsd, bs_transformation) = binding_site_superposition.superpose(reference[reference_chain],
#                                                                                     prot[chain])
#                 aligned_prots.append(prot)
#                 for lig in prot.ligands:
#                     if str(t.clustered_ligand) == str(lig.identifier.split(":")[1][0:3]):
#                         if chain == str(lig.identifier.split(":")[0]):
#                             aligned_ligands.append(lig)
#             except IndexError:
#                 print("\n        {} failed!".format(t.identifier))
#                 continue
#
#         return aligned_prots, aligned_ligands
