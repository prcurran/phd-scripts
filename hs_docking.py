import os
import tempfile

import numpy as np
from scipy import stats
from ccdc.cavity import Cavity
from ccdc.docking import Docker
from ccdc.io import MoleculeWriter, MoleculeReader
from ccdc.molecule import Molecule, Atom
from ccdc.protein import Protein
from hotspots import calculation
from hotspots import hs_io
from hotspots import hs_screening
from hotspots import best_volume
from pdb_python_api import PDBResult
from rdkit import Chem
from rdkit.Chem import AllChem


class Organiser(object):

    def __init__(self, path, hs_pdb, ligand_pdb, ligand_identifier):
        self.path = path

        self.hs_pdb = hs_pdb
        self.ligand_pdb = ligand_pdb
        self.ligand_identifier = ligand_identifier
        self.temp = tempfile.mkdtemp()

        # download PDB file using PDB python API
        self.hs_pdb_info = PDBResult(identifier=self.hs_pdb)
        self.hs_pdb_info.download(out_dir=self.temp)
        self._hs_fname = os.path.join(self.temp, hs_pdb + ".pdb")

        # calculate hotspot using Hotspots API
        self.protein = self._prepare_protein(Protein.from_file(self._hs_fname))

        # self.hr = self.calc_hr()
        # with hs_io.HotspotWriter(path=os.path.join(self.path), zip_results=True) as hw:
        #     hw.write(self.hr)
        self.hr = hs_io.HotspotReader(path=os.path.join(self.path, "out.zip")).read()


        # download ligand PDB file using PDB python API
        self.ligand_pdb_info = PDBResult(identifier=self.ligand_pdb)
        self.ligand_pdb_info.download(out_dir=self.temp)
        self._ligand_fname = os.path.join(self.temp, ligand_pdb + ".pdb")

        # align
        target = self._prepare_protein(Protein.from_file(self._ligand_fname), remove_ligands=False)
        self.ligand_protein = self.align(self.protein, "A", target, "A")
        self.ligand = self.extract_ligand()

        # substructure search of the PDB
        self.search_ligands = self.similarity_search()

        # dock search ligands into hotspot protein
        self.docked_ligands = self.dock()
        self.rescored_ligands = self.rescore()

    def rescore(self):
        """
        rescores docked molecules in hotspots
        :return:
        """
        rescored_mols = []
        mol_dic = {}
        self.docked_ligands = [self.hr.score(m) for m in self.docked_ligands]
        for m in self.docked_ligands:
            array = []
            for atm in m.atoms:
                if atm.atomic_number > 0:
                    if type(atm.partial_charge) is not float:
                        array.append(0)
                    else:
                        array.append(atm.partial_charge)
            if len(array) == 0:
                continue

            else:
                g = stats.mstats.gmean(np.array(array))
                try:
                    mol_dic[g].append(m)
                except KeyError:
                    mol_dic.update({g: [m]})

        for key in sorted(mol_dic.keys()):
            rescored_mols.extend(mol_dic[key])

        return rescored_mols

    def dock(self):
        """
        handle docking run with GOLD
        :return:
        """
        docker = Docker()

        # enables hotspot constraints
        docker.settings = hs_screening.DockerSettings()

        f = os.path.join(self.temp, self.hs_pdb + ".mol2")
        with MoleculeWriter(f) as w:
            w.write(self.protein)

        # setup
        docker.settings.add_protein_file(f)
        docker.settings.binding_site = docker.settings.BindingSiteFromPoint(protein=docker.settings.proteins[0],
                                                                            origin=self.ligand.centre_of_geometry(),
                                                                            distance=12.0)

        docker.settings.fitness_function = 'plp'
        docker.settings.autoscale = 10.
        docker.settings.output_directory = self.temp
        docker.settings.output_file = "docked_ligands.mol2"
        docker.settings.add_ligand_file(self.search_ligands, ndocks=3)

        # constraints
        # docker.settings.add_constraint(
        #     docker.settings.TemplateSimilarityConstraint(type="all", template=self.ligand, weight=150)
        #)

        # extractor = best_volume.Extractor(hr=self.hr, volume=300, mode="global", mvon=False)
        # bv = extractor.extracted_hotspots[0]
        #
        # with hs_io.HotspotWriter(path=os.path.join(self.path, "bv")) as hw:
        #     hw.write(extractor.extracted_hotspots)
        #
        # hs = docker.settings.HotspotHBondConstraint.from_hotspot(protein=docker.settings.proteins[0],
        #                                                          hr=bv,
        #                                                          weight=150,
        #                                                          max_constraints=2)
        #
        # docker.settings.add_constraint(hs)
        # docker.settings.add_apolar_fitting_points(hr=self.hr)
        #
        # mol = Molecule(identifier="constraints")
        # for a in hs.atoms:
        #     mol.add_atom(Atom(atomic_symbol="C",
        #                       atomic_number=14,
        #                       label="Du",
        #                       coordinates=a.coordinates))
        #
        # with MoleculeWriter(os.path.join(self.path, "constaints.mol2")) as w:
        #     w.write(mol)

        # dock
        docker.dock()
        return MoleculeReader(os.path.join(docker.settings.output_directory,
                                           docker.settings.output_file))

    def similarity_search(self):
        """
        similarity search on the selected PDB ligand
        :return: str, path to ligand .sdf file
        """
        for lig in self.ligand_pdb_info.ligands:
            if lig.chemical_id == self.ligand_identifier:
                ligands = lig.search(search_type="similarity", tanimoto=0.6)

        mols = [Chem.AddHs(Chem.MolFromSmiles(l.smiles)) for l in ligands]

        f = os.path.join(self.temp, "search_ligands.sdf")
        w = Chem.SDWriter(f)
        for m in mols:
            AllChem.Compute2DCoords(m)
            AllChem.MMFFOptimizeMolecule(m)
            w.write(m)

        w.close()
        return f

    def extract_ligand(self):
        """
        return a `ccdc.molecule.Molecule` object of the desired small moleculae
        :return:
        """
        for l in self.ligand_protein.ligands:
            if l.identifier.split(":")[1][0:3] == self.ligand_identifier:
                ligand = l

        return ligand

    def _prepare_protein(self, prot, remove_ligands=True):
        """
        prepares the protein for hotspot analysis
        :return:
        """
        prot.add_hydrogens()
        prot.detect_ligand_bonds()
        prot.remove_all_metals()
        prot.remove_all_waters()

        if remove_ligands:
            for lig in prot.ligands:
                prot.remove_ligand(lig.identifier)
        return prot

    def align(self, reference, reference_chain, target, target_chain):
        """
        sequence based alignment of target to reference
        :return:
        """
        target.detect_ligand_bonds()
        target.add_hydrogens()
        binding_site_superposition = Protein.ChainSuperposition()
        binding_site_superposition.superpose(reference[reference_chain],
                                             target[target_chain])

        return target

    def calc_hr(self):
        """
        runs the hotspot calculation and returns a `hotspots.calculation.Result`
        :return: `hotspots.calculation.Result`
        """
        h = calculation.Runner()
        settings = calculation.Runner.Settings(sphere_maps=True)
        cavs = Cavity.from_pdb_file(self._hs_fname)
        return h.from_protein(protein=self.protein,
                              charged_probes=True,
                              buriedness_method='ghecom',
                              cavities=cavs,
                              nprocesses=5,
                              settings=settings)


def main():
    r = Organiser(path="/home/pcurran/use_case/no_constraints",
                  hs_pdb="1hcl",
                  ligand_pdb="2vta",
                  ligand_identifier="LZ1")

    with MoleculeWriter("/home/pcurran/use_case/no_constraints/rescored.mol2") as w:
        for m in r.rescored_ligands:
            w.write(m)


if __name__ == "__main__":
    main()
