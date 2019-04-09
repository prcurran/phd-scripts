#!/usr/bin/env python
#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2018-03-01: created by the Cambridge Crystallographic Data Centre
#
'''
    binding_site_superposition.py    -   superpose protein chains or binding sites and write mol2 file

A number of target chains are superposed to a reference chain using a CCDC implementation of the Smith-Waterman algorithm. The target chains can be specified by the user
or automatically selected from a file containing a list of protein sequences based on their identity and similarity to the reference. A sequence alignment tool (e.g. FASTA) is used
to calculate the sequence similarity and identity with respect to the reference chain. 3D structural information of the top N sequences is then retrieved from the wwPDB and used to
perform a geometrical alignment of the entire chain or of the binding site around a reference ligand onto the reference chain.

This code uses biopython (https:https://biopython.org/)  as an interface to the sequence alignment tool which may need installing before use.

For more help use the -h argument.

'''
########################################################################################

from __future__ import division, absolute_import, print_function

import os
import argparse
import urllib
import tempfile
import shutil
import warnings
from subprocess import call
from csv import DictWriter

from Bio import SeqIO
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO

from ccdc.utilities import _temporary_copy
from ccdc.io import MoleculeWriter
from ccdc.protein import Protein
from ccdc.entry import Entry


class Runner(argparse.ArgumentParser):
    '''Defines and parses arguments, runs the extraction according to the arguments.'''

    def __init__(self):
        '''Defines and parses the arguments.'''

        usage_txt = "Superpose protein-ligand binding sites based on their protein chain sequence."
        super(self.__class__, self).__init__(description=usage_txt,
                                             formatter_class=argparse.RawDescriptionHelpFormatter)

        self.add_argument(
            'reference_chain',
            help='Reference protein chain, e.g. 3a7h:B',
        )
        self.add_argument(
            '-lig', '--reference_ligand', nargs='?',
            help='Reference ligand defining binding site to superpose, e.g. B:ATP400.'
                 'If None, whole chain will be superposed',
            default=None
        )
        group = self.add_mutually_exclusive_group(required=True)
        group.add_argument(
            '-c', '--target_chains', nargs='?',
            help='Protein chains to superpose, e.g. 3a7i:A,3a7j:A,4qms:A,4qo9:B,4qmv:A',
            default=None
        )
        group.add_argument(
            '-seq', '--sequence_file', nargs='?',
            help='Protein chain sequence file to search',
            default=None
        )
        self.add_argument(
            '-max', '--maximum_sequences', nargs='?',
            help='Maximum number of protein chain sequences to superpose',
            type=int,
            default=25
        )
        self.add_argument(
            '-atoms', '--superposition_atoms', nargs='?',
            help='Atoms to superpose',
            choices=['rigid', 'backbone', 'calpha'],
            default='rigid'
        )
        self.add_argument(
            '-o', '--output_file', nargs='?',
            help='Binding site superposition output file',
            default='superposition.mol2'
        )
        self.add_argument(
            '-s', '--superposition_output', nargs='?',
            help='Superposition output type',
            choices=['binding_site', 'full_protein'],
            default=None
        )
        self.add_argument(
            '-csv', '--csv_summary_file', nargs='?',
            help='CSV summary file',
            default=None
        )
        group = self.add_mutually_exclusive_group()
        group.add_argument(
            '--pdb_mirror', nargs='?',
            help='PDB mirror file path template, e.g. /local/mirror/pub/pdb/data/structures/all/pdb/pdb{}s.ent.gz',
            default=None
        )
        group.add_argument(
            '--pdb_url', nargs='?',
            help='PDB URL file template',
            default='http://www.rcsb.org/pdb/files/{}.pdb'
        )
        self.add_argument(
            '-r',
            '--radius_around_ligand',
            help='Radius around ligand to define binding site [6.0]',
            type=float,
            default=6.0
        )
        self.add_argument(
            '-search', '--sequence_search_tool',
            help='Sequence search tool, e.g. FASTA',
            default=None
        )
        self.add_argument(
            '-align', '--sequence_alignment_tool',
            help='Sequence alignment tool, e.g. FASTA',
            default=None
        )
        self.add_argument(
            '-Mi',
            '--maximum_identity',
            help='Maximum identity',
            type=float,
            default=100.0
        )
        self.add_argument(
            '-mi',
            '--minimum_identity',
            help='Minimum identity',
            type=float,
            default=0.0
        )
        self.add_argument(
            '-Ms',
            '--maximum_similarity',
            help='Maximum similarity',
            type=float,
            default=100.0
        )
        self.add_argument(
            '-ms',
            '--minimum_similarity',
            help='Minimum similarity',
            type=float,
            default=0.0
        )

        self.args = self.parse_args()
        self.output_dir = tempfile.mkdtemp()

        self.reference_id = self.args.reference_chain.split(':')[0]
        self.reference_chain_id = self.args.reference_chain.split(':')[1]
        self.reference_pdb_file = os.path.join(self.output_dir, 'pdb{}.pdb').format(self.reference_id)
        self.reference_protein = self.load_reference_protein()
        (self.reference_ligand, self.reference_binding_site) = self.load_reference_ligand()

        self.superposition_output = self.args.superposition_output

        if self.reference_ligand is None and self.superposition_output == 'binding_site':
            self.error('--reference_ligand argument required for --superposition_output binding_site')

        if self.superposition_output is None:
            if self.reference_binding_site is not None:
                self.superposition_output = 'binding_site'
            else:
                self.superposition_output = 'full_protein'

    def load_reference_protein(self):
        self.fetch_pdb(self.reference_id, self.reference_pdb_file)
        reference_protein = Protein.from_file(self.reference_pdb_file)
        print('Reference protein {} chain {} '.format(
            reference_protein.identifier,
            self.reference_chain_id))
        return reference_protein

    def load_reference_ligand(self):
        reference_binding_site = None
        reference_ligand = None
        if self.args.reference_ligand is not None:

            reference_ligand = next((ligand for ligand in self.reference_protein.ligands
                                     if ligand.identifier == self.args.reference_ligand), None)
            if reference_ligand is None:
                raise IndexError('Ligand not found for ligand identifier {}'.format(self.args.reference_ligand))

        if reference_ligand is not None:
            reference_binding_site = Protein.BindingSiteFromMolecule(
                self.reference_protein, reference_ligand, self.args.radius_around_ligand
            )
            print('Reference ligand {}'.format(reference_ligand.identifier))

        return reference_ligand, reference_binding_site

    def load_target(self, target_pdb_id):
        '''Load target protein chain

        :param target_pdb_id: The pdb identifier of the target
        '''
        target_pdb_file = os.path.join(self.output_dir, 'pdb{}.pdb').format(target_pdb_id)
        self.fetch_pdb(target_pdb_id, target_pdb_file)
        target_protein = Protein.from_file(target_pdb_file)

        return target_protein

    def fetch_pdb_from_url(self, pdb_id, filename):
        '''Fetch pdb file from PDB URL

        :param pdb_id The pdb identifier
        :param filename: The output filename
        '''
        url = self.args.pdb_url.format(pdb_id)
        pdb_file = urllib.urlopen(url).read()
        with open(filename, 'w') as fh:
            fh.write(pdb_file)

    def fetch_pdb_from_mirror(self, pdb_id, filename):
        '''Fetch pdb file from PDB mirror

        :param pdb_id The pdb identifier
        :param filename: The output filename
        '''
        full_path_input_file = self.args.pdb_mirror.format(pdb_id)
        with _temporary_copy(full_path_input_file) as pdb_file:
            shutil.copy2(pdb_file, filename)

    def fetch_pdb(self, pdb_id, filename):
        if self.args.pdb_mirror is not None:
            self.fetch_pdb_from_mirror(pdb_id, filename)
        else:
            self.fetch_pdb_from_url(pdb_id, filename)

    @staticmethod
    def fix_relibase_id(chain_id):
        '''Convert Relibase chain id to CSD Python API chain id'''
        chain_id = chain_id.replace('CHN::', '')
        chain_id = chain_id.replace('pdb', '')
        chain_id = chain_id.replace('-', ':')
        return chain_id

    @staticmethod
    def extract_ids(target_id):
        '''Extract pdb identifier and chain identifier from target identifier in a sequence file

        :param target_id: The input target identifier
        :returns: tuple of pdb identifier and chain identifier
        '''
        target_chain = Runner.fix_relibase_id(target_id)

        for separator in ['|', ' ']:
            target_chain = target_chain.split(separator)[0]

        for separator in [':', '_']:
            if separator in target_chain:
                return target_chain.split(separator)

        return target_chain, ''

    def sequence_search(self):
        '''Perform sequence search

        :returns: a list of target chains
        '''
        sequence = None
        # Write sequence file for input chain
        records = list(SeqIO.parse(self.reference_pdb_file, 'pdb-seqres'))
        for record in records:
            if ':' in record.id:
                (record_pdb_id, record_chain_id) = record.id.split(':')
            else:
                (record_pdb_id, record_chain_id) = ('', record.id)
            if record_pdb_id.lower() == self.reference_id.lower() and record_chain_id == self.reference_chain_id:
                sequence = record
                break

        if sequence is None:
            raise IndexError('Sequence not found for protein chain {}'.format(self.reference_chain_id))

        # Write query sequence file
        with open(os.path.join(self.output_dir, 'query.fasta'), 'w') as output_handle:
            SeqIO.write(sequence, output_handle, 'fasta')

        # Run FASTA '-b', str(self.args.maximum_sequences),
        args = ['-q', '-m', '10',
                os.path.join(self.output_dir, 'query.fasta'), self.args.sequence_file]

        # Read alignment output
        with open(os.path.join(self.output_dir, 'output.fasta'), 'w') as out:
            call([self.args.sequence_search_tool] + args, stdout=out)

        target_chains = []

        query_result = SearchIO.read(os.path.join(self.output_dir, 'output.fasta'), "fasta-m10")
        for hit in query_result:
            (target_pdb_id, target_chain_id) = self.extract_ids(hit.id)

            hsp = hit[0]
            if target_pdb_id.lower() == self.reference_id.lower() and target_chain_id == self.reference_chain_id:
                continue
            if hsp.ident_pct < self.args.minimum_identity or hsp.ident_pct > self.args.maximum_identity:
                continue
            if hsp.pos_pct < self.args.minimum_similarity or hsp.pos_pct > self.args.maximum_similarity:
                continue

            target_chains.append({'pdb_id': target_pdb_id,
                                  'chain_id': target_chain_id,
                                  'identity': '{:.2f}'.format(hsp.ident_pct),
                                  'similarity': '{:.2f}'.format(hsp.pos_pct)})

            if len(target_chains) >= self.args.maximum_sequences:
                break

        return target_chains

    @staticmethod
    def write_csv_summary_file(csv_summary_file, target_chains):
        '''Write a summary of the superposition results in CSV format

        :param csv_summary_file: The output csv filename
        :param target_chains: A list of target chain dictionaries
        '''
        with open(csv_summary_file, 'wb') as csv_file:
            field_names = ['pdb_id', 'chain_id', 'identity', 'similarity', 'rmsd']
            csv_writer = DictWriter(csv_file, fieldnames=field_names)
            csv_writer.writeheader()
            for target_chain in target_chains:
                csv_writer.writerow(target_chain)

    def run(self):
        '''Search for target chains and superpose target chains'''
        target_chains = []
        if self.args.sequence_file is not None:
            if self.args.sequence_search_tool is None:
                self.error('--sequence_search_tool argument required with --sequence_file')
            target_chains = self.sequence_search()
        elif self.args.target_chains is not None:
            for target_chain in self.args.target_chains.split(','):
                (target_pdb_id, target_chain_id) = target_chain.split(':')
                target_chains.append({'pdb_id': target_pdb_id,
                                      'chain_id': target_chain_id,
                                      'identity': 0,
                                      'similarity': 0})

        with MoleculeWriter(self.args.output_file) as mol_writer:

            if self.superposition_output == 'full_protein':
                mol_writer.write(self.reference_protein)
            else:
                mol_writer.write(self.reference_binding_site)

            for target_chain in target_chains:
                try:
                    self.run_one(target_chain, mol_writer)
                except RuntimeError as exception:
                    print("Superposition failed: {}".format(exception.args[0]))

        if self.args.csv_summary_file is not None:
            self.write_csv_summary_file(self.args.csv_summary_file, target_chains)

    def run_one(self, target_chain, mol_writer):
        '''Superpose a target chain onto the reference chain

        :param target_chain: The target chain dictionary
        :param mol_writer: Superposition file writer
        '''
        target_protein = self.load_target(target_chain['pdb_id'])

        print('Target protein {} chain {} identity {}% similarity {}%'.format(
            target_protein.identifier, target_protein[target_chain['chain_id']].identifier,
            target_chain['identity'], target_chain['similarity']))

        binding_site_superposition = Protein.ChainSuperposition()
        if self.args.sequence_alignment_tool is not None:
            binding_site_superposition.settings.sequence_search_tool = self.args.sequence_search_tool
            binding_site_superposition.settings.sequence_alignment_tool = self.args.sequence_alignment_tool
        binding_site_superposition.settings.superposition_atoms = self.args.superposition_atoms

        if self.reference_binding_site is not None:
            (bs_rmsd, bs_transformation) = binding_site_superposition.superpose(
                self.reference_protein[self.reference_chain_id],
                target_protein[target_chain['chain_id']],
                self.reference_binding_site)
            target_chain['rmsd'] = '{:.4f}'.format(bs_rmsd)
            print('Binding site RMSD is: {}'.format(target_chain['rmsd']))
            print('Transformation matrix is: {}'.format(bs_transformation))

        else:
            (chain_rmsd, chain_transformation) = binding_site_superposition.superpose(
                self.reference_protein[self.reference_chain_id],
                target_protein[target_chain['chain_id']])
            target_chain['rmsd'] = '{:.4f}'.format(chain_rmsd)
            print('Chain RMSD is: {}'.format(target_chain['rmsd']))
            print('Transformation matrix is: {}'.format(chain_transformation))

        if self.superposition_output == 'full_protein':
            target_molecule = target_protein
        else:
            target_molecule = Protein.BindingSiteFromMolecule(
                target_protein, self.reference_ligand, self.args.radius_around_ligand
            )

        target_entry = Entry.from_molecule(target_molecule,
                                           rmsd='{}'.format(target_chain['rmsd']),
                                           identity='{}%'.format(target_chain['identity']),
                                           similarity='{}%'.format(target_chain['similarity']))
        mol_writer.write_entry(target_entry)

########################################################################################

if __name__ == '__main__':
    Runner().run()
