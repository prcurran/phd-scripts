#!/usr/bin/env python

from pymol import cmd
import os
from os import listdir
from os.path import isfile, join, splitext

###################################################################################################################
# RUN ME IN PYMOL
# Ligand overlay protocol
# Align based on sequence similarity. Write overlayed protein and ligand to file
###################################################################################################################

def initialise_pymol(align_to):
   cmd.reinitialize()
   # load the reference structure
   cmd.fetch(align_to)
   #cmd.load(r'Z:\PycharmProjects\LigandMaps\CK2_overlay\pdb_files\ar008_K74A_002_deposit.pdb',"protein")

def ligands_where_are_you(pdb, align_to, in_dir, out_dir):
    # load the pdb and align to reference 1AQL
    pdbe = join(in_dir, pdb + ".pdb")
    cmd.load(r"{}".format(pdbe),"{}".format(pdb))
    #cmd.align("""{} and chain A""".format(pdb), """{} and chain A""".format(align_to))
    cmd.align(pdb, align_to)
    cmd.save((os.path.join(out_dir,pdb + ".pdb")), pdb)
    cmd.delete(pdb)

def codes_from_filename(struc, mypath, mylist):
    #extract pdb_code from file name
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    for i in onlyfiles:
        if i == '{}.pdb'.format(struc):
            continue
        else:
            head, tail = splitext(i)
            mylist.append(head)

#input##################################################################################################################
align_to = "1aq1"
target = "CDK2"
########################################################################################################################

pdb_codes = []
in_dir = "C:/Users/pcurran/Desktop/paper1/{}/input".format(target)
out_dir = "C:/Users/pcurran/Desktop/paper1/{}/output".format(target)
codes_from_filename(align_to, in_dir, pdb_codes)

print pdb_codes

initialise_pymol(align_to)
for i, pdb in enumerate(pdb_codes):
    ligands_where_are_you(pdb, align_to, in_dir, out_dir)
