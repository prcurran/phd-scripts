'''
Ghecom source code is available from:
    http://strcomp.protein.osaka-u.ac.jp/ghecom/download_src.html

information about the method is available from:
    1)Kawabata T. (2010) Detection of multi-scale pockets on protein surfaces using mathematical morphology.
    Proteins,78, 1195-1121
    doi: 10.1002/prot.22639

    2)Kawabata T, Go N. (2007) "Detection of pockets on protein surfaces using small and large probe spheres
    to find putative ligand binding sites". Proteins, 68,516-529
    doi: 10.1002/prot.22639

This script runs ghecom pocket detection and converts the output to a CCDC grid to facillate further manipulation
using the csd-python-api
'''


##########################################################################
from os.path import join, dirname, abspath, isabs
from os import chdir, system

#from ccdc.utilities import Grid
from ccdc import io
from ccdc_internal.interaction import Grid
from ccdc.protein import Protein

##########################################################################

class GhecomInfo():
    '''
    Information and functions required to run ghecom and output a CCDC Grid object
    '''
    def __init__(self, fname = "protein.pdb"):

        if isabs(fname):
            self.fname = fname
        else:
            self.fname = abspath(fname)

        self.prot = Protein.from_file(fname)
        self._protein_preparation()
        self.out_dir = dirname(self.fname)
        self.ghecom_out = join(self.out_dir, "ghecom_out.pdb")
        self.ghecom_grid = self._initalise_grid()

        # add a more general way of getting to the "ghecom" cmd
        self.run_dir = "/home/pcurran/src"

    def _protein_preparation(self):
        ''' Remove water, ligands, metals and write out protein ready for call'''

        self.prot.remove_all_waters()

        for lig in self.prot.ligands:
            self.prot.remove_ligand(lig.identifier)
        self.prot.remove_all_metals()

        with io.MoleculeWriter(self.fname) as w:
            w.write(self.prot)

    def run_ghecom(self):
        '''Organising function, run command and return out CCDC Grid object'''

        command = "./ghecom {} -M M -gw 0.5 -rli 2.5 -rlx 9.5  -opoc {}".format(self.fname, self.ghecom_out)
        chdir(self.run_dir)
        system(command)
        self._get_ghecom_grid()

        return self.ghecom_grid

    def _initalise_grid(self, padding = 1):
        '''Create a blank CCDC grid using protein atom coordinates'''

        x = []
        y = []
        z = []

        for atm in self.prot.atoms:
            x.append(atm.coordinates.x)
            y.append(atm.coordinates.y)
            z.append(atm.coordinates.z)

        bl = (min(x) - padding, min(y) - padding, min(z) - padding)
        tr = (max(x) + padding, max(y) + padding, max(z) + padding)

        print "bottom left: {}, top right: {}".format(bl, tr)

        g = Grid(origin= bl, far_corner=tr, spacing = 0.5)
        return g

    def _point_to_indices(self, p, g):
        '''Return the nearest grid index for a given point. '''

        gs = 0.5
        rx, ry, rz = [round(i / gs) for i in p]
        ox, oy, oz = [round(i / gs) for i in g.bounding_box[0]]
        return int(rx - ox), int(ry - oy), int(rz - oz)

    def _get_lines_from_file(self):
        '''fetch lines from pdb file'''

        f = open(self.ghecom_out)
        lines = f.readlines()
        f.close()

        for i in range(0, len(lines)):
            lines[i] = lines[i].strip()
        return lines

    def _get_ghecom_grid(self):
        '''for every grid point, fetch coords and rinacc value and set to blank grid'''

        input = self._get_lines_from_file()

        for line in input:
            if line.startswith("HETATM"):
                coords = (float(line[31:38]), float(line[39:46]), float(line[47:54]))
                rinacc = float(line[61:66])
                i, j, k = self._point_to_indices(coords, self.ghecom_grid)
                x, y, z = self.ghecom_grid.nsteps

                if i < x and i > 0 and j < y and j > 0 and k < z and k > 0:
                    self.ghecom_grid.set_value(i, j, k, 10 - rinacc)
                    # FYI rinacc score inverted to match LIGSITE output

##########################################################################

def main():
    '''
    supply relative or absoulte path
    ghecom_out.pdb will write to the same dir as import protein
    run_ghecom returns a grid
    '''

    g = GhecomInfo(fname = "protein.pdb")
    ghecom = g.run_ghecom()
    print ghecom

    ghecom.write("/home/pcurran/expt/ghecom.grd")

if __name__ == "__main__":
    main()