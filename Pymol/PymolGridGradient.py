from ccdc.utilities import Grid
from inferno import inferno_data
from os.path import join

class GridPoint(object):
    """
    class object constructing a sphere at each grid point over a fhm score of 5
    """

    def __init__(self, x, y, z, grd, color_ramp):
        self.value = grd.value(x,y,z)
        self.min = 5
        self.max = grd.extrema[1]
        self.red, self.green, self.blue = self.index(color_ramp)
        self.x, self.y, self.z = self._indices_to_point(x, y, z, grd)


    @staticmethod
    def _indices_to_point(i, j, k, g):
        """
        Return x,y,z coordinate for a given grid index

        :param i: int, indice (x-axis)
        :param j: int, indice (y-axis)
        :param k: int, indice (z-axis)
        :param g: a :class: `ccdc.utilities.Grid` instance
        :return: float(x), float(y), float(z)
        """

        ox, oy, oz, = g.bounding_box[0]
        gs = 0.5
        return ox + float(i) * gs, oy + float(j) * gs, oz + gs * float(k)

    def point_to_string(self):
        """

        :return:
        """
        str_template = '''cluster["apolar"] += [COLOR, {0},{1},{2}] + [SPHERE, {3},{4},{5}, 0.1]\n'''.format(self.red, self.green, self.blue,
                                                                         self.x, self.y, self.z)
        return str_template

    def index(self, inferno_data):
        """
        fetch color from ramp based on score
        """
        index = round(((self.value - 5) / (self.max - 5)) * (len(inferno_data)))
        if index > 255:
            index = 255
        print index
        return inferno_data[int(index)]


def header(fname):
    """
    prep output file
    :param fname:
    :return:
    """
    string = '''from pymol import cmd
from pymol.cgo import *

cluster = {"apolar": [], "acceptor": [], "donor": []}
'''

    with open(fname, "w") as w:
        w.write(string)


def footer(fname):
    """
    adjust output file
    :param fname:
    :return:
    """
    string = '''cmd.load_cgo(cluster["apolar"], "surface", 1)
cmd.set("transparency", 0.2, "surface")'''

    with open(fname, "a") as w:
        w.write(string)


def main():
    out_dir = "Z://fragment-hotspot-results//patel_set//2//out"
    grd = Grid.from_file(join(out_dir, "apolar.grd"))
    nx, ny, nz = grd.nsteps
    fname = join(out_dir, "inferno_gradient.py")
    header(fname)
    string_list = ""

    for a in range(nx):
        for b in range(ny):
            for c in range(nz):
                if grd.value(a,b,c) > 5:
                    pnt = GridPoint(a,b,c, grd, inferno_data)
                    string_list += pnt.point_to_string()
                else:
                    continue

    with open(fname, "a") as pymol_file:
        pymol_file.write(string_list)

    footer(fname)

if __name__ == "__main__":
    main()