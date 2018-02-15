"""
More information about the SuperStar method is available from:
    1) SuperStar: A knowledge-based approach for identifying interaction sites in proteins
    M. L. Verdonk, J. C. Cole and R. Taylor, J. Mol. Biol., 289, 1093-1108, 1999
    DOI: 10.1006/jmbi.1999.2809

    2)SuperStar: Improved knowledge-based interaction fields for protein binding sites
    M. L. Verdonk, J. C. Cole, P. Watson, V. Gillet and P. Willett, J. Mol. Biol., 307, 841-859, 2001
    DOI: 10.1006/jmbi.2001.4452

Tutorial for SuperStar is available from:
    https://www.ccdc.cam.ac.uk/support-and-resources/ccdcresources/superstar.pdf
    Chapter 4- Using SuperStar from the Command-line is particularly relevant

This code is based upon CCDC_internal API
"""

from template_strings import superstar_ins
from ccdc import io, utilities
import os
import sys
import glob
from ccdc.protein import Protein
from ccdc.utilities import Grid
import subprocess


class RunSuperstar:
    """
    class to handle SuperStar run
    """

    class Settings:
        """
        setting for Superstar run
        """
        def __init__(self):
            self.jobname = None
            self.probename = None
            self.moleculefile = None
            self.cavity_origin = None
            self.occulsionthreshold = 5
            self.mapbackgroundvalue = 1
            self.boxborder = 10
            self.minpropensity = 1
            self._superstar_executable = None
            self._superstar_env = None
            self.working_directory = None

    def __init__(self, **kw):
        settings = kw.get('settings')
        if settings is None:
            settings = self.Settings
        self.settings = settings

        main_dir = os.environ.get('MAINDIR')
        if main_dir:
            if sys.platform == 'win32':
                self.settings._superstar_executable = 'superstar_app.exe'
            else:
                self.settings._superstar_executable = ' '.join([
                    os.path.join(os.environ['MAINDIR'], 'run.sh'),
                    'superstar_app.x'
                ])
            self.settings._superstar_env = dict()
        else:
            if sys.platform == 'win32':
                base = os.path.dirname(io.csd_directory())
                merc = glob.glob(os.path.join(base, 'mercury*'))
                if len(merc):
                    merc = merc[0]
                self.settings._superstar_executable = os.path.join(merc, 'superstar_app.exe')
            elif sys.platform == 'darwin':
                self.settings._superstar_executable = os.path.join(
                    os.path.dirname(io.csd_directory()), 'mercury.app', 'Contents', 'MacOS', 'superstar'
                )
            else:
                self.settings._superstar_executable = os.path.join(os.path.dirname(io.csd_directory()), 'bin', 'superstar')
            self.settings._superstar_env = dict(
                SUPERSTAR_ISODIR=str(os.path.join(os.path.dirname(io.csd_directory()), 'isostar_files', 'istr')),
                SUPERSTAR_ROOT=str(os.path.join(os.path.dirname(io.csd_directory()), "Mercury"))
            )
        self.settings.working_directory = utilities._test_output_dir()
        print "WKDIR", " ", self.settings.working_directory

    def _append_cavity_info(self):
        """
        updates ins file with any cavity information

        :return: None
        """

        if self.settings.cavity_origin is not None:
            pnt = self.settings.cavity_origin
            extension = '\nCAVITY_ORIGIN {} {} {}'.format(pnt[0], pnt[1], pnt[2])
        else:
            extension = '\nSUBSTRUCTURE ALL'
        self.ins += extension

    def _get_inputs(self, out_dir):
        """
        assembles the ins files, uses a template string from template_strings.py

        :param out_dir: str, output directory
        :return: None
        """

        self.ins = superstar_ins(self.settings)
        self._append_cavity_info()
        self.fname = os.path.join(out_dir, "superstar_{}.ins".format(self.settings.jobname.split(".")[0]))
        w = open(self.fname, "w")
        w.write(self.ins)
        w.close()

    def run_superstar(self, prot, out_dir):
        """
        calls SuperStar as command-line subprocess

        :param prot: a :class:`ccdc.protein.Protein` instance
        :param out_dir: str, output directory
        :return:
        """

        with utilities.PushDir(self.settings.working_directory):
            if prot is not None:
                with io.MoleculeWriter('protein.pdb') as writer:
                    writer.write(prot)
            self._get_inputs(out_dir)
            env = os.environ.copy()
            env.update(self.settings._superstar_env)
            cmd = self.settings._superstar_executable + ' ' + self.fname
            subprocess.call(cmd, shell=sys.platform != 'win32', env=env)
        return SuperstarResult(self.settings)


class SuperstarResult:
    """
    stores a superstar result
    """

    def __init__(self, settings):
        self.settings = settings
        self.identifier = settings.jobname.split(".")[0]
        print self.identifier

        gridpath = os.path.join(self.settings.working_directory, self.identifier + ".ins.acnt")
        if os.path.exists(gridpath):
            self.grid = Grid.from_file(gridpath)
        else:
            raise AttributeError('{} superstar grid could not be found'.format(self.identifier))

        ligsitepath = os.path.join(self.settings.working_directory, self.identifier + ".ins.ligsite.acnt")
        if os.path.exists(gridpath):
            self.grid = Grid.from_file(ligsitepath)
        else:
            raise AttributeError('{} ligsite grid could not be found'.format(self.identifier))


def _run_ss(prot, out_dir, centroid = None, charged_probes=True):
    """
    initiates a SuperStar run for a given protein and probe

    :param prot: a :class:`ccdc.protein.Protein` instance
    :param out_dir: str, output directory
    :param centroid: tup, coordinates of cavity origin
    :param charged_probes: bool, if True 'positive' and 'negative' probes will be used
    :return: a :class:`SuperstarResult` instance
    """

    sites = []

    if charged_probes:
        probe_dict = dict(
            apolar='AROMATIC CH CARBON',
            donor='UNCHARGED NH NITROGEN',
            acceptor='CARBONYL OXYGEN',
            positive='CHARGED NH NITROGEN',
            negative='CHLORIDE ANION'
        )
    else:
        probe_dict = dict(
            apolar='AROMATIC CH CARBON',
            donor='UNCHARGED NH NITROGEN',
            acceptor='CARBONYL OXYGEN')

    for n, probe in probe_dict.items():
        if n == "positive":
            s = RunSuperstar()
            s.settings.jobname = "{}.ins".format(n)
            s.settings.probename = probe
            s.settings.moleculefile = "protein.pdb"
            s.settings.cavity_origin = centroid
            result = s.run_superstar(prot, out_dir)
            sites.append(result)
    return sites

if __name__ == "__main__":
    prot = Protein.from_file("protein.pdb")
    out_dir = "C:/Users/pcurran/Desktop/experiments/ss"
    sites = _run_ss(prot, out_dir)

    print sites[0].grid, sites[0].ligsite
