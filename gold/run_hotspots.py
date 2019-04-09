import os
import pandas as pd
import tarfile
from ccdc.protein import Protein
import tempfile
import time

from ccdc.utilities import Timer
from hotspots.calculation import Runner
from hotspots.hs_io import HotspotWriter


timer = Timer()


def get_protein(pdb, base, tmp):
    """
    get protein instance from zipped folder

    :param pdb:
    :param base:
    :return:
    """
    stem = os.path.join(base, pdb)
    zfile = [os.path.join(stem, f) for f in os.listdir(stem) if f.endswith(".gz")][0]
    f = tarfile.open(zfile)

    for tar_info in f:
        if tar_info.name.split("/")[1] == "receptor.pdb":
            with open(os.path.join(tmp, "receptor.pdb"), 'w') as w:
                out = f.extractfile(member=tar_info).read()
                w.write(out)
    f.close()
    return Protein.from_file(os.path.join(tmp, "receptor.pdb"))


@timer.decorate('hotspots')
def get_hotspot(protein):
    """
    generate hotspot result and time

    :param protein:
    :return:
    """
    start = time.time()
    r = Runner()
    settings = Runner.Settings()
    settings.rotations = 3000
    settings.sphere_maps = False
    result = r.from_protein(protein=protein,
                            charged_probes=False,
                            nprocesses=3,
                            cavities=None,
                            settings=settings
                            )
    return result, time.time() - start

def main():
    """
    main

    :return:
    """
    base = "/local/pcurran/GOLD"
    pdbs = set(pd.read_csv(os.path.join(base, "targets.csv"))['PDB'])
    tmp = tempfile.mkdtemp()

    pid = []
    times = []
    for pdb in pdbs:
        prot = get_protein(pdb, base, tmp)
        hs, time = get_hotspot(prot)
        pid.append(pdb)
        times.append(time)

        with HotspotWriter(os.path.join(base, pdb), zip_results=True) as w:
            w.write(hs)

    df = pd.DataFrame({'PDB': pid,
                       'Time': times})
    df.to_csv(os.path.join(base, "run_stats.csv"))
    timer.report()


if __name__ == "__main__":
    main()
