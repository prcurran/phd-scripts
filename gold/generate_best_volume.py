import os
import time
import pandas as pd

from hotspots.hs_io import HotspotWriter, HotspotReader
from hotspots.result import Extractor


def get_best_vol(path):
    """
    generate best volume

    :param str path: path to input
    :return:
    """

    hs = HotspotReader(path).read()

    start = time.time()
    extractor = Extractor(hs)

    return extractor.extract_best_volume(volume=350,
                                         pharmacophores=False)[0], time.time() - start


pdbs = ['3hl5', '2qd9', '1sj0', '3eml', '2i78', '1vso', '1njs']
base = "/local/pcurran/GOLD"

times = []
pid = []

for pdb in pdbs:
    print pdb
    path = os.path.join(base, pdb, "out.zip")
    hs, t = get_best_vol(path)
    pid.append(pdb)
    times.append(t)

    with HotspotWriter(os.path.join(os.path.dirname(path), "best_vol")) as w:
        w.write(hs)

df = pd.DataFrame({"time": times, "pdb": pid})
df.to_csv("best_volume_stats.csv")