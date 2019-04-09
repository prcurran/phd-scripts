from fragment_hotspot_maps.pharmacophore import PharmacophoreModel
from pprint import pprint

target = "CDK2"
pdb = "1aq1"


# fname = "Z:/fragment-hotspot-results/patel_set/{}/{}/gold_standard/model.cm".format(target, pdb)
# pm = PharmacophoreModel.from_file(fname)
#
# #pm.write("Z:/fragment-hotspot-results/patel_set/{}/{}/gold_standard/gold_standard.py".format(target, pdb))
#
# pm.write("Z:/fragment-hotspot-results/patel_set/{}/{}/gold_standard/test.cm".format(target, pdb))

fname = "Z:/fragment-hotspot-results/patel_set/{}/{}/hotspot_volume/0/pharmacophore.cm".format(target,pdb)
pm = PharmacophoreModel.from_file(fname)
pm.write("Z:/fragment-hotspot-results/patel_set/{}/{}/hotspot_volume/0/p.cm".format(target,pdb))