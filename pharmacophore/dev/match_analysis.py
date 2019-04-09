
from fragment_hotspot_maps.pharmacophore import PharmacophoreModel
from os.path import join
from fragment_hotspot_maps.utilities import Utilities
from pprint import pprint


base = "Z:/fragment-hotspot-results/patel_set/CDK2/1hcl"

reference = PharmacophoreModel.from_file(join(base, "gold_standard/model.cm"))
hotspot_model = PharmacophoreModel.from_file(join(base, "hotspot_volume/1/pharmacophore.cm"))

for model_feat in hotspot_model.features:
    feat_by_dist = {Utilities.get_distance(ref_feat.feature_coordinates, model_feat.feature_coordinates): ref_feat
                    for ref_feat in reference.features if ref_feat.feature_type == model_feat.feature_type}
    if len(feat_by_dist) > 1:
        s = [value for (key, value) in sorted(feat_by_dist.items(), reverse=False)]
        shortest = s[0]
    else:
        shortest = feat_by_dist.values()[0]

    rmsd = Utilities.get_distance(model_feat.feature_coordinates, shortest.feature_coordinates)

    if model_feat.feature_type == "apolar":
        if rmsd < 4.0:
            print "MATCH", "type", model_feat.feature_type, "coordinates", model_feat.feature_coordinates, "rmsd", rmsd
        else:
            print "NO_MATCH", "type", model_feat.feature_type, "coordinates", model_feat.feature_coordinates, "rmsd", rmsd
    else:
        if rmsd < 2.0:
            print "MATCH", "type", model_feat.feature_type, "coordinates", model_feat.feature_coordinates, "rmsd", rmsd
        else:
            print "NO_MATCH", "type", model_feat.feature_type, "coordinates", model_feat.feature_coordinates, "rmsd", rmsd
