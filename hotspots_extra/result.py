
    # def tractability_map(self):
    #     """
    #     generate the best volume and labels with the median value. A median > 14 is more likely to be tractable
    #
    #     :return: a :class:`hotspots.result.Results` instance
    #     """
    #     extractor_settings = Extractor.Settings()
    #     extractor_settings.cutoff = 5
    #     extractor_settings.island_max_size = 500
    #
    #     extractor = Extractor(self, settings=extractor_settings)
    #     extractor.extract_best_volume(volume=500)
    #     # hist = extractor.extracted_hotspots[0].map_values()
    #     #
    #     # all_points = []
    #     # for x in hist.values():
    #     #     all_points += x.flatten().tolist()
    #     #
    #     # all_points = all_points[all_points != 0]
    #     # print(all_points)
    #     best_vol = extractor.extracted_hotspots[0]
    #     best_vol.identifier = best_vol.score()
    #
    #     return best_vol
    #
    # def all_tractability_maps(self):
    #     """
    #     generate the best volume and labels with the median value. A median > 14 is more likely to be tractable
    #
    #     :return: a :class:`hotspots.result.Results` instance
    #     """
    #     extractor_settings = Extractor.Settings()
    #     extractor_settings.cutoff = 5
    #     extractor_settings.island_max_size = 500
    #
    #     extractor = Extractor(self, settings=extractor_settings)
    #     extractor.extract_all_volumes(volume=500)
    #     extracted = []
    #     for cav in extractor.extracted_hotspots:
    #         hist = cav.map_values()
    #         all_points = []
    #         for x in hist.values():
    #             all_points += x.flatten().tolist()
    #
    #         all_points = all_points[all_points != 0]
    #         best_vol = cav
    #         best_vol.identifier = np.median(all_points)
    #         extracted.append(best_vol)
    #
    #     return extracted


    # def histogram(self, fpath="histogram.png"):
    #     """
    #     get histogram of zero grid points for the Fragment Hotspot Result
    #
    #     :param fpath: path to output file
    #     :return: data, plot
    #
    #     >>> result
    #     <hotspots.result.Results object at 0x000000001B657940>
    #     >>> plt = result.histogram()
    #     >>> plt.show()
    #     """
    #     data, plt = Figures.histogram(self)
    #     plt.savefig(fpath)
    #     return data, plt

    # def get_2D_diagram(self, ligand, fpath="diagram.png", title=False):
    #     """
    #     broken
    #     :param ligand:
    #     :param fpath:
    #     :param title:
    #     :return:
    #     """
    #     Figures._2D_diagram(hr, ligand, title=False, output="diagram.png")

    # def _get_superstar_profile(self, feature_radius=1.5, nthreads=6, features=None, best_volume=None):
    #     """
    #     *experimental feature*
    #
    #     enable calculation to different superstar probes at hotspot features. Perhaps a better understanding
    #     of the nature of each feature can be gained from doing this or perhaps it just adds noise.
    #
    #     :return:
    #     """
    #     # set additional object properties
    #     if features:
    #         self.features = features
    #     else:
    #         self.features = self._get_features(threshold=5, min_feature_size=6)
    #
    #     if best_volume:
    #         self.best_volume = best_volume
    #     else:
    #         self.best_volume = Grid.get_single_grid(self.super_grids, mask=False)
    #
    #     self.feature_spheres = self.best_volume.copy_and_clear()
    #     for feat in self.features:
    #         self.feature_spheres.set_sphere(point=feat.feature_coordinates,
    #                                         radius=feature_radius,
    #                                         value=1,
    #                                         scaling="None"
    #                                         )
    #
    #     # superstar run
    #     centroid = [self.best_volume.centroid()]
    #     a = _AtomicHotspot()
    #     a.settings.atomic_probes = ["carbonyl_oxygen", "carboxylate", "pyramidal_r3n", "water_oxygen"]
    #
    #     self.superstar_result = a.calculate(protein=self.protein,
    #                                         nthreads=nthreads,
    #                                         cavity_origins=centroid)
    #
    #     self.ss = []
    #
    #     # find overlap
    #     for r in self.superstar_result:
    #         common_spheres, common_result = Grid.common_grid([self.feature_spheres, r.grid])
    #         r.grid = (common_spheres & common_result) * common_result
    #
    #     # assign island to Hotspot Feature
    #     feat_id = []
    #     ss_id = []
    #     score = []
    #     import pandas as pd
    #
    #     for i, feat in enumerate(self.features):
    #
    #         for r in self.superstar_result:
    #             feat_id.append(i)
    #             ss_id.append(r.identifier)
    #
    #             ss_dict = {Helper.get_distance(feat.feature_coordinates, island.centroid()): island
    #                        for island in r.grid.islands(threshold=1)
    #                        if Helper.get_distance(feat.feature_coordinates, island.centroid()) < 1}
    #
    #             if len(ss_dict) == 0:
    #                 g = r.grid.copy_and_clear()
    #
    #             else:
    #                 shortest = sorted([f[0] for f in ss_dict.items()], reverse=False)[0]
    #                 g = ss_dict[shortest]
    #
    #             feat.superstar_results.append(_AtomicHotspotResult(identifier=r.identifier,
    #                                                                grid=g,
    #                                                                buriedness=None)
    #                                           )
    #
    #             score.append(g.grid_score(threshold=1, percentile=50))
    #
    #     return pd.DataFrame({"feature_id": feat_id, "interaction": ss_id, "score": score})


    # def _ngl_widget(self, out_dir=None):
    #     """
    #     jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000000.0
    #     creates ngl widget from hotspot. For use in ipython notebooks
    #     :param str out_dir:
    #     :return:
    #     """
    #     import nglview as nv
    #     from ipywidgets import IntSlider, interact
    #
    #     color_dict = {"apolar": "yellow",
    #                   "donor": "blue",
    #                   "acceptor": "red",
    #                   "negative": "magenta",
    #                   "positive": "cyan"}
    #     if out_dir:
    #         out = Helper.get_out_dir(out_dir)
    #     else:
    #         out = tempfile.mkdtemp()
    #
    #     for p, g in self.super_grids.items():
    #         g.write(join(out, "{}.ccp4".format(p)))
    #
    #     with MoleculeWriter(join(out, "protein.pdb")) as w:
    #         w.write(self.protein)
    #
    #     view = nv.NGLWidget()
    #     view.add_component(join(out, "protein.pdb"))
    #
    #     k = self.super_grids.keys()
    #     for i, p in enumerate(k):
    #         view.add_component(join(out, "{}.ccp4".format(p)))
    #         view.add_representation('isosurface', component=i + 1)
    #         view.update_representation(component=i + 1, color=color_dict[p])
    #
    #     @interact(x=IntSlider(description="HS Score", min=0, max=30, step=1))
    #     def f(x):
    #         view.update_representation(component=1, isolevel=int(x), isoleveltype='value')
    #         view.update_representation(component=2, isolevel=int(x), isoleveltype='value')
    #         view.update_representation(component=3, isolevel=int(x), isoleveltype='value')
    #         view.update_representation(component=4, isolevel=int(x), isoleveltype='value')
    #         view.update_representation(component=5, isolevel=int(x), isoleveltype='value')
    #
    #     return view


#
# class Extractor2(object):
#     """
#     A class to handle the extraction of molecular volumes from a Fragment Hotspot Map result
#
#     :param `hotspots.HotspotResults` hr: A Fragment Hotspot Maps result
#     :param `hotspots.Extractor.Settings` settings: Extractor settings
#     """
#
#     class Settings(object):
#         """
#         Default settings for hotspot extraction
#
#         :param float volume: required volume (default = 150)
#         :param float cutoff: only features above this value are considered (default = 14)
#         :param float spacing: grid spacing, (default = 0.5)
#         :param int min_feature_gp: the minimum number of grid points required to create a feature (default = 5)
#         :param int max_features: the maximum number of features in a extracted volume (default = 10)(not recommended, control at pharmacophore)
#         :param float min_distance: the minimum distance between two apolar interaction peaks (default = 6)
#         :param int island_max_size: the maximum number of grid points a feature can take. (default = 100)(stops overinflation of polar features)
#         :param bool pharmacophore: if True, generate a Pharmacophore Model (default = True)
#
#         """
#
#         def __init__(self, volume=150, cutoff=14, spacing=0.5, min_feature_gp=5, max_features=10, min_distance=6,
#                      island_max_size=100, pharmacophore=True):
#             self.volume = volume
#             self.cutoff = cutoff
#             self.spacing = spacing
#             self.min_feature_gp = min_feature_gp
#             self.max_features = max_features
#             self.min_distance = min_distance
#             self.island_max_size = island_max_size
#             self.pharmacophore = pharmacophore
#             self.mode = None
#             self.mvon = False
#
#         @property
#         def _num_gp(self):
#             """
#             number of grid point for a given volume
#             :return:
#             """
#             return int(float(self.volume) / self.spacing ** 3)
#
#         @property
#         def _search_radius(self):
#             """
#             describes search radius around a given seed
#             :return:
#             """
#             s = 3
#             s += round((int(self.volume) / 50))
#             print('search_radius', s)
#             return s
#
#     class _Optimiser(object):
#         """
#         A class to handle the optimisation operations
#
#         :param mask:
#         :param settings:
#         :param peak:
#         """
#
#         def __init__(self, mask, settings, peak=None):
#             self.peak = peak
#             self.mask = mask
#             self.settings = settings
#
#         def _count_island_points(self, threshold):
#             """
#             For a given island, the difference between the target number of grid points and the actual number of
#              grid points is returned
#             :param threshold:
#             :return: int
#             """
#             island = self.mask.get_best_island(threshold, mode="count", peak=self.peak)
#             if island is None:
#                 return 999999
#             points = (island > threshold).count_grid()
#             return abs(self.settings._num_gp - points)
#
#         def _raw_count_island_points(self, threshold):
#             """
#             For a given island, the difference between the target number of grid points and the actual number of
#              grid points is returned
#
#              -ve = too many points
#              +ve = too few points
#
#             :param threshold:
#             :return: int
#             """
#             island = self.mask.get_best_island(threshold, mode="count", peak=self.peak)
#             if island is None:
#                 return 999999
#             points = (island > threshold).count_grid()
#             return 800 - points
#
#         def _count_grid_points(self, threshold):
#             """
#             For a given island, the difference between the target number of grid points and the actual number of
#              grid points is returned
#             :param threshold:
#             :return: int
#             """
#             points = (self.top_island > threshold).count_grid()
#             return abs(self.settings._num_gp - points)
#
#         def _grow(self, best_island, tolerance=0.1):
#             """
#             *experimental*
#             fixes joining islands problem, this has not been tested use with caution :) .
#
#
#             :param best_island:
#             :return:
#             """
#
#             inner = self.mask.common_boundaries(best_island)
#             num_gp = inner.count_grid()
#             print((self.settings._num_gp - num_gp) / self.settings._num_gp)
#             grown = Grid.grow(inner, self.mask)     # adding points above 80th percentile as default,
#                                                     # this hasn't been played with
#             while ((self.settings._num_gp - num_gp) / self.settings._num_gp) > tolerance:
#                 # iterate over the growth cycle until target is reached
#                 grown = Grid.grow(inner, self.mask)
#                 diff = grown > inner
#
#                 # if diff.count_grid() < 10:
#                 #     break
#
#                 inner = grown
#                 num_gp = inner.count_grid()
#                 print(num_gp, 'out of', self.settings._num_gp)
#
#             tmp_best_island = inner * self.mask
#             g_vals = tmp_best_island.grid_values()
#             g_vals[::-1].sort()
#
#             try:
#                 threshold = g_vals[self.settings._num_gp]
#             except IndexError:
#                 threshold = g_vals.min()
#
#             return threshold, grown
#
#         def optimize_island_threshold(self, tolerance=0.1, start_threshold=25):
#             """
#             TEST
#             1) step down from 25 find a pair of best islands
#             finds the island threshold for a grid which returns the desired volume
#
#             :return:
#             """
#             from skimage.morphology import reconstruction
#
#             # find inital fragment volume hard coded as 800
#             while self._raw_count_island_points(threshold=start_threshold) > 0:
#                 start_threshold -= 1
#
#             large_island = self.mask.get_best_island(threshold=start_threshold, mode='count', peak=self.peak)
#             small_island = self.mask.get_best_island(threshold=start_threshold + 1, mode='count', peak=self.peak)
#             reconstruct_pts = 0
#
#             while reconstruct_pts < self.settings._num_gp:
#                 # ensure same dimensions (precaution: not sure if required)
#                 c_small_island = large_island.common_boundaries(small_island)
#
#                 seed = c_small_island.get_array()
#                 mask = large_island.get_array()
#                 # flood connected peaks
#                 reconstruct = reconstruction(seed=seed, mask=mask, method='dilation')
#                 small_island = Grid.array_to_grid(reconstruct, large_island.copy_and_clear())
#
#                 reconstruct_pts = ((small_island > 0) * small_island).count_grid()
#                 print("pts:", reconstruct_pts, "target:", self.settings._num_gp, "threshold:", start_threshold)
#                 start_threshold -= 1
#                 large_island = self.mask.get_best_island(threshold=start_threshold, mode='count', peak=self.peak)
#
#             print([i.extrema[1] for i in small_island.islands(threshold=18)])
#             return start_threshold, small_island
#
#             # large_island.write("/home/pcurran/reconstruct/large.grd")
#             # c_small_island.write("/home/pcurran/reconstruct/c_small_island.grd")
#             # g.write("/home/pcurran/reconstruct/reconstruct.grd")
#
#
#
#             # ret = optimize.minimize_scalar(self._count_island_points,  method='brent', tol=0.0001,
#             #                                options={'disp': True,
#             #                                         'xtol': 0.000001}
#             #                                )
#             # threshold = ret.x
#             # if threshold > 48:
#             #     threshold = 1
#             # best_island = self.mask.get_best_island(threshold=threshold, mode='count', peak=self.peak)
#             # print("target = {}, actual = {}".format(self.settings._num_gp, best_island.count_grid()))
#             #
#             # # If threshold is close to zero, keep all grid points
#             # try:
#             #     best_island = (best_island > threshold) * best_island
#             # except TypeError:
#             #     best_island = self.mask
#             #
#             # while ((self.settings._num_gp - best_island.count_grid()) / self.settings._num_gp) < 0:
#             #     # ensure underestimation
#             #     threshold += 0.1
#             #     best_island = (best_island > threshold) * best_island
#             #
#             # if ((self.settings._num_gp - best_island.count_grid()) / self.settings._num_gp) > tolerance:
#             #     print("Percentage error=", abs((self.settings._num_gp - best_island.count_grid() ) / self.settings._num_gp) * 100)
#             #     threshold, best_island = self._grow(best_island, tolerance=tolerance)
#             #
#             # print("target = {}, actual = {}".format(self.settings._num_gp, best_island.count_grid()))
#             #
#             #
#             # return threshold, best_island
#
#         # def optimize_island_threshold(self, tolerance=0.1):
#         #     """
#         #     finds the island threshold for a grid which returns the desired volume
#         #
#         #     :return:
#         #     """
#         #     # set options={'disp': True} for debugging
#         #     ret = optimize.minimize_scalar(self._count_island_points,  method='brent', tol=0.0001,
#         #                                    options={'disp': True,
#         #                                             'xtol': 0.000001}
#         #                                    )
#         #     threshold = ret.x
#         #     if threshold > 48:
#         #         threshold = 1
#         #     best_island = self.mask.get_best_island(threshold=threshold, mode='count', peak=self.peak)
#         #     print("target = {}, actual = {}".format(self.settings._num_gp, best_island.count_grid()))
#         #
#         #     # If threshold is close to zero, keep all grid points
#         #     try:
#         #         best_island = (best_island > threshold) * best_island
#         #     except TypeError:
#         #         best_island = self.mask
#         #
#         #     while ((self.settings._num_gp - best_island.count_grid()) / self.settings._num_gp) < 0:
#         #         # ensure underestimation
#         #         threshold += 0.1
#         #         best_island = (best_island > threshold) * best_island
#         #
#         #     if ((self.settings._num_gp - best_island.count_grid()) / self.settings._num_gp) > tolerance:
#         #         print("Percentage error=", abs((self.settings._num_gp - best_island.count_grid() ) / self.settings._num_gp) * 100)
#         #         threshold, best_island = self._grow(best_island, tolerance=tolerance)
#         #
#         #     print("target = {}, actual = {}".format(self.settings._num_gp, best_island.count_grid()))
#         #     return threshold, best_island
#
#     def __init__(self, hr, settings=None):
#         if settings is None:
#             self.settings = self.Settings()
#         else:
#             self.settings = settings
#         self._single_grid = None
#         self._masked_dic = None
#         self.out_dir = None
#         self.extracted_hotspots = None
#
#         if self.settings.mvon is True:
#             hr.super_grids.update({probe: g.max_value_of_neighbours() for probe, g in hr.super_grids.items()})
#
#         try:
#             hr.super_grids["negative"] = hr.super_grids["negative"].deduplicate(hr.super_grids["acceptor"],
#                                                                                 threshold=10,
#                                                                                 tolerance=2)
#
#             hr.super_grids["positive"] = hr.super_grids["positive"].deduplicate(hr.super_grids["donor"],
#                                                                                 threshold=10,
#                                                                                 tolerance=2)
#         except KeyError:
#             pass
#
#         try:
#             hr.super_grids.update({probe: g.minimal() for probe, g in hr.super_grids.items()})
#         except RuntimeError:
#             pass
#
#         self.hotspot_result = hr
#         self._masked_dic, self._single_grid = Grid.get_single_grid(self.hotspot_result.super_grids)
#
#     def extract_best_volume(self, volume="125", pharmacophores=True):
#         """
#         from the main Fragment Hotspot Map result, the best continuous volume is returned
#
#         :param float volume: volume in Angstrom^3
#         :param bool pharmacophores: if True, generates pharmacophores
#         :return: a `hotspots.result.Results` instance
#
#
#         >>> result
#         <hotspots.result.Results object at 0x000000001B657940>
#
#         >>> from hotspots.result import Extractor
#
#         >>> extractor = Extractor(result)
#         >>> best = extractor.extract_best_volume(volume=400)
#         [<hotspots.result.Results object at 0x0000000028E201D0>]
#
#         """
#         self.settings.volume = volume
#         self.settings.pharmacophore = pharmacophores
#         self.out_dir = None
#         self._peaks = None
#         self.settings.mode = "global"
#
#         print("???????", self.single_grid.extrema[1])
#         e = self._from_hotspot(self.single_grid,
#                                self.masked_dic,
#                                self.settings,
#                                self.hotspot_result.protein,
#                                seed=None)
#
#         self.extracted_hotspots = [e]
#         self._rank_extracted_hotspots()
#
#         if self.settings.pharmacophore:
#             self._get_pharmacophores()
#
#         return self.extracted_hotspots
#
#     def extract_all_volumes(self, volume="125", pharmacophores=True):
#         """
#         from the main Fragment Hotspot Map result, the best continuous volume is calculated using peaks in the apolar
#         maps as a seed point.
#
#         :param float volume: volume in Angstrom^3
#         :param bool pharmacophores: if True, generates pharmacophores
#         :return: a `hotspots.result.Results` instance
#
#         >>> result
#         <hotspots.result.Results object at 0x000000001B657940>
#
#         >>> from hotspots.result import Extractor
#
#         >>> extractor = Extractor(result)
#         >>> all_vols = extractor.extract_all_volumes(volume=150)
#         [<hotspots.result.Results object at 0x000000002963A438>,
#          <hotspots.result.Results object at 0x0000000029655240>,
#          <hotspots.result.Results object at 0x000000002963D2B0>,
#          <hotspots.result.Results object at 0x000000002964FDD8>,
#          <hotspots.result.Results object at 0x0000000029651D68>,
#          <hotspots.result.Results object at 0x00000000296387F0>]
#
#         """
#         self.settings.volume = volume
#         self.settings.pharmacophore = pharmacophores
#         self.settings.mode = "seed"
#         self.out_dir = None
#         self.extracted_hotspots = []
#         self._peaks = self._get_peaks()
#
#         for peak in self.peaks:
#             print(peak)
#
#             e = self._from_hotspot(self.single_grid,
#                                    self.masked_dic,
#                                    self.settings,
#                                    self.hotspot_result.protein,
#                                    seed=peak)
#
#             self.extracted_hotspots.append(e)
#
#         self._rank_extracted_hotspots()
#         # generate pharmacophores
#         if self.settings.pharmacophore:
#             self._get_pharmacophores()
#
#         return self.extracted_hotspots
#
#     def _from_hotspot(self, single_grid, mask_dic, settings, prot, seed=None):
#         """
#         create a continuous volume
#
#         :param single_grid:
#         :param mask_dic:
#         :param settings:
#         :param prot:
#         :param seed:
#         :return:
#         """
#         if seed:
#             sphere = single_grid.copy_and_clear()
#             sphere.set_sphere(point=seed, radius=settings._search_radius, value=1, scaling='None')
#             mask = (sphere & single_grid) * single_grid
#         else:
#             mask = single_grid
#
#         optimiser = Extractor._Optimiser(mask=mask, settings=settings, peak=seed)
#         threshold, best_island = optimiser.optimize_island_threshold()
#         print("???????", best_island.extrema[1])
#         if best_island is not None:
#             location, features = self._get_interaction_type(mask_dic, best_island, threshold, settings)
#             grd_dict = self._get_grid_dict(location, features, settings)
#
#             hr = Results(super_grids=grd_dict, protein=prot)
#
#             hr.threshold = threshold
#             hr.best_island = best_island.minimal()
#             hr.location = location
#
#             hr.features = threshold
#             hr.score_value = hr.score()
#
#             hr.rank = hr._rank_features()
#             return hr
#
#     # def _grow_from_seed(self, single_grid, mask_dic, settings, prot, seed=None):
#     #     """
#     #     *experimental*
#     #
#     #     create a Extracted Hotspot object from HotspotResult object
#     #
#     #     :param single_grid:
#     #     :param mask_dic:
#     #     :param settings:
#     #     :param prot:
#     #     :param seed:
#     #     :return:
#     #     """
#     #
#     #     inner = single_grid.copy_and_clear()
#     #     inner.set_sphere(point=seed, radius=1, value=20, scaling='None')
#     #     # mask = (sphere & single_grid) * single_grid
#     #
#     #     # optimiser = _Results.Optimiser(mask=mask, settings=settings, peak=seed)
#     #     # threshold, best_island = optimiser.optimize_island_threshold()
#     #     num_gp = inner.count_grid()
#     #     grown = Grid.grow(inner, single_grid)
#     #     while num_gp < settings._num_gp:
#     #
#     #         grown = Grid.grow(inner, single_grid)
#     #         diff = grown > inner
#     #         if diff.count_grid() < 10:
#     #             break
#     #         inner = grown
#     #         num_gp = inner.count_grid()
#     #         print(num_gp, 'out of', settings._num_gp)
#     #
#     #     tmp_best_island = inner * single_grid
#     #     g_vals = tmp_best_island.grid_values()
#     #     g_vals[::-1].sort()
#     #     try:
#     #         threshold = g_vals[settings._num_gp]
#     #     except IndexError:
#     #         threshold = g_vals.min()
#     #
#     #     best_island = grown
#     #
#     #     if best_island is not None:
#     #         location, features = self._get_interaction_type(mask_dic, best_island, threshold, settings)
#     #         grd_dict = self._get_grid_dict(location, features, settings)
#     #
#     #         hr = Results(super_grids=grd_dict, protein=prot)
#     #
#     #         hr.threshold = threshold
#     #         hr.best_island = best_island.minimal()
#     #         hr.location = location
#     #         hr.features = features
#     #         hr.score_value = hr.score()
#     #         hr.rank = hr._rank_features()
#     #         return hr
#
#     def _get_interaction_type(self, mask_dic, best_island, threshold, settings):
#         """
#         seperates single grid into grid by interaction type
#         :return:
#         """
#         common_best_island = mask_dic["apolar"].common_boundaries(best_island)
#         features_in_vol = {p: g * (g & common_best_island) for p, g in mask_dic.items()}
#         location = features_in_vol["apolar"]
#         features = self.hotspot_result._get_features(interaction_dict=features_in_vol,
#                                                      threshold=threshold,
#                                                      min_feature_gp=settings.min_feature_gp)
#         return location, features
#
#     @staticmethod
#     def _get_grid_dict(location, features, settings):
#         """
#         Creates super grid dict from location and _features
#         :param location:
#         :param features:
#         :param settings:
#         :return:
#         """
#         grid_dic = {"apolar": location.minimal()}
#         interaction_types = set([feat.feature_type for feat in features])
#         feature_by_score = {f.score_value: f for f in features}
#         features = [feature_by_score[s]
#                     for s in sorted([f[0] for f in feature_by_score.items()], reverse=True)][:settings.max_features - 1]
#         for probe in interaction_types:
#             if settings.mode == "seed":
#                 grids = [feat.grid for feat in features
#                          if feat.feature_type == probe and
#                          feat.score_value >= settings.cutoff]
#
#             else:
#                 grids = [feat.grid for feat in features if feat.feature_type == probe]
#             if len(grids) == 0:
#                 grids = [location.minimal().copy_and_clear()]
#
#             grid_dic.update({probe: Grid.super_grid(1, *grids)})
#
#         return grid_dic
#
#     # def get_superstar_result(self, superstar_results):
#     #     """
#     #     finds the overlap between the extracted hotspot and the superstar results
#     #     :param superstar_result:
#     #     :return:
#     #     """
#     #     # TO DO: ALLOW SS RUN IN EXTRACTED_HOTSPOT CLASS
#     #
#     #     extracted_superstar = []
#     #
#     #     for result in superstar_results:
#     #         common_best_island, common_result_grid = Grid.common_grid(self.best_island, result.grid)
#     #         ss_boundary = (common_best_island & common_result_grid) * common_result_grid
#     #         new = copy.copy(result)
#     #         if len(ss_boundary.islands(threshold=2)) != 0:
#     #             g = Grid.super_grid(2, *ss_boundary.islands(threshold=2))
#     #             threshold = g.grid_score(threshold=1, percentile=50)
#     #             print threshold
#     #             new.grid = (g > threshold) * g
#     #         else:
#     #             new.grid = ss_boundary.copy_and_clear()
#     #
#     #         extracted_superstar.append(new)
#     #
#     #     return extracted_superstar
#     #
#     # def calc_feature_profile(self):
#     #     """
#     #     for each hotspot feature, the overlap between the feature sphere and superstar result is calculated
#     #     this is stored as an HotspotFeature attribute (superstar profile)
#     #     :return:
#     #     """
#     #     for feat in self._features:
#     #         super_profile = []
#     #
#     #         for result in self.superstar_results:
#     #             common_result_grid, common_sphere = Grid.common_grid(result.grid, feat.sphere)
#     #             super_sphere = (common_sphere & common_result_grid) * common_result_grid
#     #
#     #             if len(super_sphere.islands(threshold=2)) != 0:
#     #                 result.grid = Grid.super_grid(2, *super_sphere.islands(threshold=2))
#     #
#     #             else:
#     #                 result.grid = feat.sphere.copy_and_clear()
#     #
#     #             super_profile.append(result)
#     #
#     #         feat.superstar_profile = super_profile
#
#     @property
#     def single_grid(self):
#         return self._single_grid
#
#     @property
#     def masked_dic(self):
#         return self._masked_dic
#
#     @property
#     def peaks(self):
#         return self._peaks
#
#     # def grid_post_process(self, super_grids):
#     #     """
#     #     carry out post-processing of fragment hotspot maps
#     #
#     #     Limit the size of polar islands. Keep top scores upto X grid points
#     #     :return:
#     #     """
#     #     for probe, g in super_grids.items():
#     #         if probe == "apolar":
#     #             super_grids.update({probe: g.max_value_of_neighbours()})
#     #
#     #         else:
#     #             h = g.max_value_of_neighbours()
#     #             h = h.limit_island_size(self.settings.island_max_size)
#     #             if h.bounding_box != super_grids["apolar"].bounding_box:
#     #                 h = super_grids["apolar"].common_boundaries(g)
#     #
#     #             super_grids.update({probe: h})
#     #
#     #     # try:
#     #     #     super_grids["negative"] = super_grids["negative"].deduplicate(super_grids["acceptor"],
#     #     #                                                                         threshold=10,
#     #     #                                                                         tolerance=2)
#     #     #
#     #     #     super_grids["positive"] = super_grids["positive"].deduplicate(super_grids["donor"],
#     #     #                                                                         threshold=10,
#     #     #                                                                         tolerance=2)
#     #     # except KeyError:
#     #     #     pass
#     #
#     #     return super_grids
#
#     def _get_peaks(self):
#         """
#         find peak coordinates in apolar maps, used as seeds to find top volumes
#         :return:
#         """
#         apolar = self.hotspot_result.super_grids["apolar"]
#         peaks = feature.peak_local_max(apolar.get_array(),
#                                        min_distance=self.settings.min_distance,
#                                        threshold_abs=self.settings.cutoff)
#         peak_by_value = {}
#         for peak in peaks:
#             val = apolar.value(int(peak[0]), int(peak[1]), int(peak[2]))
#             if val > self.settings.cutoff:
#                 if val in peak_by_value:
#                     peak_by_value[val].append((peak[0], peak[1], peak[2]))
#                 else:
#                     peak_by_value.update({val: [(peak[0], peak[1], peak[2])]})
#
#         average_peaks = []
#         for key in peak_by_value.keys():
#             x = [point[0] for point in peak_by_value[key]]
#             y = [point[1] for point in peak_by_value[key]]
#             z = [point[2] for point in peak_by_value[key]]
#             average_peaks.append(apolar.indices_to_point(int(sum(x) / len(x)),
#                                                          int(sum(y) / len(y)),
#                                                          int(sum(z) / len(z))
#
#                                                          )
#                                  )
#         return average_peaks
#
#     def _get_extracted_hotspots(self):
#         """
#         locate peaks in apolar maps and define fragment size volume
#         :return: list of peak coordinates
#         """
#         extracted_hotspots = []
#         if self.settings.mode == "seed":
#             print(self.peaks)
#             for peak in self.peaks:
#                 print(peak)
#
#                 e = self._from_hotspot(self.single_grid,
#                                        self.masked_dic,
#                                        self.settings,
#                                        self.hotspot_result.protein,
#                                        seed=peak)
#
#                 # if e:
#                 #     if e.threshold > 0:
#                 print(e.threshold)
#                 extracted_hotspots.append(e)
#
#         elif self.settings.mode == "grow":
#             print(self.peaks)
#             for peak in self.peaks:
#                 e = self._grow_from_seed(self.single_grid,
#                                          self.masked_dic,
#                                          self.settings,
#                                          self.hotspot_result.protein,
#                                          seed=peak)
#
#                 # if e:
#                 #     if e.threshold > 0:
#                 print(e.threshold)
#                 extracted_hotspots.append(e)
#
#
#         else:
#
#             e = self._from_hotspot(self.single_grid,
#                                    self.masked_dic,
#                                    self.settings,
#                                    self.hotspot_result.protein,
#                                    seed=None)
#
#             extracted_hotspots.append(e)
#
#         return extracted_hotspots
#
#     def _rank_extracted_hotspots(self):
#         """
#         assigns rank based upon extracted hotspot score
#         :return:
#         """
#         hotspot_by_score = {hotspot.score_value: hotspot for hotspot in self.extracted_hotspots}
#         score = sorted([f[0] for f in hotspot_by_score.items()], reverse=True)
#
#         for i, key in enumerate(score):
#             hotspot_by_score[key].rank = int(i + 1)
#
#         extracted_hotspots_by_rank = {h.rank: h for h in self.extracted_hotspots}
#         self.extracted_hotspots = [value for (key, value) in sorted(extracted_hotspots_by_rank.items())]
#
#         for i, hs in enumerate(self.extracted_hotspots):
#             hs.identifier = hs.score_value  # "rank_{}".format(hs.rank)
#             print("rank", hs.rank, "score", hs.score_value)
#
#     def _select_cavity_grids(self, cavs):
#         """get empty cavity grids"""
#         grds = [Grid(origin=cav.bounding_box[0],
#                      far_corner=cav.bounding_box[1],
#                      spacing=self.settings.spacing,
#                      default=0)
#                 for cav in cavs]
#
#         if self.settings.mode == "seed":
#             filtered = set([g for seed in [p for p in self.peaks]
#                             for g in grds
#                             if g.contains_point(seed)])
#
#         else:
#             raise IOError("Currently only available in seed mode")
#
#         return filtered
#
#     def _get_pharmacophores(self):
#         """
#         generates a pharmacophore model, stores as attribute of hotspot result
#         :return:
#         """
#         for i, hotspot in enumerate(self.extracted_hotspots):
#             hotspot.pharmacophore = hotspot.get_pharmacophore_model(identifier=hotspot.identifier)
#
#     def _write(self, out_dir, mode="best_islands"):
#         """
#         write out information to aid debugging: valid modes:
#             -peaks:
#             -locations: spheres and islands at apolar peak locations
#             -_features: islands and probes at feature point locations
#         """
#
#         pymol_out = pymol_imports()
#         if mode == "peaks":
#             out_dir = Helper.get_out_dir(join(out_dir))
#             for i, peak in enumerate(self.peaks):
#                 score = "{0:.2f}".format(self.hotspot_result.super_grids["apolar"].value_at_point(peak))
#                 sphere = 'score_{0} = [COLOR, 1.00, 1.000, 0.000] + ' \
#                          '[ALPHA, 0.8] + ' \
#                          '[SPHERE, float({1}), float({2}), float({3}), float(0.5)]\n' \
#                     .format(i, peak[0], peak[1], peak[2])
#                 pymol_out += sphere
#                 pymol_out += '\ncmd.load_cgo(score_{1}, "score_{0}", 1)'.format(score, i)
#                 pymol_out += '\ncmd.group("Peaks", members="score_{0}")\n'.format(score)
#             with open(join(out_dir, "peaks.py"), "w") as pymol_file:
#                 pymol_file.write(pymol_out)
#
#         elif mode == "best_islands":
#             out_dir = Helper.get_out_dir(join(out_dir, "best_islands"))
#             thresholds = []
#             for i, extracted in enumerate(self.extracted_hotspots):
#                 extracted.best_island.write(join(out_dir, "island_{}.grd".format(i)))
#                 thresholds.append(extracted.threshold)
#             pymol_out += """
# nh = {0}
# thresholds = {1}
# for n in range(nh):
#     cmd.load(r'best_islands/island_%s.grd' % (n), 'apolar_%s' % (n))
#     cmd.isosurface('surface_apolar_%s' % (n), 'apolar_%s' % (n), thresholds[n])
#     cmd.set('transparency', 0.7, 'surface_apolar_%s' % (n))
#     cmd.color('yellow', 'surface_apolar_%s' % (n))
# for n in range(nh):
#     cmd.group('hotspot_%s'%(n), members= 'surface_apolar_%s'%(n))
#     cmd.group('hotspot_%s'%(n), members= 'apolar_%s'%(n))""" \
#                 .format(len(self.extracted_hotspots), thresholds)
#
#             with open(join(dirname(out_dir), "best_islands.py"), "w") as pymol_file:
#                 pymol_file.write(pymol_out)
#
#         else:
#             raise IOError("mode not supported")
