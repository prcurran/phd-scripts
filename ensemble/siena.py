"""

This script wraps around the Protein Plus web services
    - SIENA: automated construction of protein ensembles from the PDB

    - Bietz, S. Rarey, M.: SIENA: Efficient Compilation of Selective Protein Binding Site Ensembles.
    Journal of Chemical Information and Modeling,56(1): 248-59.

    - Bietz, S. Rarey, M.: ASCONA: Rapid Detection and Alignment of Protein Binding Site Conformations.
    Journal of Chemical Information and Modeling, 55(8):1747-1756.

"""
import shlex, subprocess, ast, time, json, os, csv, re

import urllib2
import urlparse
from collections import OrderedDict


class Ensemble(object):
    """

    """
    def __init__(self, data):
        self.data = data
        self.to_retrieve = ['result_table', 'pdb_files', 'ligands', 'alignment']

    def _write(self, url, out_dir):
        """

        :param url:
        :param out_dir:
        :return:
        """
        print url
        url = re.sub('/esults/', '/results/', url)              # fix url
        a = urlparse.urlparse(url)
        fname = os.path.basename(a.path)
        ext = fname.split(".")[1]
        path = os.path.join(out_dir, fname)

        with open(path, "wb") as w:
            if ext == "csv":
                # format csv
                stri = [[a for a in line.split(";") if a != '\n'] for line in urllib2.urlopen(url).readlines()]
                csv_writer = csv.writer(w, delimiter=',')
                for s in stri:
                    csv_writer.writerow(s)
            else:
                # otherwise no formatting requried
                w.write(urllib2.urlopen(url).read())

    def save(self, out_dir):
        """
        save the ensemble data to output directory

        :param str out_dir: path to output directory
        """
        for t in self.to_retrieve:
            new_out_dir = os.path.join(out_dir, t)
            if not os.path.exists(new_out_dir):
                os.mkdir(new_out_dir)
            urls = self.data[t]
            if type(urls) is list:
                for url in urls:
                    self._write(url, new_out_dir)
            else:
                self._write(urls, new_out_dir)


class Search(object):
    """
    two parts


    :param pdb_code:
    :param ligand:
    :param mode:
    :param data_reduction:
    :param settings:
    """
    class Settings(object):
        """


        """
        def __init__(self):
            self.url = 'https://proteins.plus/api/siena_rest'
            self.data = {"reduction_procedure":"",
                         "siena": {"bb_clustering":"",
                                   "all_atom_clustering":"",
                                   "ligand_driven_selection":"",
                                   "ligand": "",
                                   "pocket":"",
                                   "pdbCode": "",
                                   "siteRadius":"6.5",
                                   "fragment_length": "10",
                                   "flexibility_sensitivity": "0.6",
                                   "fragment_distance": "4",
                                   "minimalSiteIdentity": "0.7",
                                   "minimalSiteCoverage":"",
                                   "maximum_mutations":"",
                                   "resolution":"",
                                   "maximumBackbone":"",
                                   "depositionYear":"",
                                   "ecNumber":"",
                                   "electronDensityAvailable": "",
                                   "identicalGlobalSequence": "",
                                   "noGlobalMutations": "",
                                   "unique_sequence":"",
                                   "holo_only": "",
                                   "unique_ligands": "",
                                   "complete_residues_only": ""
                                   }
                         }


    def __init__(self, settings=None):
        """

        :param pdb_code:
        :param ligand:
        :param mode:
        :param data_reduction:
        :param settings:
        """
        if settings is None:
            self.settings = self.Settings()
        else:
            self.settings = settings

    def create_ensemble(self, pdb_code, ligand, reduction_procedure, num_members):
        """


        :param pdb_code:
        :param ligand:
        :param reduction_procedure:
        :return:
        """
        self.settings.data['siena']["pdbCode"] = pdb_code
        self.settings.data['siena']["ligand"] = ligand
        self.settings.data['siena']["ligand"] = ligand
        self.settings.data['reduction_procedure'] = reduction_procedure
        self.settings.data['siena'][reduction_procedure] = num_members  # reduction_procedure

        self.results_url = self._post()
        self._get()
        self.results_url = self._post()
        data = self._get()

        return Ensemble(data)


    def _run(self, cmd):
        """
        runs commandline procedure (urllib doesn't work for some reason)

        :param str cmd: command line str using curl
        :return:
        """
        args = shlex.split(cmd)
        return subprocess.check_output(args)

    def _post(self):
        """
        Initiate the job using a POST request

        :return:
        """
        cmd = """curl -d '{}' -H "Accept: application/json" -H "Content-Type: application/json" -X POST {}"""\
            .format(json.dumps(self.settings.data, ensure_ascii=True, indent=2, default=True),self.settings.url)

        response = ast.literal_eval(self._run(cmd))
        return response['location']

    def _get(self):
        """
        Collect the results using GET. May have to try multiple times.

        :return:
        """
        cmd = """curl {}""".format(self.results_url)

        response = ast.literal_eval(self._run(cmd))
        if response["status_code"] == 400:
            print "status code: {}".format(response['status_code'])
            raise RuntimeError()

        while response["status_code"] == 220:
            time.sleep(10)
            response = ast.literal_eval(self._run(cmd))

        # create Ensemble
        print "status code: {}".format(response['status_code'])
        return response


def main():
    searcher = Search()
    ensemble = searcher.create_ensemble(pdb_code="1ia1",
                                        ligand="TQ3_A_194",
                                        reduction_procedure="all_atom_clustering",
                                        num_members="5")

    ensemble.save(out_dir = "/home/pcurran/patel/CDK2/1aq1_ensemble")


if __name__ == "__main__":
    main()

