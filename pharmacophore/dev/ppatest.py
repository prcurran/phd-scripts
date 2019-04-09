from pdb_python_api import Query, PDB


q = Query()
q.add(query_type="UpAccessionIdQuery", query_parameters={'accessionIdList':'P24941'})
q.add(query_type="NoLigandQuery", query_parameters={'haveLigands':'yes'})
q.add(query_type="ResolutionQuery", query_parameters={'refine.ls_d_res_high.comparator': 'between',
                                                      'refine.ls_d_res_high.min': '0',
                                                      'refine.ls_d_res_high.max': '1.5'})

p = PDB()
print p.search()



#
# p = PDB()
# print len(p.search(q.query))