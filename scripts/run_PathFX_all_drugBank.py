import pickle, os, csv

rscs_dir = '../rscs/'

db2n = pickle.load(open(os.path.join(rscs_dir, 'drugbankid_to_name.pkl'), 'rb'))
dint = pickle.load(open(os.path.join(rscs_dir, 'drug_intome_targets.pkl'), 'rb'))

analysis_name = 'all_network_results'
res_dir = os.path.join('../results/', analysis_name)
if not os.path.exists(res_dir):
    os.makedirs(res_dir)

for dbid in sorted(db2n.keys()):
    if dbid in dint:
        print(dbid)
        cmd = 'python phenotype_enrichment_pathway_Pfx050120.py -d %s -a %s' % (dbid, analysis_name)
        os.system(cmd)
