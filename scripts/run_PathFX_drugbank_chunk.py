import pickle, os, sys

chunk_id = int(sys.argv[1])
total_chunks = int(sys.argv[2])

rscs_dir = '../rscs/'
db2n = pickle.load(open(os.path.join(rscs_dir, 'pfxDB050620_dbid2name.pkl'), 'rb'))
dint = pickle.load(open(os.path.join(rscs_dir, 'pfxDB050620_dint.pkl'), 'rb'))

analysis_name = 'all_network_results'
res_dir = os.path.join('../results/', analysis_name)
if not os.path.exists(res_dir):
    os.makedirs(res_dir)

drugs = sorted([d for d in db2n if d in dint])
chunk_drugs = [d for i, d in enumerate(drugs) if i % total_chunks == (chunk_id - 1)]

print('Chunk %d/%d: %d drugs' % (chunk_id, total_chunks, len(chunk_drugs)))

for dbid in chunk_drugs:
    print(dbid)
    cmd = 'python phenotype_enrichment_pathway_Pfx050120.py -d %s -a %s' % (dbid, analysis_name)
    os.system(cmd)
