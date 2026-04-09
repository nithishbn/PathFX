import pickle,os,csv

### 1. Enter the name of the analysis
### 2. Enter the name of the drug. This can be a drug bank ID or the name of drug that isn't in drug bank
### 3. (Optional) Enter the targets. This is only optional if you have entered a drug from DrugBank. For investigational drugs, you will need to enter the targets associated with that drug.


analysis_name = 'PathFX_FlurNap'
drug_name = 'Naftifine'

### call the algorithm without phenotype clustering
cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s'%(drug_name,analysis_name)
print(cmd)
os.system(cmd)

drug_name = 'Vedolizumab'
cmd = 'python phenotype_enrichment_pathway.py -d %s -a %s'%(drug_name,analysis_name)
os.system(cmd)

