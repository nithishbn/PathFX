#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=12G,h_rt=24:00:00
#$ -pe shared 1
#$ -o logs/pathfx_all_drugbank.log
#$ -M nithishbn@g.ucla.edu
#$ -m bea

# PathFX - Run enrichment analysis for all DrugBank drugs
# Usage: qsub submit_pathfx_all_drugbank.sh
# Results saved to results/all_network_results/

mkdir -p logs

. /u/local/Modules/default/init/modules.sh
module load mamba/1.5.5

mamba activate pathfx

echo "=========================================="
echo "PathFX All DrugBank Analysis"
echo "Job ID: $JOB_ID"
echo "Hostname: $HOSTNAME"
echo "Date: $(date)"
echo "=========================================="

cd scripts/

/usr/bin/time -v python run_PathFX_all_drugBank.py

EXIT_CODE=$?

echo "=========================================="
echo "Finished PathFX All DrugBank"
echo "Exit code: $EXIT_CODE"
echo "Date: $(date)"
echo "=========================================="

exit $EXIT_CODE
