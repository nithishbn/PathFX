#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_data=12G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-64
#$ -o out/pathfx_drugbank_$TASK_ID.log
#$ -M nithishbn@g.ucla.edu
#$ -m bea

# PathFX - All DrugBank drugs, chunked array job
# Usage: qsub submit_pathfx_drugbank_array.sh
# 7012 drugs split across 64 chunks (~110 drugs/chunk, ~9 min each)
# Results saved to results/all_network_results/

mkdir -p out

. /u/local/Modules/default/init/modules.sh
module load mamba/1.5.5

mamba activate pathfx

CHUNK_ID=$SGE_TASK_ID
TOTAL_CHUNKS=64

echo "=========================================="
echo "PathFX DrugBank Array - Chunk $CHUNK_ID of $TOTAL_CHUNKS"
echo "Job ID: $JOB_ID"
echo "Task ID: $SGE_TASK_ID"
echo "Hostname: $HOSTNAME"
echo "Date: $(date)"
echo "=========================================="

cd scripts/

/usr/bin/time -v python run_PathFX_drugbank_chunk.py $CHUNK_ID $TOTAL_CHUNKS

EXIT_CODE=$?

echo "=========================================="
echo "Finished Chunk $CHUNK_ID"
echo "Exit code: $EXIT_CODE"
echo "Date: $(date)"
echo "=========================================="

exit $EXIT_CODE
