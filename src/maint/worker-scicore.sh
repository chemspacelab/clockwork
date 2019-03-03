#!/bin/bash
#SBATCH --job-name=exconwork
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --error=/dev/null
#SBATCH --output=/dev/null
#SBATCH --time=00:30:00
#SBATCH --time-min=00:05:00
#SBATCH --qos=30min

module load intel/2019.01
module load Anaconda3/5.0.1

export OMP_NUM_THREADS=1

source activate /scicore/home/lilienfeld/rudorff/opt/conda/rdkit
connection=$(cat /scicore/home/lilienfeld/rudorff/.redis-credentials)
export PYTHONPATH=/scicore/home/lilienfeld/rudorff/.local/lib/python3.7/site-packages:$PYTHONPATH

python3 src/worker/worker.py \
        -f /scicore/home/lilienfeld/bexusi36/qm9-C7O2H10/full.sdf.gz \
        --torsions-file /scicore/home/lilienfeld/bexusi36/qm9-C7O2H10/list_torsions_idx \
        --connect-redis "$connection"

