#!/bin/bash
#SBATCH --job-name=exconwork
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --error=/dev/null
#SBATCH --output=/dev/null
#SBATCH --time-min=00:05:00
#SBATCH --time=02:00:00
#SBATCH --partition=normal,classic,long

export OMP_NUM_THREADS=1

. /home/vonrudorff/opt/conda/install/etc/profile.d/conda.sh
conda activate /home/vonrudorff/opt/conda/rdkit
connection=$(cat /home/vonrudorff/.redis-credentials)

python3 /home/vonrudorff/workcopies/clockwork/src/worker/worker.py \
        -f /home/vonrudorff/workcopies/qm9-C7O2H10/full.sdf.gz \
        --torsions-file /home/vonrudorff/workcopies/qm9-C7O2H10/list_torsions_idx \
        --connect-redis "$connection"

