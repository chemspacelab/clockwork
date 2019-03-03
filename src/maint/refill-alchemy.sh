#!/bin/bash
# connect to REDIS
. /home/vonrudorff/opt/conda/install/etc/profile.d/conda.sh
conda activate /home/vonrudorff/opt/conda/rdkit
connection=$(cat /home/vonrudorff/.redis-credentials)
CHEMSPACELAB_REDIS_CONNECTION="$connection" python $1/../worker/rediscomm.py --haswork PRODUCTION || exit

# work left, let's submit
CURRENT=$(/home/admin/slurm/slurm16/bin/squeue -u $(whoami) | grep exconwo | wc -l)
TARGET=2
DIFF=$(($TARGET-CURRENT))
[ $DIFF -lt 0 ] && DIFF=0

for i in $(seq 1 $DIFF)
do
	/home/admin/slurm/slurm16/bin/sbatch $1/worker-alchemy.sh &> /dev/null
done
