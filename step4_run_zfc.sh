#! /bin/bash
#SBATCH --job-name=
#SBATCH --output=arrayjob.%x.%j_%A_%a.%N.out
#SBATCH --error=arrayjob.%x.%j_%A_%a.%N.err
#SBATCH --partition=C032M0128G
#SBATCH --qos=low
#SBATCH --get-user-env
#SBATCH -n 1
#SBATCH --cpu-freq=high
#SBATCH -A hpc
#SBATCH --mail-type=end
#SBATCH --mail-user=
#SBATCH --time=0:24:00

set -e -u -o pipefail

# source activate zfc


DIR_BASE=$(pwd)
DIR_DATA=${DIR_BASE}/data
DIR_RAW_COUNT=${DIR_DATA}/count
DIR_ZFC=${DIR_DATA}/zfc

LABELs=(
K562
)


REFs=(
HeLa
)

lab=${LABELs[$SLURM_ARRAY_TASK_ID]}
ref=${REFs[$SLURM_ARRAY_TASK_ID]}

mkdir -p ${DIR_ZFC}/${lab}

echo -e "gene\tguide\tbarcode\tctrl\texp" > ${DIR_ZFC}/${lab}/${lab}.txt

awk -F "\t" 'FNR == NR && FNR != 1 && $1 != "NA" {
    idx = ($1 "\t" $2 "\t" $3);
    a[idx] = 0 + $4;
    g[idx] = $1;
}
FNR != NR && FNR != 1 && $1 != "NA" {
    idx = ($1 "\t" $2 "\t" $3);
    b[idx] = 0 + $4;
}
END {
    for (x in a) {
        print x FS a[x] FS b[x];
    }
}' ${DIR_RAW_COUNT}/${ref}.count ${DIR_RAW_COUNT}/${lab}.count | \
sort -k 1,1 >> ${DIR_ZFC}/${lab}/${lab}.txt

zfc -i ${DIR_ZFC}/${lab}/${lab}.txt -o ${DIR_ZFC}/${lab}/${lab} --plot

# source deactivate
####################
