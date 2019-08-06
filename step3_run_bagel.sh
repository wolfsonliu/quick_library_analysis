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
#SBATCH --time=120:00:00

set -e -u -o pipefail

module unload anaconda/3-4.4.0.1
module load anaconda/2-4.4.0.1

DIR_BAGEL=$(pwd)/BAGEL 

DIR_BASE=$(pwd)
DIR_DATA=${DIR_BASE}/data
DIR_RAW_COUNT=${DIR_DATA}/count
DIR_DATA_IBAR_COUNT=${DIR_DATA}/bagel/ibar_count
DIR_DATA_IBAR_FC=${DIR_DATA}/bagel/ibar_fc
DIR_DATA_IBAR_BF=${DIR_DATA}/bagel/ibar_bf
DIR_DATA_REP_FC=${DIR_DATA}/bagel/replicate_fc
DIR_DATA_REP_BF=${DIR_DATA}/bagel/replicate_bf


LABELs=(
K562
HeLa
)


REFs=(
K562
HeLa
)

lab=${LABELs[$SLURM_ARRAY_TASK_ID]}
ref=${REFs[$SLURM_ARRAY_TASK_ID]}

echo -e "gene\t${ref}\t${lab}" > ${DIR_DATA_IBAR_COUNT}/${lab}.txt

awk -F "\t" 'FNR == NR && FNR != 1 && $1 != "NA" {
    idx = ($1 ":" $2 ":" $3);
    a[idx] = $4;
    g[idx] = $1;
}
FNR != NR && FNR != 1 && $1 != "NA" {
    idx = ($1 ":" $2 ":" $3);
    b[idx] = $4;
}
END {
    for (x in a) {
        print x FS g[x] FS a[x] FS b[x];
    }
}' ${DIR_RAW_COUNT}/${ref}.count ${DIR_RAW_COUNT}/${lab}.count | \
sort -k 1,1 >> ${DIR_DATA_IBAR_COUNT}/${lab}.txt

python2 ${DIR_BAGEL}/BAGEL-calc_foldchange.py \
    -i ${DIR_DATA_IBAR_COUNT}/${lab}.txt \
    -o ${DIR_DATA_IBAR_FC}/${lab} -c 1

Rscript ${DIR_BASE}/bagel_prepare_replicate_fc.r ${lab} \
    ${DIR_DATA_IBAR_FC}/${lab}.foldchange \
    ${DIR_DATA_REP_FC}/${lab}.foldchange

python2 ${DIR_BASE}/BAGEL_rm0.py -i ${DIR_DATA_REP_FC}/${lab}.foldchange \
    -o ${DIR_DATA_REP_BF}/${lab}.txt \
    -e ${DIR_BAGEL}/training_essentials.txt \
    -n ${DIR_BAGEL}/training_nonessential.txt -c 1,2,3,4

module unload anaconda/2-4.4.0.1
