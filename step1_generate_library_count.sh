#! /bin/bash
#SBATCH --job-name=
#SBATCH --output=arrayjob.%x.%j_%A_%a.%N.out
#SBATCH --error=arrayjob.%x.%j_%A_%a.%N.err
#SBATCH --partition=C032M0128G
#SBATCH --qos=low
#SBATCH --get-user-env
#SBATCH -n 2 
#SBATCH --cpu-freq=high
#SBATCH -A hpc
#SBATCH --mail-type=end
#SBATCH --mail-user=
#SBATCH --time=120:00:00

set -u -e -o pipefail

DIR_BASE=$(pwd)
DIR_DATA=${DIR_BASE}/data
DIR_DATA_FQ=${DIR_DATA}/rawdata
DIR_DATA_COUNT=${DIR_DATA}/count

# library file
# gene <tab> guide <tab> barcode
FILE_LIBRARY=${DIR_DATA}/library.txt

#

FQ1s=(
${DIR_DATA_FQ}/K562_1.fq
${DIR_DATA_FQ}/HeLa_1.fq
)



FQ2s=(
${DIR_DATA_FQ}/K562_2.fq
${DIR_DATA_FQ}/HeLa_2.fq
)

LABELS=(
K562
HeLa
)

fq1=${FQ1s[$SLURM_ARRAY_TASK_ID]}
fq2=${FQ2s[$SLURM_ARRAY_TASK_ID]}

fq1gz=${fq1}.gz
fq2gz=${fq2}.gz

lab=${LABELS[$SLURM_ARRAY_TASK_ID]}

# Unzip file

echo -e "${fq1gz}\n${fq2gz}" | xargs -n 1 -P 2 -I {} gunzip {}

# Count

${DIR_BASE}/library_count_sgrna_with_barcode \
    -a "ACCG([ATGC]{20})GTTT[ATGC]{5,35}TGGA([ATCG]{6})AACA" \
    -p "ACCG([ATGC]{20})GTTT" \
    -b "TGGA([ATCG]{6})AACA" \
    -f ${fq1} -r ${fq2} -l ${FILE_LIBRARY} \
> ${DIR_DATA_COUNT}/${lab}.count

# Count fq reads

echo -e "${lab}\tfq_reads\t$(wc -l ${fq1} | awk -F " " '{print $1 / 4.0;}')" >> ${DIR_DATA}/count_stats.txt

# Count library raw count

echo -e "${lab}\tlibrary_counts\t"$(awk -F " " 'FNR != 1 {a += $4;} END {print a;}' ${DIR_DATA_COUNT}/${lab}.count ) \
>> ${DIR_DATA}/count_stats.txt

# Zip file

echo -e "${fq1}\n${fq2}" | xargs -n 1 -P 2 -I {} gzip {}

