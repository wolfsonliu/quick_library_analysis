#! /bin/bash

LABELs=(
)

for x in ${LABELs[*]}; do
    Rscript scatter_zlfc.r \
        data/zfc/${x}_R1/${x}_R1_gene.txt \
        data/zfc/${x}_R2/${x}_R2_gene.txt \
        fig/scatter_zlfc_${x}.pdf
    convert -density 300 fig/scatter_zlfc_${x}.pdf \
        fig/scatter_zlfc_${x}.png
done


ONE=(
)

TWO=(
)

for i in `seq 0 9`; do
    Rscript scatter_zlfc.r \
        data/zfc/${ONE[${i}]}/${ONE[${i}]}_gene.txt \
        data/zfc/${TWO[${i}]}/${TWO[${i}]}_gene.txt \
        fig/scatter_zlfc_${ONE[${i}]}-${TWO[${i}]}.pdf
    convert -density 300 fig/scatter_zlfc_${ONE[${i}]}-${TWO[${i}]}.pdf \
        fig/scatter_zlfc_${ONE[${i}]}-${TWO[${i}]}.png
done
