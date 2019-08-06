#! /bin/bash

Rscript count_barplot.r data/count_stats.txt fig/count_barplot.pdf

LABELs=(
)

for x in ${LABELs[*]}; do
    Rscript scatter_rpm.r \
        data/count/${x}_R1.count \
        data/count/${x}_R2.count \
        fig/scatter_rpm_${x}.pdf

    convert -density 300 fig/scatter_rpm_${x}.pdf \
        fig/scatter_rpm_${x}.png
done
