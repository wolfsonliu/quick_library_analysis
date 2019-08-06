library(VennDiagram)
library(grid)

labels <- c(
    'HeLa_IBAR_MOI3_D15_Rm',
    'HeLa_IBAR_MOI3_D21_Rm',
    'HeLa_IBAR_MOI3_Exp_Rm',
    'HeLa_MOI3_D15_R1',
    'HeLa_MOI3_D15_R2',
    'HeLa_MOI3_D15_Rm',
    'HeLa_MOI3_D21_R1',
    'HeLa_MOI3_D21_R2',
    'HeLa_MOI3_D21_Rm',
    'HeLa_MOI3_Exp_Rm',
    'K562_MOI10_D30_Rm',
    'K562_MOI3_D30_Rm'
)

data <- list()

for (x in labels) {
    data[[x]] <- read.table(
        file.path(
            'venn', paste(x, 'txt', sep='.')
        ),
        header=FALSE, stringsAsFactors=FALSE
    )$V1
}

# HeLa
data[['2015_Cell']] <- read.table(
    file.path('venn', 'hart.txt'),
    header=FALSE, stringsAsFactors=FALSE
)$V1
# K562
data[['2015_Science']] <- read.table(
    file.path('venn', '2015Science-K562+essential+gene+list.csv'),
    header=FALSE, stringsAsFactors=FALSE
)$V1
data[['2016_NBT']] <- read.table(
    file.path('venn', '2016NBT-K562+essential+gene+list.csv'),
    header=FALSE, stringsAsFactors=FALSE
)$V1

plotpairs <- list(
    c('HeLa_IBAR_MOI3_D15_Rm', 'HeLa_IBAR_MOI3_D21_Rm', 'HeLa_IBAR_MOI3_Exp_Rm'),
    c('HeLa_IBAR_MOI3_D15_Rm', 'HeLa_IBAR_MOI3_D21_Rm'),
    c('HeLa_MOI3_D15_R1', 'HeLa_MOI3_D15_R2'),
    c('HeLa_MOI3_D21_R1', 'HeLa_MOI3_D21_R2'),
    c('HeLa_MOI3_D15_R1', 'HeLa_MOI3_D15_R2', 'HeLa_MOI3_D15_Rm'),
    c('HeLa_MOI3_D21_R1', 'HeLa_MOI3_D21_R2', 'HeLa_MOI3_D21_Rm'),
    c('HeLa_MOI3_D15_Rm', 'HeLa_MOI3_D21_Rm', 'HeLa_MOI3_Exp_Rm'),
    c('HeLa_MOI3_D15_Rm', 'HeLa_IBAR_MOI3_D15_Rm', '2015_Cell'),
    c('HeLa_MOI3_D21_Rm', 'HeLa_IBAR_MOI3_D21_Rm', '2015_Cell'),
    c('HeLa_MOI3_Exp_Rm', 'HeLa_IBAR_MOI3_Exp_Rm', '2015_Cell'),
    c('K562_MOI10_D30_Rm', 'K562_MOI3_D30_Rm'),
    c('2015_Science', '2016_NBT'),
    c('K562_MOI10_D30_Rm', 'K562_MOI3_D30_Rm', '2015_Science', '2016_NBT'),
    c('K562_MOI3_D30_Rm', '2015_Science', '2016_NBT'),
    c('K562_MOI10_D30_Rm', '2015_Science', '2016_NBT')
)

for (x in plotpairs) {
    pdf(
        file.path('venn', paste(paste(x, collapse='-'), 'pdf', sep='.'))
    )
    plotdata <- data[x]
    names(plotdata) <- gsub(
        '(K562_|HeLa_)', '', names(plotdata)
    )
    vplot <- venn.diagram(
        x=plotdata, filename=NULL,
        fill=c("cornflowerblue", "green", "yellow", "darkorchid1", "brown1")[seq(length(x))],
        col="transparent",
        alpha=0.5
    )
    grid.draw(vplot)
    dev.off()
}
