args <- commandArgs(trailingOnly = TRUE)

args.label <- args[1]
args.input <- args[2]
args.output <- args[3]

library(dplyr)
library(tidyr)

ibarfc <- read.table(
    args.input,
    header=TRUE, sep='\t',
    stringsAsFactors=FALSE
)

colnames(ibarfc) <- c('REAGENT_ID', 'GENE', args.label)

ibarfc$REAGENT_ID <- gsub(
    '_gene', '-gene', ibarfc$REAGENT_ID
)

ibarfc$REAGENT_ID <- gsub(
    'Control_', 'Control-', ibarfc$REAGENT_ID
)

ibarfc$GENE <- gsub(
    '_gene', '-gene', ibarfc$GENE
)

ibarfc$GENE <- gsub(
    'Control_', 'Control-', ibarfc$GENE
)


ibarfc$Guide <- unlist(
    lapply(
        strsplit(
            ibarfc$REAGENT_ID,
            split=':'
        ),
        '[',
        2
    )
)

ibarfc$Barcode <- unlist(
    lapply(
        strsplit(
            ibarfc$REAGENT_ID,
            split=':'
        ),
        '[',
        3
    )
)

ibarfc$REAGENT_ID <- paste(
    ibarfc$GENE, ibarfc$Guide, sep=':'
)

ibarfc$Guide <- NULL


ibarfc <- ibarfc%>% group_by(REAGENT_ID) %>% mutate(id=row_number(REAGENT_ID))

ibarfc$Barcode <- NULL

repfc <- eval(
    substitute(
        spread(d, 'id', a),
        list(d=ibarfc,a=args.label)
    )
)

colnames(repfc) <- c(
    'REAGENT_ID', 'GENE', paste(args.label, seq(4), sep='.')
)

nalist <- as.list(rep(0, 4))

names(nalist) <- paste(args.label, seq(4), sep='.')
repfc <- replace_na(repfc, nalist)[
    c('REAGENT_ID', 'GENE', paste(args.label, seq(4), sep='.'))
]

write.table(
    repfc,
    args.output,
    row.names=FALSE, quote=FALSE, sep='\t'
)
