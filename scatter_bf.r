#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)

data <- list()

for (x in seq(2)) {
    label <- paste0('R', x)
    data[[x]] <- read.table(args[x], header=TRUE, stringsAsFactors=FALSE)[c('GENE', 'BF')]
    colnames(data[[x]]) <- c('gene', label)
    data[[x]]$type <- 'Genes'
    data[[x]][grepl('^aavs1', data[[x]]$gene), 'type'] <- 'AAVS1'
    data[[x]][grepl('^negCon|^negative', data[[x]]$gene), 'type'] <- 'Non-targeting Control'
}

pdata <- merge(data[[1]], data[[2]])

pdata$type <- factor(
    pdata$type,
    levels=c('Genes', 'AAVS1', 'Non-targeting Control')
)

pdata <- pdata[order(pdata$type),]

## plot figure


alldata <- c(pdata$R1, pdata$R2)
alldata <- na.omit(alldata[is.finite(alldata)])

lim1 <- min(alldata, na.omit=TRUE)

lim2 <- max(alldata, na.omit=TRUE)

p <- ggplot(
    pdata, aes(x=R1, y=R2, color=type)
) + geom_point(
    ) + scale_color_manual(
            values=c(
                'Non-targeting Control'='black',
                'AAVS1'='red',
                'Genes'='lightgray'
            )
        ) + theme(
                aspect.ratio=1,
                panel.background=element_rect(fill='transparent', color='black'),
                panel.grid.major.x=element_line(
                    color="#d9d9d9", linetype='dashed'
                ),
                panel.grid.minor.x=element_blank(),
                panel.grid.major.y=element_line(
                    color="#d9d9d9", linetype='dashed'
                ),
                panel.grid.minor.y=element_blank()
            ) + xlim(
                    lim1, lim2
                ) + ylim(
                    lim1, lim2
                    ) + xlab(
                        paste(args[1], '(Bayesian Factor)', sep=' ')
                    ) + ylab(
                            paste(args[2], '(Bayesian Factor)', sep=' ')
                        ) + ggtitle(
                                paste('Corr:', round(cor(pdata$R1, pdata$R2), 2))
                            ) + coord_fixed(ratio=1)

ggsave(args[length(args)], p, units='cm')
