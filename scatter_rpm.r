#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)


data <- list()
data.norm <- list()

for (x in seq(2)) {
    label <- unlist(strsplit(basename(args[x]), '[.]'))[1]
    data[[x]] <- read.table(args[x], header=TRUE, stringsAsFactors=FALSE)
    colnames(data[[x]]) <- c('gene', 'guide', 'barcode', 'count')
    data[[x]][[paste0('R',x)]] <- data[[x]][['count']] / sum(data[[x]][['count']], na.rm=TRUE) * 10^6
    data[[x]]$count <- NULL
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

alldata <- c(pdata$R1, pdata$R2)
alldata <- na.omit(alldata[is.finite(alldata)])

lim1 <- min(log10(alldata + 1), na.omit=TRUE)

lim2 <- max(log10(alldata + 1), na.omit=TRUE)

p <- ggplot(
    pdata, aes(x=log10(R1 + 1), y=log10(R2 + 1), color=type)
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
                            paste(basename(args[1]), '(log10(RPM + 1))', sep=' ')
                        ) + ylab(
                                paste(basename(args[2]), '(log10(RPM + 1))', sep=' ')
                            ) + ggtitle(
                                    paste('Corr:', round(cor(pdata$R1, pdata$R2), 2))
                                ) + coord_fixed(ratio=1) 


ggsave(args[length(args)], p, units='cm')
