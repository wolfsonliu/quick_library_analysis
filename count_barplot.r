#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)

data <- read.table(args[1], header=FALSE, stringsAsFactors=FALSE)
colnames(data) <- c('Label', 'Type', 'Count')
data <- data[order(data$Label),]

data <- data %>% group_by(
                     Label
                 ) %>% mutate(
                           Rate=Count/max(Count)
                       )

data$Percentage <- ifelse(
    data$Rate < 1,
    paste0(round(data$Rate * 100, 2), '%'),
    ''
)

data$Type <- as.factor(data$Type)

typen <- length(levels(data$Type))

p <- ggplot(
    data, aes(x=Label, y=Count, group=Type, fill=Type)
) + geom_bar(
        stat='identity', position='dodge'
    ) + geom_text(
            mapping=aes(
                x=Label, y=Count, label=Percentage
            ),
            nudge_x=(as.numeric(data$Type) - 1) * (0.5 / typen),
            nudge_y=rep(0.05 * max(data$Count), dim(data)[1])
        ) + theme(
                panel.background=element_rect(
                    fill='transparent', color='black'
                ),
                panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                panel.grid.major.y=element_line(
                    color="#d9d9d9", linetype='dashed'
                ),
                panel.grid.minor.y=element_blank()
            ) + xlab('Experiment') + ylab('Counts')

ggsave(
    args[2], p,
    width=4 * length(unique(data$Label)),
    height=10,
    units='cm'
)
