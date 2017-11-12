# Make figure of nitrogen marker genes by lake
library(ggplot2)
library(cowplot)

genes <- c("nitrogenase", "nitrate reductase", "nitrite reductase", "nitric oxide reductase", "arginase", "glutamine synthetase", "put/sperm transport")
me.counts <- c(27, 19, 53, 7, 65, 149, 168)/101 * 100
tb.counts <- c(128, 69, 80, 43, 24, 141, 151)/101 * 100

# Counts obtained by searching for key genes in IMG function search. Can have more than one hit per genome. 

ME.plot <- data.frame(genes, me.counts)
TB.plot <- data.frame(genes, tb.counts)

ggplot(data = ME.plot, aes(x = genes, y = me.counts)) + geom_bar(stat = "identity", fill = "#5ab4ac", color = "black") + labs(x = "", y = "") + coord_flip() + scale_y_reverse(limits = c(170, 0))

ggplot(data = TB.plot, aes(x = genes, y = tb.counts)) + geom_bar(stat = "identity", fill = "#d8b365", color = "black") + labs(x = "", y = "") + coord_flip() + ylim(0, 170)
