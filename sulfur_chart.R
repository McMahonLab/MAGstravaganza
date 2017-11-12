# Make plot of # genomes with sulfur genes
library(ggplot2)
library(cowplot)

types <- c("assimilatory sulfate reduction", "dissimilatory sulfate reduction", "sulfide oxidation")
me <- c(45, 0, 6)
tb <- c(45, 3, 18)
me <- me/101 * 100
tb <- tb/102 * 100

ME.plot <- data.frame(types, me)
TB.plot <- data.frame(types, tb)

me_sulfur <- ggplot(data = ME.plot, aes(x = types, y = me)) + geom_bar(stat = "identity", fill = "#5ab4ac", color = "black") + labs(x = "", y = "") + coord_flip() + scale_y_reverse(limits = c(50, 0))

tb_sulfur <- ggplot(data = TB.plot, aes(x = types, y = tb)) + geom_bar(stat = "identity", fill = "#d8b365", color = "black") + labs(x = "", y = "") + coord_flip() + ylim(0, 50)

save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/me_sulfur.pdf", me_sulfur, base_height = 4, base_aspect_ratio = 1.5)
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/tb_sulfur.pdf", tb_sulfur, base_height = 4, base_aspect_ratio = 1.5)

