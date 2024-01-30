#!/usr/bin/env Rscript

args = commandArgs(T)
route_file = unlist(strsplit(args[1], "/"))
route = paste(route_file[1:(length(route_file)-1)], collapse="/")
setwd(route)
file_name = route_file[length(route_file)]

library(ggplot2)
data = read.table(file_name, header=T, sep="\t")

data_sort = data[order(data[,3], data[,2], decreasing=F),]
result = ggplot(data_sort, aes(x = func, y = Count, fill=level_1)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    #theme(axis.text.x=element_text(angle=45, hjust=1)) +
    #angle：调整横轴标签倾斜角度
    labs(x = "COG level 2", y = "Count", fill = "COG level 1") +
    scale_x_discrete(limits=factor(data_sort[,1])) +
    #scale_y_continuous(expand=c(0, 0)) +
    theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent'))
ggsave(result, filename=args[2], height = 7, width = 7)
