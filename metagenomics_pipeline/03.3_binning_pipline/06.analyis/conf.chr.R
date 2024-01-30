#!/usr/bin/env Rscript

data=read.table("bin.contig", header=T, sep="\t")
data$chr=rep("chr", length(data[,1]))
data$bar=rep("-", length(data[,1]))
data$rank=seq(1:length(data[,1]))
data$col=rep("purple", length(data[,1]))  # 添加新列

data=data[, c("chr", "bar", "contig", "rank", "start", "end", "col")]# 指定列顺序

colnames(data)[1]="#chr"
# 注销第一列

write.table(data, file="bin.karyotype", quote=F, row.names=F, sep="\t")
