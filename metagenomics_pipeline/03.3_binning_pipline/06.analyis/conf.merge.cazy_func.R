#!/usr/bin/env Rscript

cazy=read.table("bin.cazy", header=T, sep=" ")
func=read.table("bin.func", header=T, sep=" ")
data=merge(func, cazy, by="id", all=T)
data2=data[, c(2,3,4,5)]
write.table(data2, file="bin.cazy.conf", quote=F, row.names=F, sep="\t")

