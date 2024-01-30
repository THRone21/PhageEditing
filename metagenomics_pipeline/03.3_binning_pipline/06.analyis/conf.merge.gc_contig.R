#!/usr/bin/env Rscript

contig=read.table("bin.contig", sep="\t", header=T)
gc=read.table("bin.gc", sep="\t", header=T)
data=merge(contig, gc, by="contig", all=T)
write.table(data, file="bin.gc.conf", quote=F, row.names=F, sep="\t")
