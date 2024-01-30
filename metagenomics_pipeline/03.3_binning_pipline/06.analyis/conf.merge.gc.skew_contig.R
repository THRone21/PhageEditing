#!/usr/bin/env Rscript

contig=read.table("bin.contig", header=T, sep="\t")
gcskrew=read.table("bin.gc.skew", header=T, sep="\t")
data=merge(contig, gcskrew, by="contig", all=T)
write.table(data, file="bin.gc.skew.conf", quote=F, row.names=F, sep="\t")
