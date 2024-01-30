#!/usr/bin/env Rscript

args = commandArgs(T)
route_file = unlist(strsplit(args[1], "/"))
route = paste(route_file[1:(length(route_file)-1)], collapse="/")
setwd(route)
file_name = route_file[length(route_file)]

data = read.table(file_name, sep="\t", header=T)
db = read.table(args[2], sep="\t", header=T)
result = merge(data, db, by="func", all.x=T)

write.table(result, file=args[3], quote=F, row.names=F, sep="\t")
