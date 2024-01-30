#!/usr/bin/env Rscript

library("png")

legend=readPNG("legend.png")

files=list.files(pattern="bin.*.png")
id=vector()
for(i in 1:length(files))
{
    id[i]=as.character(strsplit(files[i], split=".png"))
}

img=list()
for(i in 1:length(files))
{
    img[[i]]=readPNG(files[i])
}

for(i in 1:length(files))
{
    pdf(paste(id[i], "circos.pdf", sep="."), width=7, height=14)
    
    opar=par(no.readonly=TRUE)
    par(mfrow=c(2, 1))  # 填充：2行，1列
    par(mar=rep(0, 4))  # 去掉边空(英分)
    
    plot(NA, xlim=0:1, ylim=0:1, xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    rasterImage(img[[i]], 0, 0, 1, 1)
    plot(NA, xlim=0:1, ylim=0:1, xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
    rasterImage(legend, 0, 0, 1, 1)

    par(opar)
    dev.off()
}
