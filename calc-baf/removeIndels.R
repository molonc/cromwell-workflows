#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

bed <- as.data.frame(read.table(args[1],header=FALSE, sep="\t",stringsAsFactors=FALSE))

a <- c("A", "C", "G", "T")
newbed <- subset(bed, V3%in%a & V4%in%a)

metrics <- c(dim(bed)[1], dim(newbed)[1], dim(newbed)[1]/dim(bed)[1])

write.table(newbed,"tmp.baf.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(metrics,"metrics.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)



