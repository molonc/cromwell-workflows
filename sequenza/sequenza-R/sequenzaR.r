#!/usr/bin/env Rscript

library(sequenza)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#data.file<-system.file("data", args[1], package="sequenza")
test<-sequenza.extract(args[1],verbose=FALSE)
CP<-sequenza.fit(test)

sequenza.results(sequenza.extract=test,cp.table=CP,sample.id=args[2],out.dir="TEST")