#!/usr/bin/env Rscript

#this script merge annotation outputs from Eggnog and Interposcan
#USAGE
# Rscript --vanilla merge_annotation.R eggnog.output interpro_Pased.output output/directory



args = commandArgs(trailingOnly=TRUE)

eggnog<- read.delim(args[1], header = T)
interpro<- read.delim(args[2], header = T)
colnames(interpro)[1]<- "gene_id"

merge(x =eggnog , y =interpro , by = "gene_id", all = TRUE)-> tabla_final

b<-strsplit(args[1],"\\/")[[1]]
c <- strsplit(tail(b, n=1), "\\.")[[1]]
name<-paste(c[1],"annotation.tsv",sep="_")

write.table(tabla_final,
            file= paste(args[3],name, sep= "/"), 
            sep="\t",
            row.names = F, 
            col.names = T) 
