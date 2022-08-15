#!/usr/bin/env Rscript
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1])
outfn = args[2]

x = x[complete.cases(x), ]

insertRow <- function(existingDF, new_dataframe, index) {
  colnames(new_dataframe) = c("variant",
                              "gene",
                              "pval")
  existingDF = rbind(existingDF,new_dataframe)
  existingDF
}

outdf = data.frame(
  variant=complex(),
  gene=complex(),
  pval=numeric()
)
  
current_index = 1

all_variants_genes = unique(x[,c(1,2)])


for (k in rownames(all_variants_genes)) {
  tmp = all_variants_genes[c(k),]
  variant = tmp$V1
  gene = tmp$V2
  tmp = subset(x,V1 == variant)
  y = subset(tmp,V2 == gene)
  y = dcast(y, V1 + V2 + V3 ~ V4, value.var="V5")
  rownames(y) = y$V3
  y = subset(y, select=-c(V1,V2,V3))
  if (ncol(y) == 1) {
    next
  }
  if (nrow(y) == 1) {
    next
  }
  print(variant)
  print(gene)
  fish <- fisher.test(y,workspace = 2e8)
  variant = as.character(variant)
  gene = as.character(gene)
  add_df = data.frame(variant=c(variant),gene=c(gene),pval=c(fish$p.value))
  outdf = insertRow(outdf,add_df,current_index)
  current_index = current_index + 1
}

outdf = apply(outdf,2,as.character)
write.table(outdf,file=outfn,quote=FALSE,sep="\t",row.names=FALSE)
