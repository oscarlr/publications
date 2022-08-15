library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

setwd(args[1])

allele_count_pvals = read.table(args[2],header=T)

panel_a = function(x) {
    x = subset(x, select = -c(variant))
    x = unique(x)
    x$sign = "Yes"
    x$sign[x$pval > 0.05] = "No"
    x$pval_log10 = -log10(as.numeric(as.character(x$pval)))
    ggplot(x,aes(reorder(gene,pval),pval_log10,fill=sign)) +
        geom_bar(stat='identity') +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_fill_manual(name ="Significant",
                          values=c("#84D7E1FF","#FF6348FF"))
    ggsave("figures/allele_diff.pdf",height=3,width=8)
}

panel_a(allele_count_pvals)
