library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

setwd(args[1])

enrichment = read.table(args[2],sep="\t",header=FALSE)

panel_a = function(x) {
    x = subset(x,V8 == "ccre")
    x = subset(x,V1 != "ccre")
    x$guqtl_frac = x$V2 / x$V3
    x$all_frac = x$V4 / x$V5
    y = melt(x, id.vars=c("V1", "V2","V3","V4","V5","V6","V7","V8"))
    ggplot(y,aes(reorder(V1,V7),value,fill=variable)) +
        geom_bar(stat='identity',position="dodge") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_fill_manual(name ="Variants",
                          labels=c("guQTL","non-guQTL"),
                          values=c("#84D7E1FF","#FF6348FF"))
    ggsave("figures/ccre_enrich.pdf",height=4,width=5)
}

panel_b = function(x) {
    x = subset(x,V8 == "tfbs")
    x = subset(x,V1 != "tfbs")
    x = subset(x,V7 < .1)
    x$guqtl_frac = x$V2 / x$V3
    x$all_frac = x$V4 / x$V5
    y = melt(x, id.vars=c("V1", "V2","V3","V4","V5","V6","V7","V8"))
    ggplot(y,aes(reorder(V1,V7),value,fill=variable)) +
        geom_bar(stat='identity',position="dodge") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        scale_fill_manual(name ="Variants",
                          labels=c("guQTL","non-guQTL"),
                          values=c("#84D7E1FF","#FF6348FF"))
    ggsave("figures/tfbs_enrich.pdf",height=4,width=5)
}

panel_a(enrichment)
panel_b(enrichment)
