library(ggplot2)
library("ggsci")
library(reshape2)
library(stringr)
library(ggpmisc)

args = commandArgs(trailingOnly=TRUE)

setwd(args[1])

lead_guqtl_gt_usage = read.table(args[2],sep="\t",header=FALSE)
M_IGHV1_2_132_0 = read.table(args[3],header=TRUE)
M_IGHV1_2_132_0_gt_usage = read.table(args[4],header=FALSE,sep="\t")
M_IGHV3_66_10_0 = read.table(args[5],header=TRUE)
gene_count_per_variant = read.table(args[6])
gutql = read.table(args[7],header=T)
genotype_usage = read.table(args[8])

panel_a_1 = function(x) {
    x = subset(x,V5 != "NA")
    y = subset(x,V2 == "IGHV1-2")
    ggplot(y,aes(V5,V3)) +
        geom_boxplot() +
        facet_wrap(V2 ~ .,scales = "free",nrow=1) +
        scale_y_continuous(lim=c(0,NA)) +
        theme_bw() 
    ggsave("figures/ighv_1_2_guqtl_boxplot.pdf",width=2,height=2)
}

panel_a_2 = function(x) {
    y = x
    y$variant = as.numeric(y$variant)
    ggplot(y,aes(variant,-log10(pvalue),color=sign_corr)) +
        geom_point(size=1.25) +
        theme_bw() +
        scale_color_manual(name ="P-value < 1e-5",
                           values=c("#FF6F00FF","#C71000FF"))
  ggsave("figures/ighv1_2_cond_manh_plot.pdf",width=8,height=2)
}

panel_a_3 = function(x) {
    x = subset(x,V1 == "115056")
    x = subset(x,V5 != "NA")
    y = subset(x,V2 == "IGHV1-2")
    y$V3 = as.numeric(as.character(y$V3))
    y$V5 = as.character(y$V5)
    ggplot(y,aes(V5,V3)) +
        geom_boxplot() +
        facet_wrap(V2 ~ .,scales = "free",nrow=1) +
        scale_y_continuous(lim=c(0,NA)) +
        theme_bw() 
    ggsave("figures/ighv1_2_cond_guqtl_boxplot.pdf",width=1.6,height=2)
}

panel_b = function(y) {
    y$variant = as.numeric(y$variant)
    ggplot(y,aes(variant,-log10(pvalue),color=sign_corr)) +
        geom_point(size=1.25) +
        theme_bw() +
        scale_color_manual(name ="P-value < 1e-5",
                           values=c("#FF6F00FF","#C71000FF"))
    ggsave("figures/ighv3_66_cond_manh_plot.pdf",width=8,height=2)
}


panel_c = function(x) {
    ggplot(x,aes(V2)) +
        geom_histogram() +
        theme_bw()
    ggsave("figures/snp_per_gene.pdf",width = 2.5,height=2.5)
}

panel_f = function(x) {
    y = subset(x,sign_corr == "True")
    y$variant = as.numeric(as.character(y$variant))
    y = subset(y,gene %in% c("IGHV1-69/-69D","IGHV3-53","IGHV4-31",
                             "IGHV3-64","IGHV3-66","IGHV4-59","IGHV4-61"))
    ggplot(y,aes(variant,-log10(pvalue),color=gene)) +
        geom_point(size=1,alpha=0.8) +
        ylim(0,NA) +
        geom_vline(xintercept = 905253,linetype = 2) +
        theme_bw() +
        scale_color_futurama(name="QTL Genes")
    ggsave("figures/snps_locus.pdf",height=2.5,width=8)
}

panel_g = function(x) {
    x = subset(x,V1 == 905253)
    x = subset(x,V5 != "NA")
    x$V2 = factor(x$V2, levels=c("IGHV4-31","IGHV3-53","IGHV4-59","IGHV4-61",
                                 "IGHV3-64","IGHV3-66","IGHV1-69/-69D"))
    x = subset(x,V2 != "NA")
    #x = subset(x,V2 != NA)
    ggplot(x,aes(as.factor(V5),V3)) +
        geom_boxplot() +
        facet_wrap(.~ V2,nrow=1,scales="free") +
        scale_y_continuous(lim=c(0,NA)) +
        theme_bw()
    ggsave("figures/snps_genes_multi_box_plot.pdf",height=2.5,width=12)
}

panel_a_1(lead_guqtl_gt_usage)
panel_a_2(M_IGHV1_2_132_0)
panel_a_3(M_IGHV1_2_132_0_gt_usage)
panel_b(M_IGHV3_66_10_0)
panel_c(gene_count_per_variant)
panel_f(gutql)
panel_g(genotype_usage)
