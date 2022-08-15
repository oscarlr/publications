library(ggplot2)
library("ggsci")
library(reshape2)
library(stringr)
library(ggpmisc)

args = commandArgs(trailingOnly=TRUE)

supp = read.table(args[1],sep="\t",header=TRUE)
supp = subset(supp,isotype == "M")

#x =  "/home/o0rodr03/projects/Genetic_effect_on_AbRep/scratch/2022-08-04_fixing_eqtl_code/parse_eqtl/M/all_rank.txt"
guqtl = read.table(args[2],header=TRUE)

setwd(args[3])

panel_a_1 = function(x) {
    x = x[,c("gene","gene_start","num_sign_variants"),]
    x = unique(x)
    ggplot(x,aes(reorder(gene,gene_start),num_sign_variants,fill="FALSE")) +
        geom_bar(stat="identity") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6)) +
        scale_fill_manual(name ="Source",values=c("#3F4041FF")) 
    ggsave("figures/eqtl_variant_count.pdf",width=10,height=2.5)
}

panel_a_2 = function(x) {
    x = x[,c("gene","gene_start","anova_pval","sign_corr"),]
    x = unique(x)
    ggplot(x,aes(reorder(gene,gene_start),-log10(anova_pval),fill=sign_corr)) +
        geom_bar(stat="identity") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6)) +
        scale_fill_futurama(name = "P-value < 9e-6")
    ggsave("figures/top_pval_gene.pdf",width=10,height=2.5)
}

panel_a_3 = function(x) {
    x = x[,c("gene","gene_start","linear_r","sign_corr"),]
    x = unique(x)
    x$r_square = as.numeric(as.character(x$linear_r))
    x$r_square[x$sign_corr == "False"] = NA
    ggplot(x,aes(reorder(gene,gene_start),r_square,fill="FALSE")) +
        geom_bar(stat="identity") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6)) +
        scale_fill_manual(name ="Source",values=c("#3F4041FF"))
    ggsave("figures/top_pval_r_square.pdf",width=10,height=2.5)
}

panel_a_4 = function(x) {
    x = x[,c("gene","gene_start","max_min_fold_changes","sign_corr"),]
    x = unique(x)
    x$fc = as.numeric(as.character(x$max_min_fold_changes))
    x$fc[x$max_min_fold_changes == "inf"] = 15    
    x$fc[x$sign_corr == "False"] = NA
    x$fc_cat = "[1,2)"
    x$fc_cat[x$fc >=2] = "[2,15)"
    x$fc_cat[x$fc >= 15] = "[15,inf)"
    x$fc[x$fc >= 15] = 15
    ggplot(x,aes(reorder(gene,gene_start),fc,fill=fc_cat)) +
        geom_bar(stat="identity") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6)) +
        scale_fill_manual(name ="Fold-change",
                          values=c("#008EA0FF","#FF6F00FF","#C71000FF"),
                          breaks=c("[1,2)","[2,15)","[15,inf)"),
                          labels=c("[1,2)","[2,15)","[15,inf)"))
    ggsave("figures/top_pval_fc.pdf",width=10,height=2.5)
}

panel_a_5 = function(x) {
    x = x[,c("gene","gene_start","variant_type","sign_corr"),]
    x = unique(x)
    x$variant_c = 1
    x$variant_c[x$sign_corr == "False"] = NA
    ggplot(x,aes(reorder(gene,gene_start),variant_c,fill=variant_type)) +
        geom_bar(stat="identity") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6)) +
        scale_fill_futurama(name = "Variant type",alpha=.9)
    ggsave("figures/variant_type_1.pdf",width=10,height=1.75) 
    ggsave("figures/variant_type_2.pdf",width=10,height=4)
}

panel_e = function(x) {
    y = subset(x,gene == "IGHV3-66")
    y$variant = as.numeric(as.character(y$variant))
    ggplot(y,aes(variant,-log10(pvalue),color=sign_corr)) +
        geom_point(size=1.25,alpha=.5) +
        theme_bw() +
        scale_color_manual(name="P-value < 9e-6",
                           values=c("#FF6F00FF","#C71000FF"))
    ggsave("figures/ighv3_66.pdf",width=6,height=1.75)
}


panel_a_1(supp)
panel_a_2(supp)
panel_a_3(supp)
panel_a_4(supp)
panel_a_5(supp)
panel_e(guqtl)
