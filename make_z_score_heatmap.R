setwd("~/Desktop/PhD_Project_related/CR705_WT_P7ko_RNAseq_July2024")

library(dplyr)
`%ni%`=Negate(`%in%`)

parp7ko=read.csv(file="PARP7KO_vs_WT.csv",header = T,stringsAsFactors = F)

parp7ko_up=subset(parp7ko,parp7ko$padj < 0.01 & parp7ko$log2FoldChange > 1)
parp7ko_dn=subset(parp7ko,parp7ko$padj < 0.01 & parp7ko$log2FoldChange < -1)

parp7ko_up=parp7ko_up[order (-parp7ko_up$log2FoldChange),]
parp7ko_dn=parp7ko_dn[order (parp7ko_dn$log2FoldChange),]


parp7ko_up_dn=c(parp7ko_dn$Ensembl,parp7ko_up$Ensembl)

cpm=read.table(file="Normalised_CPM_count.txt",header = T,sep="\t",stringsAsFactors = F)

cpm_parp7ko=cpm[,grep("K",colnames(cpm))]
cpm_wt=cpm[,grep("W",colnames(cpm))]


colnames(cpm_parp7ko)[1:9]=c(paste("PARP7KO",1:9,sep="_"))
colnames(cpm_wt)[1:7]=c(paste("WT",1:7,sep="_"))



library(pheatmap)
library(RColorBrewer)

cpm_parp7ko_plot=cpm_parp7ko[parp7ko_up_dn,]
cpm_wt_plot=cpm_wt[parp7ko_up_dn,]

cpm_parp7ko_wt=cbind(cpm_wt_plot,cpm_parp7ko_plot)


heatmap_parp7ko_wt=pheatmap(cpm_parp7ko_wt,scale = "row",show_rownames = FALSE,main="Z score heatmap of DEG's in PARP7KO vs WT in CR705 mice (N=1671)",cluster_cols = FALSE,cluster_rows = FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))


pdf(file="PARP7KO_vsWT_CR705_DEG_Zscore_heatmap_ordered.pdf",height = 8,width = 8)
grid.draw(heatmap_parp7ko_wt)
dev.off()


