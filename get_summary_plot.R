setwd("~/Desktop/PhD_Project_related/CR705_WT_P7ko_RNAseq_July2024")


norm_counts=read.table(file="Normalised_CPM_count.txt",header = T,stringsAsFactors = F,sep="\t")
parp7ko_vs_wt=read.csv(file="PARP7KO_vs_WT_use.csv",header = T,stringsAsFactors = F,row.names = 7)
adjp=0.01
parp7ko_vs_wt=parp7ko_vs_wt[,c(1:8)]

#parp7ko_vs_wt=parp7ko_vs_wt[complete.cases(parp7ko_vs_wt),]
parp7ko_vs_wt$baseMean_log=log2(parp7ko_vs_wt$baseMean+1)
parp7ko_vs_wt$Ensembl=rownames(parp7ko_vs_wt)

library(org.Mm.eg.db)

parp7ko_vs_wt$Entrez <- mapIds(org.Mm.eg.db, parp7ko_vs_wt$Ensembl,keytype="ENSEMBL", column="ENTREZID")
parp7ko_vs_wt$Symbol <- mapIds(org.Mm.eg.db, parp7ko_vs_wt$Entrez,keytype="ENTREZID", column="SYMBOL")
parp7ko_vs_wt$Genename <- mapIds(org.Mm.eg.db, parp7ko_vs_wt$Entrez,keytype="ENTREZID", column="GENENAME")

parp7ko_vs_wt$Entrez=as.character(parp7ko_vs_wt$Entrez)
parp7ko_vs_wt$Symbol=as.character(parp7ko_vs_wt$Symbol)
parp7ko_vs_wt$Genename=as.character(parp7ko_vs_wt$Genename)


parp7ko_vs_wt_noGm_riken=parp7ko_vs_wt[- grep("RIKEN",parp7ko_vs_wt$Genename),]
#ahrko_noGm_riken=ahrko[- grep("predicted",ahrko$Gene_Name),]
#parp7ko_vs_wt_noGm_riken=parp7ko_vs_wt_noGm_riken[complete.cases(parp7ko_vs_wt_noGm_riken),]
parp7ko_vs_wt_noGm_riken=parp7ko_vs_wt_noGm_riken[- grep("predicted",parp7ko_vs_wt_noGm_riken$Genename),]
parp7ko_vs_wt_noGm_riken=parp7ko_vs_wt_noGm_riken[- grep("Riken",parp7ko_vs_wt_noGm_riken$Genename),]
parp7ko_vs_wt_noGm_riken=parp7ko_vs_wt_noGm_riken[- grep("pseudogene",parp7ko_vs_wt_noGm_riken$Genename),]



library(dplyr)
library(biomaRt)

#ensembl_symbol=read.table(file="GRCm38_gene_symbol.txt",header = T,sep="\t",stringsAsFactors = F)

#ensembl_symbol=subset(ensembl_symbol,ensembl_symbol$Gene.name!="")
#ensembl_symbol=ensembl_symbol[,c(1,3)]
#ensembl_symbol=ensembl_symbol[complete.cases(ensembl_symbol),]


#parp7ko_vs_wt=left_join(parp7ko_vs_wt,ensembl_symbol,by=c("Ensembl"="Gene.stable.ID"))
#colnames(parp7ko_vs_wt)[9]="Symbol"

#all_combined=rbind(up,down,no_change)
library(ggrepel)
#ggplot(all_combined, aes(avg_exp_log, log2FoldChange)) +
#  geom_point(color = all_combined$Color) +
# theme_classic(base_size = 16)+
# geom_point(data=up_logfc_4, aes(x=avg_exp_log, y=log2FoldChange), colour="red", size=5)

parp7ko_vs_wt_noGm_riken$Neg_log_p_val=-log10(parp7ko_vs_wt_noGm_riken$padj)

parp7ko_vs_wt_sig_up=subset(parp7ko_vs_wt_noGm_riken,parp7ko_vs_wt_noGm_riken$log2FoldChange>1 & parp7ko_vs_wt_noGm_riken$padj < 0.01)
parp7ko_vs_wt_sig_dn=subset(parp7ko_vs_wt_noGm_riken,parp7ko_vs_wt_noGm_riken$log2FoldChange < -1 & parp7ko_vs_wt_noGm_riken$padj < 0.01)

parp7ko_vs_wt_sig_up$parp7ko_vs_wt_Direction="parp7ko_vs_wt_Up"
parp7ko_vs_wt_sig_dn$parp7ko_vs_wt_Direction="parp7ko_vs_wt_Down"

parp7ko_vs_wt_sig=rbind(parp7ko_vs_wt_sig_up,parp7ko_vs_wt_sig_dn)
parp7ko_vs_wt_sig$Ensembl=rownames(parp7ko_vs_wt_sig)

parp7ko_vs_wt_sig=parp7ko_vs_wt_sig[,c("Ensembl","parp7ko_vs_wt_Direction")]

parp7ko_vs_wt_noGm_riken$Ensembl=rownames(parp7ko_vs_wt_noGm_riken)

parp7ko_vs_wt_noGm_riken=left_join(parp7ko_vs_wt_noGm_riken,parp7ko_vs_wt_sig,by=c("Ensembl"))
parp7ko_vs_wt_noGm_riken$parp7ko_vs_wt_Direction[is.na(parp7ko_vs_wt_noGm_riken$parp7ko_vs_wt_Direction)]="No sig change"




#table(parp7ko_vs_wt=parp7ko_vs_wt_sig_up$parp7ko_vs_wt_Direction,`parp7ko_vs_wt Up`=parp7ko_vs_wt_sig_up$parp7ko_vs_wt_Direction)


parp7ko_vs_wt_noGm_riken$baseMean_log=log2(parp7ko_vs_wt_noGm_riken$baseMean+1)



parp7ko_vs_wt_sig_up=subset(parp7ko_vs_wt_noGm_riken,parp7ko_vs_wt_noGm_riken$log2FoldChange>1 & parp7ko_vs_wt_noGm_riken$padj < 0.01)
parp7ko_vs_wt_sig_dn=subset(parp7ko_vs_wt_noGm_riken,parp7ko_vs_wt_noGm_riken$log2FoldChange < -1 & parp7ko_vs_wt_noGm_riken$padj < 0.01)

parp7ko_vs_wt_sig_up$parp7ko_vs_wt_Direction=ifelse(parp7ko_vs_wt_sig_up$log2FoldChange > 0,"Up")
parp7ko_vs_wt_sig_dn$parp7ko_vs_wt_Direction=ifelse(parp7ko_vs_wt_sig_dn$log2FoldChange < 0,"Down")


parp7ko_vs_wt_sig_up=parp7ko_vs_wt_sig_up[order(parp7ko_vs_wt_sig_up$padj),]
parp7ko_vs_wt_sig_dn=parp7ko_vs_wt_sig_dn[order(parp7ko_vs_wt_sig_dn$padj),]

top_15_up=parp7ko_vs_wt_sig_up[1:15,]
top_15_dn=parp7ko_vs_wt_sig_dn[1:15,]

top_15_up_dn=rbind(top_15_up,top_15_dn)

View(parp7ko_vs_wt_new)


library(grid)
library(ggrepel)

p=ggplot(parp7ko_vs_wt_noGm_riken, aes(baseMean_log, log2FoldChange)) +
  theme_classic(base_size = 16)+
  geom_point(data=parp7ko_vs_wt_noGm_riken, aes(x=baseMean_log, y=log2FoldChange), colour="grey", size=2)
p1 <- p +  geom_point(data = parp7ko_vs_wt_sig_up, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="red")
p2 <- p1 +  geom_point(data = parp7ko_vs_wt_sig_dn, aes(x=baseMean_log, y=log2FoldChange) ,size=3,color="blue")
maplot=p2+ggtitle("MA plot for PARP7KO vs WT CR705")+theme(plot.title = element_text(hjust = 0.5))+annotate(geom="text", x=14, y=-17, label="1169 genes up in PARP7KO vs WT CR705",color="red",size=5)+annotate(geom="text", x=14, y=-19, label="461 genes down in PARP7KO vs WT CR705",color="blue",size=5)+xlab("Log2(Mean+1)")+ylim(-27,15)+geom_text_repel(data=top_15_up_dn,aes(x=baseMean_log, y=log2FoldChange,label=Symbol),color="black",arrow=arrow(ends="last",type="open"))+xlim(0,18)


tiff(file="PARP7KO_vsWT_MA_plot.tiff",res=300,height = 3000,width = 3000)
grid.draw(maplot)
dev.off()

pdf(file="PARP7KO_vsWT_MA_plot.pdf",height = 8,width = 9)
grid.draw(maplot)
dev.off()



top_15_up_dn=subset(top_15_up_dn,top_15_up_dn$Ensembl!="ENSMUSG00000121725")

p=ggplot(parp7ko_vs_wt_noGm_riken, aes(log2FoldChange, Neg_log_p_val)) +
  theme_classic(base_size = 16)+
  geom_point(data=parp7ko_vs_wt_noGm_riken, aes(x=log2FoldChange, y=Neg_log_p_val), colour="grey", size=2)
p1 <- p +  geom_point(data = parp7ko_vs_wt_sig_up, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="red")
p2 <- p1 +  geom_point(data = parp7ko_vs_wt_sig_dn, aes(x=log2FoldChange, y=Neg_log_p_val) ,size=3,color="blue")
volcano=p2+ggtitle("Volcano plot for PARP7KO vs WT CR705")+theme(plot.title = element_text(hjust = 0.5))+annotate(geom="text", x=7, y=75, label="1169 genes up in PARP7KO vs WT CR705",color="red",size=5)+annotate(geom="text", x=7, y=80, label="461 genes down in PARP7KO vs WT CR705",color="blue",size=5)+xlab("Log2 Fold Change")+ylab("-Log10(adj. p val)")+ylim(0,96)+geom_text_repel(data=top_15_up_dn,aes(x=log2FoldChange, y=Neg_log_p_val,label=Symbol),color="black",arrow=arrow(ends="last",type="open"))+xlim(-15,12)

volcano

tiff(file="PARP7KO_vsWT_Volcano_plot.tiff",res=300,height = 3000,width = 3500)
grid.draw(volcano)
dev.off()

pdf(file="PARP7KO_vsWT_Volcano_plot.pdf",height = 8,width = 10)
grid.draw(volcano)
dev.off()



#boxplot(parp7ko_vs_wt_sig_up_parp7ko_vs_wt_up$log2FoldChange,parp7ko_vs_wt_sig_up_parp7ko_vs_wt_no_change$log2FoldChange,outline=FALSE)

#dn_tab=parp7ko_vs_wt_sig_dn

#top_15_dn=dn_tab[1:15,]

#top_15_up_dn=rbind(top_15_up,top_15_dn)



#rownames(norm_counts)=norm_counts$X
#norm_counts=norm_counts[,-c(1)]
avg_wt=apply(norm_counts[,grep("WT",colnames(norm_counts))],1,mean)
avg_parp7ko=apply(norm_counts[,grep("AHRKO",colnames(norm_counts))],1,mean)

avg_wt=as.data.frame(avg_wt)
avg_parp7ko=as.data.frame(avg_parp7ko)



rownames(parp7ko_vs_wt_noGm_riken)=parp7ko_vs_wt_noGm_riken$Ensembl

parp7ko_vs_wt_with_exp=cbind(parp7ko_vs_wt_noGm_riken,avg_wt[rownames(parp7ko_vs_wt_noGm_riken),])
parp7ko_vs_wt_with_exp=cbind(parp7ko_vs_wt_with_exp,avg_parp7ko[rownames(parp7ko_vs_wt_with_exp),])
colnames(parp7ko_vs_wt_with_exp)[c(14:15)]=c("wt","parp7ko")

parp7ko_vs_wt_with_exp$wt_log=log2(parp7ko_vs_wt_with_exp$wt+1)
parp7ko_vs_wt_with_exp$parp7ko_log=log2(parp7ko_vs_wt_with_exp$parp7ko+1)



parp7ko_vs_wt_sig_up=subset(parp7ko_vs_wt_with_exp,parp7ko_vs_wt_with_exp$log2FoldChange>1 & parp7ko_vs_wt_with_exp$padj < 0.01)
parp7ko_vs_wt_sig_dn=subset(parp7ko_vs_wt_with_exp,parp7ko_vs_wt_with_exp$log2FoldChange < -1 & parp7ko_vs_wt_with_exp$padj < 0.01)


parp7ko_vs_wt_sig_up$parp7ko_vs_wt_Direction=ifelse(parp7ko_vs_wt_sig_up$log2FoldChange > 0,"Up")
parp7ko_vs_wt_sig_dn$parp7ko_vs_wt_Direction=ifelse(parp7ko_vs_wt_sig_dn$log2FoldChange < 0,"Down")


parp7ko_vs_wt_sig_up=parp7ko_vs_wt_sig_up[order(parp7ko_vs_wt_sig_up$padj),]
parp7ko_vs_wt_sig_dn=parp7ko_vs_wt_sig_dn[order(parp7ko_vs_wt_sig_dn$padj),]

top_15_up=parp7ko_vs_wt_sig_up[1:15,]
top_15_dn=parp7ko_vs_wt_sig_dn[1:15,]

top_15_up_dn=rbind(top_15_up,top_15_dn)



p=ggplot(parp7ko_vs_wt_with_exp, aes(wt_log, parp7ko_log)) +
  theme_classic(base_size = 16)+
  geom_point(data=parp7ko_vs_wt_with_exp, aes(x=wt_log, y=parp7ko_log), colour="grey", size=2)
p1 <- p +  geom_point(data = parp7ko_vs_wt_sig_up, aes(x=wt_log, y=parp7ko_log) ,size=3,color="red")
p2 <- p1 +  geom_point(data = parp7ko_vs_wt_sig_dn, aes(x=wt_log, y=parp7ko_log) ,size=3,color="blue")
correlation_plot=p2+ggtitle("Correlation plot for PARP7KO vs WT in CR705")+theme(plot.title = element_text(hjust = 0.5))+xlab("Log2(Mean+1) WT")+ylab("Log2(Mean+1) PARP7KO")+xlim(0,15)+ylim(0,15)+geom_abline(slope=1, intercept=0,linetype="dotted")+annotate(geom="text", x=5, y=15, label="1169 genes up after PARP7KO vs wt WT",color="red",size=6)+annotate(geom="text", x=5, y=14, label="461 genes down after PARP7KO vs WT",color="blue",size=6)+geom_text_repel(data=top_15_up_dn,aes(x=wt_log, y=parp7ko_log,label=Symbol),color="black",arrow=arrow(ends="last",type="open"),fontface="bold")


tiff(file="PARP7KO_vsWT_Correlation_plot.tiff",res=300,height = 3000,width = 3500)
grid.draw(correlation_plot)
dev.off()


pdf(file="PARP7KO_vsWT_Correlation_plot.pdf",height = 8,width = 8)
grid.draw(correlation_plot)
dev.off()




DEG_list=rbind(parp7ko_vs_wt_sig_dn,parp7ko_vs_wt_sig_up)

write.csv(DEG_list,file="PARP7KO_vs_WT_CR705_DEG_list.csv",col.names = T,row.names = F,quote = F)


