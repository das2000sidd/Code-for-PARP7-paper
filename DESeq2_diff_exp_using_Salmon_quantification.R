setwd("~/Desktop/PhD_Project_related/CR705_WT_P7ko_RNAseq_July2024")


library(tximportData)

library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(BSgenome.Mmusculus.UCSC.mm10)

txdb = TxDb.Mmusculus.UCSC.mm39.refGene
keytypes(txdb)
columns(txdb)

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(biomaRt)
#ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#head(listAttributes(ensembl))
#myAttributes <- c("ensembl_gene_id","ensembl_transcript_id_version")

#res <- getBM(attributes =  myAttributes,
     #        mart = ensembl)


library(tibble)

res=read.table(file="GRCm39_gene_and_transcript_stable_ID_version.txt",header = T,sep="\t",stringsAsFactors = F)
res=res[,c(1,4)]
colnames(res)=c("GENEID","TXNAME")
res=res[,c(2,1)]



res_tibble=as_tibble(res)
#sample_id=c("KO_1","KO_2","KO_3","KO_4","KO_5","KO2W_1","KO2W_2","KO2W_3","KO2W_4","KO2W_5","WT_1","WT_2","WT_3","WT_4","WT_5","WT2W_1","WT2W_2","WT2W_3","WT2W_4","WT2W_5")
files=list.files(path="/Users/siddhaduio.no/Desktop/PhD_Project_related/CR705_WT_P7ko_RNAseq_July2024",full.names = T)



files_G12=files[grep("G12",files)]
files_G13=files[grep("G13",files)]

all_files=c(files_G12,files_G13)

files=as.data.frame(all_files)



files_quant=paste(files$all_files,"/quant.sf",sep="")


library(tximport)

sample_order=c(paste("WT",1:7,sep=""),paste("PARP7KO",1:9,sep=""))

txi = tximport(files_quant,type = "salmon",tx2gene = res)
names(txi)
head(txi$counts)
colnames(txi$counts)=sample_order
head(txi$abundance)
colnames(txi$abundance)=sample_order


write.table(txi$counts,file="Salmon_Counts_Object.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(txi$abundance,file="Salmon_Abundance_Object.txt",col.names = T,row.names = T,sep="\t",quote = F)


library(DESeq2)

genotype_info=c(rep("WT",7),rep("PARP7KO",9))
genotype_info=as.factor(genotype_info)
genotype_info=as.data.frame(genotype_info)
model_mat_geno_treat=model.matrix( ~ 0 + genotype_info,data=genotype_info)

dds <- DESeqDataSetFromTximport(txi, genotype_info, model_mat_geno_treat)



## Exploratory analysis and visualisation
nrow(dds)


# *** variance stabilizing transformation and the rlog***

## Showing effect on some simulated data

lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)

log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)


vsd <- rlog(dds, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)



library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)


library("pheatmap")
library("RColorBrewer")


###PCA plot***
tiff(file="PCA plot.tiff",res=300,height = 2000,width = 2000)
plotPCA(vsd, intgroup = c("genotype_info")) +ggtitle("")+theme(axis.text = element_text(size=12))+theme(axis.title = element_text(size=12))
dev.off()


pdf(file="PCA plot.pdf",height = 8,width = 8)
plotPCA(vsd, intgroup = c("genotype_info")) +ggtitle("")+theme(axis.text = element_text(size=12))+theme(axis.title = element_text(size=12))
dev.off()



dds$group <- factor(paste0(dds$genotype_info))
design(dds) = ~ group
dds <- DESeq(dds)

contrast_wt_vs_parp7ko <- c("group", "PARP7KO", "WT")

res_wt_vs_parp7ko <- results(dds, contrast=contrast_wt_vs_parp7ko,pAdjustMethod = "BH",format = "DataFrame")
res_wt_vs_parp7ko_df=as.data.frame(res_wt_vs_parp7ko)
res_wt_vs_parp7ko_df$Ensembl=rownames(res_wt_vs_parp7ko_df)
summary(res_wt_vs_parp7ko)
res_wt_vs_parp7ko_df=res_wt_vs_parp7ko_df[complete.cases(res_wt_vs_parp7ko_df),]
View(res_wt_vs_parp7ko_df)



library(org.Mm.eg.db)

res_wt_vs_parp7ko_df$Entrez <- mapIds(org.Mm.eg.db, res_wt_vs_parp7ko_df$Ensembl,keytype="ENSEMBL", column="ENTREZID")
res_wt_vs_parp7ko_df$Symbol <- mapIds(org.Mm.eg.db, res_wt_vs_parp7ko_df$Entrez,keytype="ENTREZID", column="SYMBOL")
res_wt_vs_parp7ko_df$Genename <- mapIds(org.Mm.eg.db, res_wt_vs_parp7ko_df$Entrez,keytype="ENTREZID", column="GENENAME")

#res_wt_vs_parp7ko_df$Entrez=as.character(res_wt_vs_parp7ko_df$Entrez)
res_wt_vs_parp7ko_df$Symbol=as.character(res_wt_vs_parp7ko_df$Symbol)
res_wt_vs_parp7ko_df$Genename=as.character(res_wt_vs_parp7ko_df$Genename)


parp7ko_vs_wt_sig=subset(res_wt_vs_parp7ko_df,abs(res_wt_vs_parp7ko_df$log2FoldChange) > 1 & res_wt_vs_parp7ko_df$padj < 0.01)

View(parp7ko_vs_wt_sig)


write.csv(parp7ko_vs_wt_sig,file="PARP7KO_vs_WT_sig_Genes.csv",col.names = T,row.names = F,quote = F)
write.csv(res_wt_vs_parp7ko_df,file="PARP7KO_vs_WT.csv",col.names = T,row.names = F,quote = F)



library(edgeR)
cpm.nor.count=cpm(dds,normalized.lib.sizes = TRUE,log=FALSE)

write.table(cpm.nor.count,file="Normalised_CPM_count.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(rownames(dds),file="Genes_used.txt",col.names = T,row.names = T,sep="\t",quote = F)






