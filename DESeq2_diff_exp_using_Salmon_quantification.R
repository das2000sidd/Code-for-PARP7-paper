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




#library(tximport)
#txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
#names(txi)

files_G12=files[grep("G12",files)]
files_G13=files[grep("G13",files)]

all_files=c(files_G12,files_G13)

files=as.data.frame(all_files)



files_quant=paste(files$all_files,"/quant.sf",sep="")




library(tximport)

sample_order=c(paste("WT",1:7,sep=""),paste("PARP7KO",1:9,sep=""))

#files=c("./ID86V4/quant.sf","./ID86V0/quant.sf","./ID82V4/quant.sf","./ID75V4/quant.sf","./ID75V0/quant.sf")
#files=all_Samples_ordered$files
#res_tibble=res_tibble[,c(2,1)]
txi = tximport(files_quant,type = "salmon",tx2gene = res)
names(txi)
head(txi$counts)
colnames(txi$counts)=sample_order
head(txi$abundance)
colnames(txi$abundance)=sample_order


write.table(txi$counts,file="Salmon_Counts_Object.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(txi$abundance,file="Salmon_Abundance_Object.txt",col.names = T,row.names = T,sep="\t",quote = F)


#sampleTable=data.frame(condition = factor(rep(c("WT_2_week","KO_30_days","KO_2_weeks","WT_3_week"),each=5)))
#sampleTable=sampleTable[-c(15),]
#sampleTable=as.data.frame(sampleTable)
#rownames(sampleTable)=colnames(txi$counts)

library(DESeq2)
#library(edgeR)


#counts_name_order=colnames(txi$counts)

#sample_info=read.csv(file="Sample_info_Table.csv",header = T,stringsAsFactors = F)

#sample_info[2:12,1]="PyMT Px459 empty"
#sample_info[14:24,1]="PyMT parp7ko c17"

#wt=sample_info[1:12,]
#parp7ko=sample_info[13:24,]

#wt$Genotype="WT"
#parp7ko$Genotype="parp7ko"

#sample_info=rbind(wt,parp7ko)

#sample_info$Genotype_treatment = paste(sample_info$Treatment,sample_info$Genotype,sep="_")



#model_mat_reduced=model.matrix( ~ Genotype + Treatment,data=sample_info)
#sample_order=colnames(txi$counts)
#rownames(sample_info)=sample_info$ID.for.seq
#sample_info_ordered=sample_info[sample_order,]
genotype_info=c(rep("WT",7),rep("PARP7KO",9))
genotype_info=as.factor(genotype_info)
genotype_info=as.data.frame(genotype_info)
model_mat_geno_treat=model.matrix( ~ 0 + genotype_info,data=genotype_info)

#rownames(model_mat_geno_treat)=counts_name_order


library(DESeq2)
library(edgeR)
dds <- DESeqDataSetFromTximport(txi, genotype_info, model_mat_geno_treat)

idx <- rowSums( cpm(dds) > 2) >= 2 ## This is not necessary

#head(all_counts_combined[,c(6:8)])

# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
#meanSdPlot(assay(vsd))
#meanSdPlot(assay(rld))

## Exploratory analysis and visualisation
nrow(dds)

#keep = rowSums( cpm(dds) >2) >= 2 ## at least two samples with count of 1 or higher
#dds = dds[keep,]
#nrow(dds) ## 20563

#keep <- filterByExpr(dds, group = condition_df$condition)
#dds = dds[keep,]


#fpm_all_genes=fpm(dds)
#avg_fpm=apply(fpm_all_genes, 1, mean)
#length(which(avg_fpm>2))

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


rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)


library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  
# the VST has a upward shift for the smaller values



#Sample distances#


sampleDists <- dist(t(assay(vsd)))
sampleDists


library("pheatmap")
library("RColorBrewer")


sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$genotype_info, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)



library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$genotype_info)
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)




###PCA plot***
tiff(file="PCA plot.tiff",res=300,height = 2000,width = 2000)
plotPCA(vsd, intgroup = c("genotype_info")) +ggtitle("")+theme(axis.text = element_text(size=12))+theme(axis.title = element_text(size=12))
dev.off()


pdf(file="PCA plot.pdf",height = 8,width = 8)
plotPCA(vsd, intgroup = c("genotype_info")) +ggtitle("")+theme(axis.text = element_text(size=12))+theme(axis.title = element_text(size=12))
dev.off()



## genotype is primary driver of difference in gene expression. Not treatment.
## Therefore we ask what is the difference between WT + DMSO vs parp7ko + DMSO
# WT + Bay vs parp7ko+Bay
# WT + Kyn vs parp7ko+Kyn

###MDS plot***

mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = genotype_info)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS Plot")  


## MDS plot using the VST data


mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color =genotype_info)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")



### SVASeq
#### Running the differential expression pipeline
## we can run the differential expression pipeline on the raw counts with a single call to the function DESeq



What above command does:
  using pre-existing size factors
estimating dispersions
found already estimated dispersions, replacing these
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("genotype_info")])
top_20_genes=assay(ntd)[select,]
colnames(df)="Condition"
rownames(df)=colnames(top_20_genes)
pheatmap(top_20_genes, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


library(Glimma)
#dds <- DESeq(dds,test = "LRT",reduced = model_mat_reduced)
## Building the results table
## 

dds$group <- factor(paste0(dds$genotype_info))
design(dds) = ~ group
dds <- DESeq(dds)
#plotMDS(dds)
#glimmaMDS(dds)


#Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula

#res <- results(dds)
#res
# 1. Create a contrast vector called contrast_kd.
#contrast_bay_dmso <- c("group", "WTBay", "WTDMSO")

#res_Bay_DMSO <- results(dds, contrast=contrast_bay_dmso,pAdjustMethod = "BH",format = "DataFrame")
#mcols(res_Bay, use.names = TRUE)
#res_Bay_DMSO_df=as.data.frame(res_Bay_DMSO)
#res_Bay_DMSO_df$Ensembl=rownames(res_Bay_DMSO_df)
#summary(res_Bay_DMSO)
#res_Bay_DMSO_df=res_Bay_DMSO_df[complete.cases(res_Bay_DMSO_df),]


#contrast_kyn_dmso <- c("group", "WTKyn", "WTDMSO")

#res_Kyn_DMSO <- results(dds, contrast=contrast_kyn_dmso,pAdjustMethod = "BH",format = "DataFrame")
#res_Kyn_DMSO_df=as.data.frame(res_Kyn_DMSO)
#res_Kyn_DMSO_df$Ensembl=rownames(res_Kyn_DMSO_df)
#summary(res_Bay_DMSO)
#res_Kyn_DMSO_df=res_Kyn_DMSO_df[complete.cases(res_Kyn_DMSO_df),]

## Use 1
contrast_wt_vs_parp7ko <- c("group", "PARP7KO", "WT")

res_wt_vs_parp7ko <- results(dds, contrast=contrast_wt_vs_parp7ko,pAdjustMethod = "BH",format = "DataFrame")
res_wt_vs_parp7ko_df=as.data.frame(res_wt_vs_parp7ko)
res_wt_vs_parp7ko_df$Ensembl=rownames(res_wt_vs_parp7ko_df)
summary(res_wt_vs_parp7ko)
res_wt_vs_parp7ko_df=res_wt_vs_parp7ko_df[complete.cases(res_wt_vs_parp7ko_df),]
View(res_wt_vs_parp7ko_df)

# Use 2


## Ninni suggested the cutoff


library(org.Mm.eg.db)

res_wt_vs_parp7ko_df$Entrez <- mapIds(org.Mm.eg.db, res_wt_vs_parp7ko_df$Ensembl,keytype="ENSEMBL", column="ENTREZID")
res_wt_vs_parp7ko_df$Symbol <- mapIds(org.Mm.eg.db, res_wt_vs_parp7ko_df$Entrez,keytype="ENTREZID", column="SYMBOL")
res_wt_vs_parp7ko_df$Genename <- mapIds(org.Mm.eg.db, res_wt_vs_parp7ko_df$Entrez,keytype="ENTREZID", column="GENENAME")

#res_wt_vs_parp7ko_df$Entrez=as.character(res_wt_vs_parp7ko_df$Entrez)
res_wt_vs_parp7ko_df$Symbol=as.character(res_wt_vs_parp7ko_df$Symbol)
res_wt_vs_parp7ko_df$Genename=as.character(res_wt_vs_parp7ko_df$Genename)


parp7ko_noGm_riken=res_wt_vs_parp7ko_df[- grep("RIKEN",res_wt_vs_parp7ko_df$Genename),]
#parp7ko_noGm_riken=parp7ko_noGm_riken[complete.cases(parp7ko_noGm_riken),]
parp7ko_noGm_riken=parp7ko_noGm_riken[- grep("predicted",parp7ko_noGm_riken$Genename),]
parp7ko_noGm_riken=parp7ko_noGm_riken[- grep("Riken",parp7ko_noGm_riken$Genename),]
parp7ko_noGm_riken=parp7ko_noGm_riken[- grep("pseudogene",parp7ko_noGm_riken$Genename),]

parp7ko_vs_wt_sig=subset(parp7ko_noGm_riken,abs(parp7ko_noGm_riken$log2FoldChange) > 1 & parp7ko_noGm_riken$padj < 0.01)

View(parp7ko_vs_wt_sig)


write.csv(parp7ko_vs_wt_sig,file="PARP7KO_vs_WT_1630_sig_Genes.csv",col.names = T,row.names = F,quote = F)



library(edgeR)
cpm.nor.count=cpm(dds,normalized.lib.sizes = TRUE,log=FALSE)

write.table(cpm.nor.count,file="Normalised_CPM_count.txt",col.names = T,row.names = T,sep="\t",quote = F)
write.table(rownames(dds),file="Genes_used.txt",col.names = T,row.names = T,sep="\t",quote = F)

#library(biomaRt)

#ensembl_symbol=read.table(file="hg38_ensembl_symbol_entrez_id.txt",header = T,sep="\t",stringsAsFactors = F)

#ensembl_symbol=subset(ensembl_symbol,ensembl_symbol$Gene.name!="")
#ensembl_symbol=ensembl_symbol[,c(1,3)]
#ensembl_symbol=subset(ensembl_symbol,ensembl_symbol$HGNC.symbol!="")






write.table(rownames(dds),file="Background_genes.txt",col.names = F,row.names = F,sep="\t",quote = F)



