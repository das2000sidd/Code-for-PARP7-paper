setwd("~/Desktop/PhD_Project_related/CR705_WT_P7ko_RNAseq_July2024")


library(dplyr)
`%ni%` = Negate(`%in%`)


parp7ko_vs_wt=read.csv(file="PARP7KO_vs_WT.csv",header = T,stringsAsFactors = F,row.names = 7)


parp7ko_vs_wt_up=subset(parp7ko_noGm_riken,parp7ko_noGm_riken$log2FoldChange > 1 & parp7ko_noGm_riken$padj < 0.01)

parp7ko_vs_wt_up_g=parp7ko_vs_wt_up$Symbol

parp7ko_vs_wt_dn=subset(parp7ko_noGm_riken,parp7ko_noGm_riken$log2FoldChange < -1 & parp7ko_noGm_riken$padj < 0.01)

parp7ko_vs_wt_dn_g=parp7ko_vs_wt_dn$Symbol

library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(magrittr)

msigdbr_species()

mm_msigdb_df <- msigdbr(species = "Mus musculus")

head(mm_msigdb_df)

#Filter the human data frame to the KEGG pathways that are included in the
# curated gene sets
hs_GO_df <- mm_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C5", # This is to filter only to the C2 curated gene sets
    gs_subcat %in% c("GO:BP","GO:CC","GO:MF") # This is because we only want KEGG pathways
  )

hs_KEGG_df <- mm_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C2", # This is to filter only to the C2 curated gene sets
    gs_subcat %in% c("CP:KEGG") # This is because we only want KEGG pathways
  )


hs_Reactome_df <- mm_msigdb_df %>%
  dplyr::filter(
    gs_cat == "C2", # This is to filter only to the C2 curated gene sets
    gs_subcat %in% c("CP:REACTOME") # This is because we only want KEGG pathways
  )


parp7ko_vs_wt_sig_genes=c(parp7ko_vs_wt_dn_g,parp7ko_vs_wt_up_g)



## parp7ko Reactome enrichment
Reactome_ora_results <- enricher(
  gene = parp7ko_vs_wt_sig_genes, # A vector of your genes of interest
  pvalueCutoff = 0.05, # Can choose a FDR cutoff
  pAdjustMethod = "BH",
  universe = parp7ko_noGm_riken$Symbol,# Method to be used for multiple testing correction
  # The pathway information should be a data frame with a term name or
  # identifier and the gene identifiers
  TERM2GENE = dplyr::select(
    hs_Reactome_df,
    gs_name,
    gene_symbol
  )
)

View(Reactome_ora_results@result)
enrich_plot <- enrichplot::dotplot(GO_ora_results, showCategory=15,font.size=8,title="Reactome terms for signifcant DEG in PARP7KO vs WT for CR705 tumor",orderBy= "p.adjust", decreasing = FALSE)
enrich_plot


write.table(Reactome_ora_results@result,file="GO term enrichment for genes up in M2 relative to M0 for parp7ko terms.txt",col.names = T,row.names = T,sep="\t",quote = F)

DEG_list=subset(parp7ko_vs_wt,abs(parp7ko_vs_wt$log2FoldChange) > 1 & parp7ko_vs_wt$padj < 0.01)


library(grid)
pdf(file="Reactome_terms_for_significant_DEGs_in_CR705_PARP7KO_vs_WT.pdf",height = 8,width = 9)
grid.draw(enrich_plot)
dev.off()





