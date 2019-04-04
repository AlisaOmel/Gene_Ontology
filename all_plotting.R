library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(Rgraphviz)
library(topGO)
library(DOSE)
library(ReactomePA)
library(EnrichedHeatmap)
library(UpSetR)


df <- fortify(gene_data)

print("Setting up ColDEG sheet...")

require(clusterProfiler)
data(geneList, package="DOSE")
gene_data = read.csv("Desktop/ColonDEGs_new.csv")
geneList = gene_data[,2]
names(geneList) = as.character(gene_data[,1])
geneList = sort(geneList, decreasing = TRUE)
gene_list <- gene_data[,1]
#save <- write.csv(geneList, "Desktop/genelist.csv")

print("Running enrichGO Set...")
enrich_go_all <- enrichGO(gene = gene_list,
                      OrgDb= org.Mm.eg.db,
                      keyType = 'ENSEMBL',
                      ont= "ALL",
                      pAdjustMethod= "BH",
                      pvalueCutoff=0.05,
                      qvalueCutoff = 0.01,
                      readable = TRUE
                      )

enrich_go_summary_all <-as.data.frame(enrich_go_all)
print("Saving enrichGO Set...")
write.csv(enrich_go_summary_all, "Desktop/runs/enrichgo/enrich_go_summary_all.csv")

print("Making plots...")
pdf("Desktop/runs/enrichgo/cnet_plot_enrich_all.pdf", width= 15, height=20)
cnetplot(enrich_go_all)
dev.off()

pdf("Desktop/runs/enrichgo/cnet_plot_enrich_all_circular.pdf", width = 23, height=25)
cnetplot(enrich_go_all, circular= TRUE, colorEdge = TRUE)
dev.off()

pdf("Desktop/runs/enrichgo/dotplot_enrich.pdf")
dotplot(enrich_go_all)
dev.off()

pdf("Desktop/runs/enrichgo/heatplot.pdf", width = 45, height=13)
heatplot(enrich_go_all,showCategory = 50)
dev.off()

print("Running gse_go...")
#Run gse_go on all data sets:
gse_go <- gseGO(geneList, OrgDb=org.Mm.eg.db, keyType = "ENSEMBL",
                exponent = 1, nPerm = 10000, minGSSize = 10, maxGSSize = 500,
                pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = TRUE,
                seed = FALSE, by = "fgsea")

#head(summary(gse_go))
print("Changing IDs to Readable...")
#Change ID's from ENSEMBLE to SYMBOL:
name= setReadable(gse_go,'org.Mm.eg.db',"ENSEMBL")
gse_go_summary_all <-as.data.frame(name)

print("Saving gsego summary ...")
write.csv(gse_go_summary_all, "Desktop/runs/gsego/gse_go_summary_all.csv")

print("Making gse_go plots....")
#Plot gse_go 
pdf("Desktop/runs/gsego/cnet_plot_gse_all_circular.pdf", width = 23, height=25)
cnetplot(name, foldChange = geneList, circular= TRUE, colorEdge = TRUE)
dev.off()

pdf("Desktop/runs/gsego/cnet_plot_gse_all.pdf", width = 15, height = 20)
cnetplot(name, foldChange = geneList)
dev.off()

#emapplot(name, showCategory = 30, color = "p.adjust", layout = "kk")

pdf("Desktop/runs/gsego/dotplot_gse_all.pdf")
dotplot(gse_go, x = "GeneRatio", color = "p.adjust", showCategory = 50, split = NULL, font.size = 8)
dev.off()

pdf("Desktop/runs/gsego/heatplot_gse_all.pdf", width = 45, height=13)
heatplot(name, foldChange=geneList)
dev.off()

print("Setting up for KEGG analysis...")

gene_enterez <- read.csv("Desktop/Enterez_genes.csv")
id_list <- inner_join(gene_enterez, gene_data, by="Gene_Name") %>% tbl_df()
id_list <- as.data.frame(id_list)
id_list_enterez <- id_list[,-1]
id_list_enterez <- as.data.frame(id_list_enterez)
gene_names <- id_list_enterez[,1]

require(clusterProfiler)
data(geneList, package="DOSE")
geneList = id_list_enterez[,2]
names(geneList) = as.character(id_list_enterez[,1])
geneList = sort(geneList, decreasing = TRUE)

print("Running enrichKEGG analysis...")
enrichKEGG <- enrichKEGG(gene= gene_names,
              organism= 'mmu',
              pvalueCutoff = 0.05)

enrichKEGG_summary <- as.data.frame(enrichKEGG)

pdf("Desktop/runs/enrichKEGG/heatplot_KEGG.pdf", width = 45, height=13)
heatplot(enrichKEGG)
dev.off()


print("Saving enrichKEGG analysis")
write.csv(enrichKEGG_summary, "Desktop/runs/enrichKEGG/enrichKEGG_summary.csv")

print("DONE!")
