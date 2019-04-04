library(biomaRt)
library(dplyr)
library(magicfor)

column_names <- c("gene_id", "go_pathway", "interaction_type", "process")
gene_interactions <- read.csv("Desktop/annotations_10090.tsv", sep = "\t", header= FALSE) %>% tbl_df()
colnames(gene_interactions) <- column_names
gene_data <- read.csv("Desktop/cdx_related_genes_only.csv")
Col_Deg <- read.csv("Desktop/ColonDEGs_symbol.csv") %>% tbl_df()
go_id <- read.csv("Desktop/go_id.csv") %>% tbl_df()
write.csv(full_anotation_symbol, "Desktop/full_annotation_symbol.tsv")

out <- vector("list", length(go_id$go_pathway))
for (i in go_id$go_pathway) {
  out[[i]] <- full_anotation_symbol %>%
    filter(go_pathway == i)
}
  
a = "GO:0010558"
sim_1 <- full_anotation_symbol %>%
  filter(go_pathway == a)

LC_1 <- Col_Deg %>%
  semi_join(sim_1, by = "mgi_symbol")

mart = useMart("ensembl",dataset="mmusculus_gene_ensembl")

full_anotation <- separate(gene_interactions,col='gene_id', into= c('id','gene'), sep='10090.')
full_anotation <- full_anotation[,-1]
full_anotation <- full_anotation[,-4]
genes <- full_anotation$gene
G_list <- getBM(filters='ensembl_peptide_id', attributes= c ('ensembl_peptide_id','mgi_symbol'),values=genes,mart= mart)
full_anotation_symbol <- merge( full_anotation, G_list,by.x='gene', by.y='ensembl_peptide_id')
write.csv( full_anotation_symbol,'Desktop/Cdx2_gene_lists/full_anotation.csv')

full_annotation_symbol <- read.csv('Desktop/full_anotation.csv') %>% tbl_df()



lancl2_gopathways <- full_anotation_symbol %>%
  filter(mgi_symbol =="Lancl2")

write.csv(lancl2_gopathways, "Desktop/Lancl2_gopathways.csv")
cdx2_gopathways <- full_anotation_symbol %>%
  filter(mgi_symbol == "Cdx2")
write.csv(cdx2_gopathways, "Desktop/cdx2_gopathways.csv")

lancl2_1 <- lancl2_gopathways %>%
  inner_join(cdx2_gopathways, by = 'go_pathway') %>%
  select(go_pathway)

glimpse(lancl2_1)


GO0045892_genes  <- read.csv("Desktop/GO0045892_genes.csv") %>% tbl_df()
common_genes <- Col_Deg %>%
  filter(ID %in% GO0045892_genes$ID) %>%
  
glimpse(common_genes)

common_genes <- Col_Deg %>%
  inner_join(GO0045892_genes, by = "ID")
