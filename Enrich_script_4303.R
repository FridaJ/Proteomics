# Script for proteomics project 4303, the enrichment part using both ReactomePA and enrichGO

if (!requireNamespace("BiocManager", quietly = TRUE))
  + install.packages("BiocManager")
#BiocManager::install() #To install latest version if multiple versions are installed
BiocManager::install(c("tidyr", "xlsx", "ReactomePA", "clusterProfiler", "readxl", "dplyr", "org.Hs.eg.db", "ggplots2", "forcats", "enrichplot"))

library(ReactomePA) 
library(clusterProfiler)
library(readxl)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(forcats)
library(enrichplot)
library(xlsx)
library(tidyr)


setwd("/Users/xjacfr/Documents/Proteomics/4303/")

my_data_clean <- read_xlsx("my_data_clean.xlsx")
my_data_clean <- my_data_clean[, -1]
row.names(my_data_clean) <- my_data_clean$Accession # For getting info from specific rows later

sig_data <- my_data_clean %>% dplyr::filter(fdr<0.05) 
#sig_data_uni <- dplyr::select('Accession') %>% na.omit
sig_data_uni <- sig_data$Accession %>% na.omit

# Add columns with entrezid and ensembl:

sig_data_ids <- bitr(sig_data[[1]], fromType="UNIPROT", toType=c("ENTREZID", "ENSEMBL"), OrgDb="org.Hs.eg.db")
# Warning message:
# In bitr(sig_data[[1]], fromType = "UNIPROT", toType = c("ENTREZID",  :
# 1.4% of input gene IDs are fail to map...

# Set row names to Uniprot IDs to easily map them to Entrez IDs:
#row.names(sig_data_ids) <- sig_data_ids$UNIPROT

# Remove duplicates of Uniprot
sig_data_ids <- sig_data_ids[!duplicated(sig_data_ids$UNIPROT),]

gsub(" ", "", sig_data_ids$ENTREZID)
react <- enrichPathway(gene=sig_data_ids$ENTREZID, organism="human", pvalueCutoff=0.05, readable=TRUE)

# Plot
dotplot(react, showCategory = 20, font.size=11, title="Reactome pathway analysis")


##### GO-terms, grouping and enrichment analysis, plus plots

# grouping
ggo <- groupGO(gene     = sig_data_ids[[2]],
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",              # Biological processes
               level    = 3,
               readable = TRUE)

# enrichment
ego_bp <- enrichGO(gene          = sig_data_ids[[2]],
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# Calculate richFactor for each GO-term (divide nominators of GeneRatio and BgRatio)
ego_bp <- mutate(ego_bp, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) # \\d+ matches any number of digits

# Plot
ggplot(ego_bp, showCategory = 20, 
        aes(richFactor, fct_reorder(Description, richFactor))) + 
        geom_segment(aes(xend=0, yend = Description)) +
        geom_point(aes(color=p.adjust, size = Count)) +
        scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
        scale_size_continuous(range=c(2, 10)) +
        theme_minimal() + 
        theme(axis.text = element_text(size = 15)) +
        xlab("over-representation") +
        ylab(NULL) + 
        ggtitle("Enriched Gene Ontology - Biological Process")


# Do same with MF (Molecular Function):

ego_mf <- enrichGO(gene       = sig_data_ids[[2]],
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# Calculate richFactor for each GO-term (divide nominators of GeneRatio and BgRatio)
ego_mf <- mutate(ego_mf, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) # \\d+ matches any number of digits

ggplot(ego_mf, showCategory = 20, 
        aes(richFactor, fct_reorder(Description, richFactor))) + 
        geom_segment(aes(xend=0, yend = Description)) +
        geom_point(aes(color=p.adjust, size = Count)) +
        scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
        scale_size_continuous(range=c(2, 10)) +
        theme_minimal() + 
        theme(axis.text = element_text(size = 15)) +
        xlab("over-representation") +
        ylab(NULL) + 
        ggtitle("Enriched Gene Ontology - Molecular Function")


##################

# Make dataframe with GO-terms and related genes

genesets_mf <- as.data.frame(cbind(ego_mf@result$ID, ego_mf@result$geneID))
colnames(genesets_mf)[1] <- "GO-term"
colnames(genesets_mf)[2] <- "GeneID"
genesets_mf <- genesets_mf %>% separate_rows(GeneID, sep = "/") # One gene on each row

# Add proteomics data for each gene

peptide_info_etc <- dplyr::select(my_data_clean, 1:6, GeneID="Gene Symbol", pwelch, fdr, abs_log2_FC) %>% na.omit # There is one NA in peptide_info but not in genesets_mf. (Nr 1641)
all_info <- as.data.frame(matrix(NA, nrow = 3862, ncol = 11))

all_info <- inner_join(genesets_mf, peptide_info_etc, by = 'GeneID')

xlsx::write.xlsx(all_info, file = "enrich_mf.xlsx")
