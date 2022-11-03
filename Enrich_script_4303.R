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

my_data_wo23_clean <- read_xlsx("my_data_wo23_clean.xlsx")
my_data_wo23_clean <- my_data_wo23_clean[, -1]

# Filter data on fdr and FC, remove NAs

sig_data_wo23 <- my_data_wo23_clean %>% dplyr::filter(fdr<0.05) %>% dplyr::filter(pwelch<0.05) %>% dplyr::filter(abs_log2_FC>abs(log2(1.8)))
# sig_data_wo23 <- dplyr::select(sig_data_wo23$Accession) %>% na.omit
sig_data_na <- which((sig_data_wo23$Accession == "NA") == TRUE) #no NAs
#sig_data_wo23 <- sig_data_wo23[-sig_data_na,]
  
# Add columns with entrezid and ensembl:

sig_data_wo23_ids <- bitr(sig_data_wo23[[1]], fromType="UNIPROT", toType=c("ENTREZID", "ENSEMBL"), OrgDb="org.Hs.eg.db")
# "Warning message:
# In bitr(sig_data_all[[1]], fromType = "UNIPROT", toType = c("ENTREZID",  :
# 2.18% of input gene IDs are fail to map..."

# Remove duplicates of Uniprot
sig_data_wo23_ids <- sig_data_wo23_ids[!duplicated(sig_data_wo23_ids$UNIPROT),] # Use as input in enrichment below

gsub(" ", "", sig_data_wo23_ids$ENTREZID)
react <- enrichPathway(gene=sig_data_wo23_ids$ENTREZID, organism="human", pvalueCutoff=0.05, qvalueCutoff = 0.05, readable=TRUE)

# Plot
dotplot(react, showCategory = 20, font.size=9, label_format=45, title="Reactome pathway analysis excluding sample 23, p < 0.05, q < 0.05, FC > 1.8")

### Do the same pathway analysis for data including sample 23:

my_data_all_clean <- read_xlsx("my_data_all_clean.xlsx")
my_data_all_clean <- my_data_all_clean[, -1]
row.names(my_data_all_clean) <- my_data_all_clean$Accession # For getting info from specific rows later

sig_data_all <- my_data_all_clean %>% dplyr::filter(fdr<0.05) %>% dplyr::filter(pwelch<0.05) %>% dplyr::filter(abs_log2_FC>abs(log2(1.8)))
# sig_data_all <- dplyr::select(sig_data_all$Accession) %>% na.omit # not working!
sig_data_na <- which((sig_data_all$Accession == "NA") == TRUE) #no NAs
#sig_data_all <- sig_data_all[-sig_data_na,]

# Add columns with entrezid and ensembl:

sig_data_all_ids <- bitr(sig_data_all[[1]], fromType="UNIPROT", toType=c("ENTREZID", "ENSEMBL"), OrgDb="org.Hs.eg.db")
# Paste warning message


# Remove duplicates of Uniprot
sig_data_all_ids <- sig_data_all_ids[!duplicated(sig_data_all_ids$UNIPROT),]

gsub(" ", "", sig_data_all_ids$ENTREZID)
react <- enrichPathway(gene=sig_data_all_ids$ENTREZID, organism="human", pvalueCutoff=0.05, qvalueCutoff = 0.05, readable=TRUE)

# Plot
dotplot(react, showCategory = 20, font.size=9, title="Reactome pathway analysis including sample 23, FC = 1.8")



##### GO-terms, grouping and enrichment analysis, plus plots. (Excluding sample 23!)

# grouping example
#ggo <- groupGO(gene     = sig_data_wo23_ids[[2]],
#               OrgDb    = org.Hs.eg.db,
#               ont      = "BP",              # Biological processes
#                level    = 3,
#                readable = TRUE)

# enrichment, creating the enrichResult object
ego_bp <- enrichGO(gene          = sig_data_wo23_ids[[2]],
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# Calculate richFactor for each GO-term (divide nominators of GeneRatio and BgRatio)
ego_bp <- mutate(ego_bp, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) # \\d+ matches any number of digits

# Plot enrichment results
ggplot(ego_bp, showCategory = 20, 
        aes(richFactor, fct_reorder(Description, richFactor))) + 
        geom_segment(aes(xend=0, yend = Description)) +
        geom_point(aes(color=p.adjust, size = Count)) +
        scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
        scale_size_continuous(range=c(2, 10)) +
        theme_minimal() + 
        theme(axis.text = element_text(size = 12)) +
        xlab("over-representation") + # This is the richFactor
        ylab(NULL) + 
        ggtitle("Enriched Gene Ontology - Biological Process")


# Do same with MF (Molecular Function):

ego_mf <- enrichGO(gene       = sig_data_wo23_ids[[2]],
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
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
        theme(axis.text = element_text(size = 11)) +
        scale_y_discrete(labels = scales::label_wrap(45)) + # To divide long rows into multiple rows with a certain width
        xlab("over-representation") + # This is the richFactor
        ylab(NULL) + 
        ggtitle("Enriched Gene Ontology - Molecular Function")


################## GO-terms with one gene per row ###############

# Make dataframe with GO-terms and related genes (mf)

genesets_mf <- as.data.frame(cbind(ego_mf@result$ID, ego_mf@result$Description, ego_mf@result$geneID, ego_mf@result$richFactor, ego_mf@result$p.adjust))
colnames(genesets_mf)[1:5] <- c("GO-term", "Description", "Gene Symbol", "richFactor", "p.adjust")
genesets_mf <- genesets_mf %>% separate_rows("Gene Symbol", sep = "/") # One gene on each row

# Add proteomics and statistics info:
#genesets_mf <- left_join(genesets_mf, 
#                     my_data_wo23_clean, 
#                     by="Gene Symbol")
#genesets_mf_clean <- genesets_mf[,-(11:72)] #keep only columns with proteomics data and statistics
#> colnames(genesets_mf_clean)
#[1] "GO-term"                            "Gene Symbol"                        "Accession"                         
#[4] "Description"                        "Coverage [%]"                       "# Peptides"                        
#[7] "# PSMs"                             "# Unique Peptides"                  "MW [kDa]"                          
#[10] "Score Mascot: Mascot"               "Gene ontology (biological process)" "Gene ontology (cellular component)"
#[13] "Gene ontology (molecular function)" "one Gene"                           "non_missing_gr1"                   
#[16] "non_missing_gr2"                    "non_missing_gr2_wo23"               "mean_gr1"                          
#[19] "mean_gr2"                           "FC_2vs1"                            "log2FC"                            
#[22] "abs_log2_FC"                        "p"                                  "pwelch"                            
#[25] "fdr" 

### Do the same with biological processes (bp)

genesets_bp <- as.data.frame(cbind(ego_bp@result$ID, ego_bp@result$Description, ego_bp@result$geneID, ego_bp@result$richFactor, ego_bp@result$p.adjust))
colnames(genesets_bp)[1:5] <- c("GO-term", "Description", "Gene Symbol", "richFactor", "p.adjust")
genesets_bp <- genesets_bp %>% separate_rows("Gene Symbol", sep = "/") # One gene on each row

# Add proteomics and statistics info:
#genesets_bp <- left_join(genesets_bp, 
#                         my_data_wo23_clean, 
#                         by="Gene Symbol")
#genesets_bp_clean <- genesets_bp[,-(11:72)] #keep only columns with proteomics data and statistics


# Add proteomics data for each gene
#peptide_info_etc <- dplyr::select(my_data_clean, 1:6, GeneID="Gene Symbol", pwelch, fdr, abs_log2_FC) %>% na.omit # There is one NA in peptide_info but not in genesets_mf. (Nr 1641)
#all_info <- as.data.frame(matrix(NA, nrow = 3862, ncol = 11))
#all_info <- inner_join(genesets_mf, peptide_info_etc, by = 'GeneID')

xlsx::write.xlsx(genesets_mf, file = "enrich_mf.xlsx")

### Do the same with biological processes (bp)
#peptide_info_etc <- dplyr::select(my_data_clean, 1:6, GeneID="Gene Symbol", pwelch, fdr, abs_log2_FC) %>% na.omit # There is one NA in peptide_info but not in genesets_mf. (Nr 1641)
#all_info <- as.data.frame(matrix(NA, nrow = 3862, ncol = 11))
#all_info <- inner_join(genesets_bp, peptide_info_etc, by = 'GeneID')

xlsx::write.xlsx(genesets_bp, file = "enrich_bp.xlsx")


######## SOME ANALYSIS ON THE NUMBER OF PROTEINS IN OBJECTS ########

# to better understand the enrichGO function

genes_GO <- str_split(ego_bp@result$geneID, "/")
genes_GO <- unlist(genes_GO, recursive = FALSE)
genes_GO <- genes_GO[!duplicated(genes_GO)]
setdiff(ego_bp@gene2Symbol, genes_GO) # 53 proteins

original_names_notGO <- setdiff(sig_data_wo23$`Gene Symbol`, genes_GO) # 71 proteins
enrich_names_notGO <- setdiff(ego_bp@gene2Symbol, genes_GO) # 53
setdiff(original_names_notGO, enrich_names_notGO) # 19 proteins in original that are not in enriched. These have changed names
#[1] "MYL12B"   "USP11"    "H4C1"     "MT-CO2"   "GK3P"     NA         "IGHM"     "POM121C"  "FAU"      "CHCHD2P9"
#[11] "MT-ND4"   "FAM160B1" "TUBA4B"   "C11orf98" "MT-CYB"   "IGHG3"    "IGKC"     "GPR89A"   "IGLC2"  
setdiff(enrich_names_notGO, original_names_notGO) # 1 protein in gene2symbol that is not used for enrichment. No GO-terms?
#[1] "FHIP2A"



######## Save enrichment GO-terms per gene in list: ########

# Make list of GO-term lists for each gene
GO_bp <- list()
for (name in ego_bp@gene2Symbol) {
  indices <- which(str_detect(ego_bp@result$geneID, name) == TRUE)
  GO_bp[length(GO_bp)+1] <- list(ego_bp@result$ID[indices])
}

# How many genes do not map to GO-terms at all?
sum(lapply(GO_bp, length) == 0) # 52 genes (538-52 = 486...)

#Filter GO-term results, p.adjust < 0.05
ego_bp_filtered <- ego_bp %>% dplyr::filter(ego_bp@result$p.adjust < 0.05)

GO_bp_filtered <- list()
for (name in ego_bp_filtered@gene2Symbol) {
  indices <- which(str_detect(ego_bp@result$geneID, name) == TRUE)
  GO_bp[length(GO_bp)+1] <- list(ego_bp@result$ID[indices])
}
