if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install() #To install latest version if multiple versions are installed
BiocManager::install(c("ggrepel", "matrixStats", "tidyverse", "pheatmap", "RColorBrewer","readxl", "EnhancedVolcano", "ggfortify", "magrittr"))

library(matrixStats)
library(tidyverse)
library(readxl)
library(EnhancedVolcano)
library(ggfortify)
library(magrittr)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)

setwd("/Users/xjacfr/Documents/Proteomics/4303/")

my_data <- read_xlsx("PCF_4303_220505.xlsx", sheet='raw', skip=2)
my_data_bak <- my_data #backup
uniprot <- read_xlsx("../Accession_uniprot-human_2021feb.xlsx")
uniprot <- uniprot[, 1:9]

# Join 'one gene' protein names and GO-term info to data by merging the dataframes
my_data <- left_join(my_data, 
                     uniprot, 
                     by="Accession")

# Fix missing gene names in my_data: 
i <- which(is.na(my_data$`Gene Symbol`) == TRUE)
my_data$`Gene Symbol`[i] <- my_data$`one Gene`[i]
# if no gene names in data file, just use the 'one gene' column from uniprot.
my_data$`Gene Symbol` <- my_data$`one Gene`
# check if missing values in 'one gene'
print("These protein entries have no gene names: ")
print(which(is.na(my_data$`one Gene`))) # For 4303, seven proteins are missing a gene name


# Check if any NAs left in Gene Symbol column, if so, fix manually:
needs_fix <- which(is.na(my_data$`Gene Symbol`) == TRUE)
if (length(needs_fix) > 0) {
  print("Manually change the gene names of entries with accession number:")
  print(my_data$Accession[needs_fix])
}
# Add relevant filter for GO-terms to use in plots
#go_terms <- "GO:0009617|GO:0009615|GO:0042742|GO:0050830|GO:0051607|GO:0019730|GO:0043312|GO:0090024|GO:0006954|GO:0006955|GO:0002682|GO:0045087|GO:0006953|GO:0019886"
#sign_genes_go <- my_data %>% 
#              filter(abs(my_data$log2FC_mean)>log2(1.2) & str_detect(my_data$`Gene ontology (biological process)`, go_terms))

# Select expression data columns
expr_data <- cbind(my_data[, 9:25], my_data[, 27:55]) # column 29 is the special sample to exclude
expr_data_wo23 <- expr_data[, -20] # removed sample 23

#Extract group labels and make index lists for t-test
group_labels <- str_extract(colnames(expr_data), "A[A-B]")
group_labels_wo23 <- str_extract(colnames(expr_data_wo23), "A[A-B]")
indices_gr1 <- grep("AA", group_labels) # For use in t-test below
indices_gr2 <- grep("AB", group_labels)
indices_gr2_wo23 <- indices_gr2[-3] # removed sample 23, the 3rd entry

# Count non-missing values per group 
my_data$non_missing_gr1 <- apply(expr_data, 1, function(x) sum(!is.na(x[indices_gr1])))
my_data$non_missing_gr2 <- apply(expr_data, 1, function(x) sum(!is.na(x[indices_gr2])))
my_data$non_missing_gr2_wo23 <- apply(expr_data, 1, function(x) sum(!is.na(x[indices_gr2_wo23])))
# Which proteins have more than 70% missing data?
too_few <- grep(TRUE, (my_data$non_missing_gr1 < 0.7*length(indices_gr1))|(my_data$non_missing_gr2 < 0.7*length(indices_gr2))) 
too_few_wo23 <- grep(TRUE, (my_data$non_missing_gr1 < 0.7*length(indices_gr1))|(my_data$non_missing_gr2_wo23 < 0.7*length(indices_gr2_wo23))) 
setdiff(too_few, too_few_wo23) # No 3122 cannot be used with all samples, but is ok if sample 23 is not used.

# Remove proteins with too many NAs for analysis
expr_data_wo23 <- expr_data_wo23[-too_few_wo23, ]
expr_data <- expr_data[-too_few_wo23, ] # Include protein 3122 to get same df rows

### Use expr_data_wo23 below:

# Group means, medians and fold changes    
my_data <- my_data[, -(which(colnames(my_data) == "23 AB"))]
my_data_clean <- my_data[-too_few_wo23, ]
my_data_clean$mean_gr1 <- rowMeans(expr_data_wo23[, indices_gr1], na.rm = TRUE)
my_data_clean$mean_gr2 <- rowMeans(expr_data_wo23[, -indices_gr1], na.rm = TRUE)

my_data_clean$FC_2vs1 <- my_data_clean$mean_gr2/my_data_clean$mean_gr1
my_data_clean$log2FC <- log2(my_data_clean$FC_2vs1)
my_data_clean$abs_log2_FC <- abs(my_data_clean$log2FC)

#my_data$median_AKI0 <- rowMedians(as.matrix(expr_data[, c(group_labels == "AKI0")]), na.rm = TRUE)
#my_data$median_AKI23 <- rowMedians(as.matrix(expr_data[, c(group_labels == "AKI23")]), na.rm = TRUE)
#my_data$FC_23_vs_0_med <- my_data$median_AKI23/my_data$median_AKI0
#my_data$log2FC_median <- log2(my_data$FC_23_vs_0_med)

#Plot data before t-test to see distributions, looks ok
hist(log2(t(expr_data[,indices_gr1])))
hist(log2(t(expr_data[,indices_gr2_wo23])))

# P-values from t-test
my_data_clean$p <- apply(expr_data, 1, function(x) {t.test(log2(x[indices_gr1]), log2(x[indices_gr2_wo23]), var.equal=TRUE)$p.value}) # log 2 to get normal distribution
my_data_clean$pwelch <- apply(expr_data, 1, function(x) {t.test(log2(x[indices_gr1]), log2(x[indices_gr2_wo23]))$p.value}) # log 2 to get normal distribution

# Plot p-values
hp<-hist(my_data_clean$p, breaks=30 ,plot=FALSE)
hpw<-hist(my_data_clean$pwelch, breaks=30 ,plot=FALSE)
cuts<-cut(hp$breaks, c(-Inf,0.04,Inf)) #Funkar inte med 0.05!!
plot(hp, col=c("red","white")[cuts], main="p-values, equal variance assumed", xlab="pvalue")
legend("topright","p < 0.05", fill="red")
plot(hpw, col=c("red","white")[cuts], main="p-values, welch correction", xlab="pvalue")
legend("topright","p < 0.05", fill="red")

my_data_clean$fdr <- p.adjust(my_data_clean$pwelch, "fdr")
h<-hist(my_data_clean$fdr, breaks=30 ,plot=FALSE)
cuts<-cut(h$breaks, c(-Inf,0.04,Inf)) #Funkar inte med 0.05!!
plot(h, col=c("red","white")[cuts], main="fdr distribution", xlab="fdr")
legend("topright","fdr < 0.05", fill="red")

fdr_list <- which((my_data_clean$fdr <= 0.05) == TRUE)
fdr <- length(fdr_list)
# Dataframe with proteins fdr ≤ 0.05:
#my_data_fdr <- my_data_clean[fdr_list, ]






#####----- LISTS OF SIGNIFICANTS -----#####

#signprot1_5 <- my_data_clean %>% 
#  filter((my_data_clean$FC_2vs1 > 1.5 | 1/my_data_clean$FC_2vs1 > 1.5) & my_data_clean$fdr <= 0.01) %>%
#  as.data.frame %>%
#  #na.omit(my_data$`Accession`) %>%
#  set_rownames(.$`Accession`)




#####----- PCA -----#####

#pdf("pca_plot_from_R.pdf") 

expr_data_wo23 %>%
  filter(my_data_clean$fdr<0.05 & my_data_clean$abs_log2_FC<abs(log2(1.5))) %>%
  drop_na %>%
  log2 %>% #always needed??
  t %>%
  prcomp(center=TRUE, scale=TRUE) %>%
  autoplot(data=data.frame(group=group_labels_wo23), col='group', size=2) +
  theme_classic() +
  geom_text_repel(aes(label = colnames(expr_data_wo23)), box.padding = unit(0.45, "lines"), max.overlaps = 15) +
  scale_color_manual(values=c("#D56115", "#0456A3")) +
  theme(aspect.ratio=1)

#dev.off() 


#####----- VOLCANO PLOT -----#####

# Colouring genes with certain GO-terms in red
# Labeling top 15 over expressed and under expressed, respectively
# Sort on absolute value of log2FC and use the top 30 entries. 

#Using mean values:

EnhancedVolcano(my_data_clean,                 
                x = 'log2FC', 
                y = 'fdr', 
                #lab = NA,
                lab = my_data_clean$`Gene Symbol`,
                #selectlab = my_data_clean %>% filter((abs_log2_FC)<abs(log2(1.5))) %>%
                #  as.data.frame %>% pull(`Gene Symbol`),
                xlim = c(-4, 4), 
                ylim= c(0, 6.5), 
                pCutoff = 0.01,
                FCcutoff = log2(3.0),
                labSize = 4, 
                title="Significant genes AA vs AB", 
                subtitle = "fdr < 0.05, FC > 1.5",
                col = c("grey30", "grey30", "grey30", "#6887E5"), #red2 for labeled points, #6887E5 lightblue
                legendPosition = NULL,
                caption = NULL,
                axisLabSize = 11,
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                arrowheads = FALSE,
                maxoverlapsConnectors = 30) +  
  guides(color = "none", size = "none") +
  scale_x_continuous(breaks=seq(-4, 4, by=0.5)) +
  scale_y_continuous(breaks=seq(0, 7)) +
  theme_classic() # +
  #geom_text_repel(aes(label = sign_prot_FC$`Gene Symbol`, col= sign_prot_FC$`Gene Symbol`), box.padding = unit(0.45, "lines")) +
  #geom_point(data = <>, aes(x = log2(FC_23_vs_0), -log10(p)), col="red2", alpha=1/1.2)
# In the line above, the data in geom_point function become red, the other above cutoffs become blue.

EnhancedVolcano(my_data_clean, 
                x = 'log2FC', 
                y = 'fdr', 
                lab = my_data_clean$`Gene Symbol`,
                selectLab = signprot1_5  %>% pull(`Gene Symbol`),
                xlim = c(-3, 3), 
                ylim= c(0, 6.5), 
                pCutoff = 0.05,
                FCcutoff = log2(1.5),
                labSize = 2.8, 
                title = NULL, 
                subtitle = NULL, 
                col = c("grey30", "grey30", "grey30", "#6887E5"), #red 2 for labeled points, #6887E5 lightblue
                legendPosition = NULL,
                caption = NULL,
                axisLabSize = 11,
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                arrowheads = FALSE,
                maxoverlapsConnectors = Inf) +  
  guides(color = "none", size = "none") +
  scale_x_continuous(breaks=seq(-3, 3, by=0.5)) +
  scale_y_continuous(breaks=seq(0, 7)) +
  theme_classic() # +
  #geom_point(data = my_data_clean %>% filter(Accession %in% all_srav_prot$`T: Accession`), aes(x = log2FC_high_vs_low, -log10(p)), col="red2", alpha=1/1.2)


EnhancedVolcano(my_data,                 
                x = 'log2FC_median', 
                y = 'p', 
                lab = my_data$`Gene Symbol`,
                #lab = NA,
                #selectLab = order(sign_genes_go$log2FC_mean, decreasing = TRUE) %>% pull(`Gene Symbol`),
                xlim = c(-4, 4), 
                ylim= c(0, 6.5), 
                pCutoff = 0.05,
                FCcutoff = log2(1.2),
                labSize = 2.8, 
                title = "Significant genes AKI23 vs AKI0, median values", 
                subtitle = "p < 0.05, FC > 1.2", 
                col = c("grey30", "grey30", "grey30", "#6887E5"), #red 2 for labeled points, #6887E5 lightblue
                legendPosition = NULL,
                caption = NULL,
                axisLabSize = 11,
                drawConnectors = TRUE,
                widthConnectors = 0.3,
                arrowheads = FALSE,
                maxoverlapsConnectors = Inf) +  
  guides(color = "none", size = "none") +
  scale_x_continuous(breaks=seq(-4, 4, by=0.5)) +
  scale_y_continuous(breaks=seq(0, 7)) +
  theme_classic() +
  geom_point(data = sign_prot_median, aes(x = log2(FC_23_vs_0_med), -log10(p)), col="red2", alpha=1/1.2)


#####----- HEATMAP -----#####


# What are breaks? /FJ
inner_breaks <-seq(-2, 2, length=98)
my_breaks <- c(-3, inner_breaks, 3)
my_col <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(inner_breaks))
my_col <- c(my_col[1], my_col, my_col[length(my_col)])

#go_term <- "GO:0043312|GO:0045087|GO:0006954"
#go_title <- "GO neutrofile degran, innate immune resp, inflam response, fdr<0.05"

# Select the most significant genes for the heatmap (p<0.05, log2FC>1.2)

#Using mean value difference:
my_data_clean %>% 
  filter(my_data_clean$fdr<0.05 & my_data_clean$abs_log2_FC<abs(log2(1.8))) %>%
  as.data.frame %>%
  #na.omit(my_data$`Accession`) %>%
  #set_rownames(.$`Gene Symbol`) %>%
  select(matches(c(" AA", " AB"))) %>% # Obs! Space at beginning!
  drop_na %>%
  log2 %>% 
  apply(., 1, function(x) {x[is.na(x) & group_labels_wo23 == "AA"] <- mean(x[group_labels_wo23 == "AA"], na.rm=TRUE)
  x[is.na(x) & group_labels_wo23 == "AB"] <- mean(x[group_labels_wo23 == "AB"], na.rm=TRUE) 
  x}) %>% # replace NA with mean
  t %>%
  pheatmap(scale = "row", 
           show_rownames = TRUE,
           #cluster_rows=FALSE,
           annotation_colors = list(group=c(`AA`="#9B3C8E", `AB`="#D56115")),
           angle_col = 90,
           main="Heatmap using FC = 1.8 as cutoff", 
           fontsize_row = 4.5, 
           fontsize_col = 7, 
           breaks=my_breaks,
           color = my_col
           )          

# Save dataframe as xlsx for enrichment analysis:

xlsx::write.xlsx(my_data_clean, file = "my_data_clean.xlsx")

        
           







