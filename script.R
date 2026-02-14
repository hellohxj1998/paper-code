# ==============================================================================
# Setup and Data Loading
# ==============================================================================
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggsci)
library(RColorBrewer)

# Source CIBERSORT function
source("Cibersort.R")

# File paths
LM22_file   <- "LM22.txt"
exp_file    <- "allCOAD_surv_risk_clinc_expr.txt"
group_file  <- "group.csv"

# Read and clean expression data
expr_df <- read.delim(exp_file, sep = "\t", header = TRUE, 
                      row.names = 1, na.strings = "#NAME?")
write.table(expr_df, "exprSet_clean.txt", sep = "\t", quote = FALSE)

# ==============================================================================
# 1. CIBERSORT Analysis
# ==============================================================================
results <- CIBERSORT(LM22_file, "exprSet_clean.txt", perm = 100, QN = FALSE)
results_df <- as.data.frame(results) %>% rownames_to_column("Sample")

# Save raw results
write.table(results_df, "CIBERSORT_Results.txt", sep = "\t", 
            quote = FALSE, row.names = FALSE)

# ==============================================================================
# 2. Visualization: Immune Composition
# ==============================================================================
# Prepare data for plotting
group_info <- read.csv(group_file)
colnames(group_info) <- c("Sample", "Group")

plot_data <- results_df %>%
  select(1:23) %>% # Keep Sample + 22 Cell types
  left_join(group_info, by = "Sample") %>%
  pivot_longer(cols = 2:23, names_to = "Celltype", values_to = "Composition")

# Boxplot with Statistical Test
p1 <- ggplot(plot_data, aes(x = Celltype, y = Composition, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test", hide.ns = TRUE) +
  scale_fill_npg() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Immune Cell Infiltration Difference", x = "", y = "Relative Proportion")

ggsave("CIBERSORT_Boxplot.pdf", p1, width = 12, height = 7)

# Stacking Barplot
colors <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(4, "Set3"))
p2 <- ggplot(plot_data, aes(x = Sample, y = Composition, fill = Celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Cell Proportion per Sample", x = "Samples", y = "Proportion")

ggsave("CIBERSORT_Barplot.pdf", p2, width = 10, height = 6)

# ==============================================================================
# 3. Correlation: Genes vs. Immune Cells
# ==============================================================================
target_genes <- c("PYGL", "MMP3", "FABP4")
exp_matrix <- as.matrix(expr_df)

# Validate genes
matched_genes <- intersect(rownames(exp_matrix), target_genes)
imm_matrix <- results_df %>% column_to_rownames("Sample") %>% select(1:22)

cor_results <- list()

for (gene in matched_genes) {
  gene_vec <- exp_matrix[gene, ]
  common <- intersect(names(gene_vec), rownames(imm_matrix))
  
  gene_sub <- gene_vec[common]
  imm_sub <- imm_matrix[common, ]
  
  res <- lapply(colnames(imm_sub), function(cell) {
    test <- cor.test(gene_sub, imm_sub[[cell]], method = "pearson")
    data.frame(Gene = gene, Celltype = cell, Cor = test$estimate, Pval = test$p.value)
  })
  cor_results[[gene]] <- bind_rows(res)
}

final_cor <- bind_rows(cor_results) %>%
  mutate(Sign = case_when(Pval < 0.001 ~ "***", Pval < 0.01 ~ "**", Pval < 0.05 ~ "*", TRUE ~ ""))

# Correlation Heatmap
p3 <- ggplot(final_cor, aes(x = Gene, y = Celltype)) +
  geom_tile(aes(fill = Cor), color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1,1)) +
  geom_text(aes(label = Sign), color = "black", size = 5) +
  theme_minimal() +
  labs(title = "Gene-Immune Correlation", fill = "Pearson R")

ggsave("MultiGene_Immune_Correlation.pdf", p3, width = 8, height = 10)
write.csv(final_cor, "Correlation_Results_Full.csv", row.names = FALSE)

message("Done. All plots and tables generated.")