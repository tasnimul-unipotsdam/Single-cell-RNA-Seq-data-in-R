# script to visualize gene expression data
# Updated with generic variable names and automated statistics

# Load libraries
library(tidyverse)
library(ggplot2)
library(ggpubr) # For statistical brackets

# 1. Data Loading ---------------------------------------------------------
# exp_data to be used (generated from previous long-format processing)
# exp_data <- read.delim('../data/GSE183947_long_format.txt', header = T)


# 2. Barplot: Expression per Sample ---------------------------------------
exp_data %>%
  filter(gene == 'BRCA1') %>%
  ggplot(aes(x = samples, y = FPKM, fill = tissue)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "BRCA1 Expression per Sample", y = "Expression (FPKM)")


# 3. Density Plot: Distribution by Tissue ---------------------------------
exp_data %>%
  filter(gene == 'BRCA1') %>%
  ggplot(aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = 0.4) +
  theme_classic() +
  labs(title = "Distribution of BRCA1 Expression")


# 4. Boxplot & Violin: Statistical Significance ---------------------------
# Define comparisons for the brackets
# Ensure these names match your data exactly (e.g., "yes" vs "no")
stat_comparisons <- list( c("yes", "no") )

exp_data %>%
  filter(gene == 'BRCA1') %>%
  ggplot(aes(x = metastasis, y = FPKM, fill = metastasis)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # Add p-values and brackets
  stat_compare_means(comparisons = stat_comparisons, 
                     method = "t.test", 
                     label = "p.signif") + 
  theme_bw() +
  labs(title = "BRCA1 Expression: Metastasis Comparison",
       x = "Metastasis Status",
       y = "FPKM")



# 5. Scatterplot: Gene Correlation ----------------------------------------
exp_data %>%
  filter(gene %in% c('BRCA1', 'BRCA2')) %>%
  # Updated spread() to pivot_wider()
  pivot_wider(names_from = gene, values_from = FPKM) %>%
  ggplot(aes(x = BRCA1, y = BRCA2, color = tissue)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_light() +
  labs(title = "Correlation: BRCA1 vs BRCA2")


# 6. Heatmap: Multiple Genes of Interest ----------------------------------
target_genes <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')

# Save heatmap to PDF
pdf("gene_expression_heatmap_updated.pdf", width = 10, height = 8)

exp_data %>%
  filter(gene %in% target_genes) %>%
  ggplot(aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Gene Expression Heatmap", x = "Samples", y = "Gene")

dev.off()
