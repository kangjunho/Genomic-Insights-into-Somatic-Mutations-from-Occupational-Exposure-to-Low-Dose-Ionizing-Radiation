# Visualization Script for Mutation Analysis

# Load Required Libraries
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggbreak)
library(tidyr)

# ---- Figure 2: Variant Distribution ----

# 1. Variant Classification Frequency (Coding)
variant_counts_coding <- df_coding %>% count(Variant_Classification, sort = TRUE) %>% filter(n > 0)
variant_df_coding <- as.data.frame(variant_counts_coding)
p_variant_classification_coding <- ggplot(variant_df_coding, aes(x = reorder(Variant_Classification, n), y = n, fill = Variant_Classification)) +
  geom_bar(stat = "identity") + coord_flip() + theme_minimal() +
  theme(text = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 16), axis.text = element_text(size = 14)) +
  labs(title = "Variant Classification (Coding)", x = "Classification", y = "Count") +
  scale_fill_brewer(palette = "Set3", guide = FALSE)

# 2. Variant Classification Frequency (Non-coding)
variant_counts_noncoding <- df_noncoding %>% count(Variant_Classification, sort = TRUE) %>% filter(n > 0)
variant_df_noncoding <- as.data.frame(variant_counts_noncoding)
p_variant_classification_noncoding <- ggplot(variant_df_noncoding, aes(x = reorder(Variant_Classification, n), y = n, fill = Variant_Classification)) +
  geom_bar(stat = "identity") + coord_flip() + theme_minimal() +
  theme(text = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 16), axis.text = element_text(size = 14)) +
  labs(title = "Variant Classification (Non-coding)", x = "Classification", y = "Count") +
  scale_fill_brewer(palette = "Set3", guide = FALSE)

# 3. Variants per Sample (Coding)
variants_per_sample_coding <- df_coding %>% group_by(Tumor_Sample_Barcode, Variant_Classification) %>% summarise(count = n()) %>% ungroup()
sample_totals_coding <- variants_per_sample_coding %>% group_by(Tumor_Sample_Barcode) %>% summarise(total_count = sum(count)) %>% arrange(desc(total_count))
variants_per_sample_coding <- variants_per_sample_coding %>% left_join(sample_totals_coding, by = "Tumor_Sample_Barcode")
median_value_coding <- median(sample_totals_coding$total_count)
p_variants_per_sample_coding <- ggplot(variants_per_sample_coding, aes(x = reorder(Tumor_Sample_Barcode, -total_count), y = count, fill = Variant_Classification)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = median_value_coding, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = median_value_coding + 1, label = paste("Median:", median_value_coding), hjust = 0, size = 5, fontface = "bold") +
  theme_minimal() +
  theme(text = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Variants per Sample (Coding)", x = "Sample", y = "Variant Count") +
  scale_fill_brewer(palette = "Set3", guide = FALSE)

# 4. Variants per Sample (Non-coding)
variants_per_sample_noncoding <- df_noncoding %>% group_by(Tumor_Sample_Barcode, Variant_Classification) %>% summarise(count = n()) %>% ungroup()
sample_totals_noncoding <- variants_per_sample_noncoding %>% group_by(Tumor_Sample_Barcode) %>% summarise(total_count = sum(count)) %>% arrange(desc(total_count))
variants_per_sample_noncoding <- variants_per_sample_noncoding %>% left_join(sample_totals_noncoding, by = "Tumor_Sample_Barcode")
median_value_noncoding <- median(sample_totals_noncoding$total_count)
p_variants_per_sample_noncoding <- ggplot(variants_per_sample_noncoding, aes(x = reorder(Tumor_Sample_Barcode, -total_count), y = count, fill = Variant_Classification)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = median_value_noncoding, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = median_value_noncoding + 1, label = paste("Median:", median_value_noncoding), hjust = 0, size = 5, fontface = "bold") +
  theme_minimal() +
  theme(text = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Variants per Sample (Non-coding)", x = "Sample", y = "Variant Count") +
  scale_fill_brewer(palette = "Set3", guide = FALSE)

# 5-6. Variant Type Frequency (Coding & Non-coding)
plot_variant_type <- function(df, title) {
  variant_type_counts <- df %>% count(Variant_Type, sort = TRUE)
  ggplot(variant_type_counts, aes(x = reorder(Variant_Type, n), y = n, fill = Variant_Type)) +
    geom_bar(stat = "identity") + coord_flip() + theme_minimal() +
    theme(text = element_text(face = "bold", size = 14),
          axis.title = element_text(size = 16), axis.text = element_text(size = 14)) +
    labs(title = title, x = "Variant Type", y = "Count") +
    scale_fill_brewer(palette = "Pastel1", guide = FALSE)
}
p_variant_type_coding <- plot_variant_type(df_coding, "Variant Type (Coding)")
p_variant_type_noncoding <- plot_variant_type(df_noncoding, "Variant Type (Non-coding)")

# 7-8. SNV Class Frequency (Top 10)
plot_snv_class <- function(df, title) {
  snv_class_counts <- df %>%
    mutate(SNV_Class = paste0(Tumor_Seq_Allele1, ">", Tumor_Seq_Allele2)) %>%
    count(SNV_Class, sort = TRUE) %>%
    filter(SNV_Class != "__UNKNOWN__") %>%
    top_n(10, n)
  ggplot(snv_class_counts, aes(x = reorder(SNV_Class, n), y = n, fill = SNV_Class)) +
    geom_bar(stat = "identity") + coord_flip() + theme_minimal() +
    theme(text = element_text(face = "bold", size = 14),
          axis.title = element_text(size = 16), axis.text = element_text(size = 14)) +
    labs(title = title, x = "SNV Class", y = "Count") +
    geom_text(aes(label = n), hjust = 0.7, size = 5, fontface = "bold") +
    scale_fill_brewer(palette = "Set3")
}
p_snv_class_coding <- plot_snv_class(df_coding, "Top 10 SNV Class (Coding)")
p_snv_class_noncoding <- plot_snv_class(df_noncoding, "Top 10 SNV Class (Non-coding)")

# Combine Plots
combined_plot_coding <- (p_variants_per_sample_coding | p_variant_classification_coding) / (p_variant_type_coding | p_snv_class_coding)
combined_plot_noncoding <- (p_variants_per_sample_noncoding | p_variant_classification_noncoding) / (p_variant_type_noncoding | p_snv_class_noncoding)
combined_plot <- (combined_plot_coding | combined_plot_noncoding) +
  plot_annotation(title = '', theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))

# Save Figure 2
ggsave("Summary_plot.png", plot = combined_plot, width = 48, height = 24, dpi = 600)

# Figure 3 Visualization
# Load Required Libraries
library(ggplot2)
library(ggpubr)
library(ggbreak)
library(dplyr)
library(tidyr)

# Thresholds for y-axis breaks
y_threshold_low_age <- 40
y_threshold_high_age <- 1500
y_threshold_low_coding <- 10
y_threshold_high_coding <- 200
y_threshold_low_noncoding <- 30
y_threshold_high_noncoding <- 1300
y_threshold_low_total <- 50
y_threshold_high_total <- 1300

# 1. Age vs Total Variants
plot_age <- ggplot(test_df, aes(x = age2, y = Total)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "navy", fill = "lightblue", se = TRUE) +
  stat_cor(method = "pearson", label.x = 20, label.y = 1950, color = "darkred", size = 6, fontface = "bold.italic") +
  theme_minimal() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black")) +
  labs(title = "Age", x = "Age", y = "Number of total variants") +
  scale_x_continuous(limits = c(19.9, 51.15)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_y_break(c(y_threshold_low_age, y_threshold_high_age), scales = 1.5, space = 0.3)

# 2. Dose vs Coding Frequency
plot_Coding <- ggplot(test_df, aes(x = Dose_T, y = Coding_Frequency)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "navy", fill = "lightblue", se = TRUE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 310, color = "darkred", size = 6, fontface = "bold.italic") +
  theme_minimal() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black")) +
  labs(title = "Coding variants", x = "Cumulative dose", y = "Number of coding variants") +
  scale_x_continuous(limits = c(0, 85)) +
  scale_y_break(c(y_threshold_low_coding, y_threshold_high_coding), scales = 1.5, space = 0.3)

# 3. Dose vs Noncoding Frequency
plot_Noncoding <- ggplot(test_df, aes(x = Dose_T, y = Noncoding_Frequency)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "navy", fill = "lightblue", se = TRUE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 1750, color = "darkred", size = 6, fontface = "bold.italic") +
  theme_minimal() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black")) +
  labs(title = "Noncoding variants", x = "Cumulative dose", y = "Number of noncoding variants") +
  scale_x_continuous(limits = c(0, 85)) +
  scale_y_break(c(y_threshold_low_noncoding, y_threshold_high_noncoding), scales = 1.5, space = 0.3)

# 4. Dose vs Total Variants
plot_Total <- ggplot(test_df, aes(x = Dose_T, y = Total)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "navy", fill = "lightblue", se = TRUE) +
  stat_cor(method = "pearson", label.x = 0, label.y = 2100, color = "darkred", size = 6, fontface = "bold.italic") +
  theme_minimal() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black")) +
  labs(title = "Total variants", x = "Cumulative dose", y = "Number of Total variants") +
  scale_x_continuous(limits = c(0, 85)) +
  scale_y_break(c(y_threshold_low_total, y_threshold_high_total), scales = 1.5, space = 0.3)

# 5. Clinical Variable Boxplot
pre_Correlation_df$Smoking <- ifelse(pre_Correlation_df$Smoking == "Smoking", 1, 0)
pre_Correlation_df$Drinking <- ifelse(pre_Correlation_df$Drinking == "Drinking", 1, 0)
pre_Correlation_df <- pre_Correlation_df[, c("Total", "Smoking", "Drinking")]

long_df <- pivot_longer(pre_Correlation_df, cols = -Total, names_to = "Variable", values_to = "Value")
p_values <- long_df %>% group_by(Variable) %>% summarize(p_value = t.test(Total ~ Value)$p.value)
long_df <- left_join(long_df, p_values, by = "Variable") %>% mutate(P_value_label = paste("p =", format(p_value, digits = 2)))

clinical_box <- ggplot(long_df, aes(x = as.factor(Value), y = Total, fill = Variable)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 2, aes(color = Variable)) +
  facet_wrap(~ Variable, scales = "free_x") +
  theme_minimal() +
  labs(y = "Total", x = "Value", title = "Comparison of Total Variants by Variables") +
  theme(strip.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_text(data = p_values, aes(x = Inf, y = -Inf, label = paste("p =", format(p_value, digits = 2))),
            hjust = 1.1, vjust = -1.1, size = 5, inherit.aes = FALSE)

# Combine all plots
first_column <- (plot_Coding | plot_Noncoding) / (plot_Total | plot_age)
final_combined_plot <- first_column + clinical_box + plot_layout(ncol = 2, widths = c(3, 1))

# Save plot
ggsave("/mnt/data/Final_Combined_Plot.png", plot = final_combined_plot, width = 16, height = 12, dpi = 300)

