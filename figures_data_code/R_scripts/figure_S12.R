###Figure S12
#plotting chr11 duplications
#-----------------------------------------------------------------------
# SCRIPT TO PLACE SIGNIFICANCE LETTERS AT THE TOP OF THE PLOT
#-----------------------------------------------------------------------
# Load the necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
# library(multcomp)
# --- Step 1: Define Data and Apply Initial Filters (No change) ---
base_data <- read.csv("../R_data/adaptive_haploid_lineages.csv") %>%
    filter(chr11_duplication == "yes") %>%
    filter(!if_any(starts_with("gene annotation_"), ~coalesce(str_detect(., "FZF1"), FALSE)))
# --- Step 2: Reshape Data and Apply Final Filter (No change) ---
plot_data_unfactored <- base_data %>%
    dplyr::select(matches("^(M3|M05)_.*(ancestor|error)$")) %>%
    pivot_longer(
        cols = everything(),
        names_to = c("treatment_group", ".value"),
        names_pattern = "(.*_mean)_(Fitness|Error).*"
    ) %>%
    rename(fitness = Fitness, error = Error) %>%
    filter(error < 5) %>%
    mutate(
        treatment_label = treatment_group %>%
            str_remove("_mean$") %>%
            str_replace_all(c("M05" = "Gly/Eth", "M3" = "Glucose")) %>%
            str_replace_all("_", " ")
    )
# --- Step 3: Control the X-Axis Order (No change) ---
unique_labels <- unique(plot_data_unfactored$treatment_label)
gly_eth_labels <- unique_labels[str_detect(unique_labels, "Gly/Eth")]
glucose_labels <- unique_labels[str_detect(unique_labels, "Glucose")]
sorted_gly_eth <- str_sort(gly_eth_labels, numeric = TRUE)
sorted_glucose <- str_sort(glucose_labels, numeric = TRUE)
final_sorted_labels <- c(sorted_gly_eth, sorted_glucose)
separator_position <- length(sorted_gly_eth) + 0.5
plot_data <- plot_data_unfactored %>%
    mutate(treatment_label = factor(treatment_label, levels = final_sorted_labels))
# --- Step 4: Tukey HSD test (commented out) ---
# anova_result <- aov(fitness ~ treatment_label, data = plot_data)
# glht_result <- glht(anova_result, linfct = mcp(treatment_label = "Tukey"))
# cld_result <- cld(glht_result)
# stat_labels <- data.frame(letters = cld_result$mcletters$Letters) %>%
#     mutate(treatment_label = rownames(.))
# plot_max_y <- max(plot_data$fitness, na.rm = TRUE)
# stat_labels$y_pos <- plot_max_y + (abs(plot_max_y) * 0.05)
# --- Step 5: Create the Final Plot ---
p <- ggplot(plot_data, aes(x = treatment_label, y = fitness)) +
    geom_violin(aes(fill = treatment_label), show.legend = FALSE, trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, height = 0) +
    geom_vline(xintercept = separator_position, linetype = "dashed", color = "grey20", linewidth = 0.8) +

    labs(
        x = "Assay Environment",
        y = "Mean Fitness Change Per Cycle"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )
ggsave("figures/fig_S12.pdf", p, width = 10, height = 6)
