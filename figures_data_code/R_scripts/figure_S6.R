##Figure S6
###FZF1 diploids with ploidy colored
#-----------------------------------------------------------------------
# FINAL SCRIPT WITH UPDATED "FZF1 Diploid" LABEL
#-----------------------------------------------------------------------
# Load the necessary libraries
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
# --- Step 1: Add a New, Combined Status Column ---
original_df <- read.csv("../R_data/adaptive_lineages.csv")
df_with_plot_group <- original_df %>%
    dplyr::mutate(dplyr::across(starts_with("gene.annotation_"), as.character)) %>%
    dplyr::mutate(
      has_fzf = dplyr::if_else(
        condition = dplyr::if_any(starts_with("gene.annotation_"), ~str_detect(., regex("FZF", ignore_case = TRUE))),
        true      = "Yes",
        false     = "No",
        missing   = "No"
      )
    ) %>%
    dplyr::mutate(
      plot_group = case_when(
        # *** CHANGED: The value is now "FZF1 Diploid" ***
        has_fzf == "Yes" ~ "FZF1 Diploid",
        
        ploidy_consensus == "diploid" ~ "Diploid",
        TRUE ~ "Haploid"
      )
    )
# --- Step 2: Reshape Data and Apply Filter ---
plot_data_unfactored <- df_with_plot_group %>%
    dplyr::select(plot_group, matches("_mean_.*(Fitness|Error)_.*_from_ancestor")) %>%
    
    pivot_longer(
        cols = -plot_group,
        names_to = c("treatment_group", ".value"),
        names_pattern = "(.*_mean)_(Fitness|Error).*"
    ) %>%
    rename(fitness = Fitness, error = Error) %>%
    # *** CHANGED: The filter now looks for the new label ***
    filter(plot_group == "FZF1 Diploid" | error < 5) %>%
    
    mutate(
        treatment_label = treatment_group %>%
            str_remove("_mean$") %>%
            str_replace_all(c("M05" = "Gly/Eth", "M3" = "Glucose")) %>%
            str_replace_all("_", " ")
    ) %>%
    filter(!is.na(fitness))
# --- Step 3: Prepare Data for Plotting (Set Factor Levels) ---
unique_labels <- unique(plot_data_unfactored$treatment_label)
gly_eth_labels <- unique_labels[str_detect(unique_labels, "Gly/Eth")]
glucose_labels <- unique_labels[str_detect(unique_labels, "Glucose")]
sorted_gly_eth <- str_sort(gly_eth_labels, numeric = TRUE)
sorted_glucose <- str_sort(glucose_labels, numeric = TRUE)
final_sorted_labels <- c(sorted_gly_eth, sorted_glucose)
# *** CHANGED: The group order now uses the new label ***
group_order <- c("FZF1 Diploid", "Diploid", "Haploid")
plot_data <- plot_data_unfactored %>%
  mutate(
    treatment_label = factor(treatment_label, levels = final_sorted_labels),
    plot_group = factor(plot_group, levels = group_order)
  ) %>%
  arrange(desc(plot_group))
# --- Step 4: Create the Final Plot with Updated Legend ---
p <- ggplot(
  data = plot_data, 
  aes(x = treatment_label, y = fitness, color = plot_group, size = plot_group)
) +
  
  geom_jitter(
    width = 0.25, 
    height = 0,
    alpha = 0.7
  ) +
  
  # *** CHANGED: Updated legend title and label name ***
  scale_color_manual(
    name = "Genotype",
    values = c("FZF1 Diploid" = "red", "Diploid" = "blue", "Haploid" = "black"),
    labels = c("FZF1 Diploid" = expression(paste(italic("FZF1"), " Diploid")), "Diploid" = "Diploid", "Haploid" = "Haploid")
  ) +

  # *** CHANGED: Updated legend title and label name ***
  scale_size_manual(
    name = "Genotype",
    values = c("FZF1 Diploid" = 3.5, "Diploid" = 2.5, "Haploid" = 1.5),
    labels = c("FZF1 Diploid" = expression(paste(italic("FZF1"), " Diploid")), "Diploid" = "Diploid", "Haploid" = "Haploid")
  ) +
  labs(
    x = "Assay Environment",
    y = "Mean Fitness Change Per Cycle"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    legend.position = "top"
  )
ggsave("figures/fig_S6.pdf", p, width = 10, height = 6)
