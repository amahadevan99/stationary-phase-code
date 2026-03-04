###Figure S5
###plotting SMF2 and chr11 duplications in the 6-day condition
#-----------------------------------------------------------------------
# FINAL SCRIPT: PLOTTING "Gly/Eth 6 day" WITH LEGEND ON THE RIGHT
#-----------------------------------------------------------------------
# Load the necessary libraries
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
# --- Step 1: Create a Grouping Column Based on Your Criteria (Unchanged) ---
original_df <- read.csv("../R_data/adaptive_lineages.csv")
df_with_groups <- original_df %>%
    dplyr::mutate(dplyr::across(starts_with("gene.annotation_"), as.character)) %>%
    dplyr::mutate(
      has_chr11_dup = (chr11_duplication == "yes"),
      has_smf2 = dplyr::if_any(starts_with("gene.annotation_"), ~str_detect(., regex("SMF2", ignore_case = TRUE)), .na = FALSE)
    ) %>%
    dplyr::mutate(
      plot_group = case_when(
        has_chr11_dup & has_smf2   ~ "chr11_and_SMF2",
        has_chr11_dup              ~ "chr11_dup",
        has_smf2                   ~ "SMF2_only",
        TRUE                       ~ "Other"
      )
    )
# --- Step 2: Filter for ONLY the Clones of Interest (Unchanged) ---
filtered_clones <- df_with_groups %>%
  filter(plot_group %in% c("chr11_and_SMF2", "chr11_dup", "SMF2_only"))
# --- Step 3: Reshape the Data and Apply All Filters (Unchanged) ---
plot_data_unfactored <- filtered_clones %>%
    dplyr::select(plot_group, matches("_mean_.*(Fitness|Error)_.*_from_ancestor")) %>%
    pivot_longer(
        cols = -plot_group,
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
    ) %>%
    filter(!is.na(fitness)) %>%
    filter(treatment_label == "Gly/Eth 6 day")
# --- Step 4: Prepare Data for Plotting (Unchanged) ---
group_order <- c("chr11_and_SMF2", "chr11_dup", "SMF2_only")
plot_data <- plot_data_unfactored %>%
  mutate(
    plot_group = factor(plot_group, levels = group_order)
  )
# --- Step 5: Create the Jitter Plot ---
p <- ggplot(
  data = plot_data, 
  aes(x = treatment_label, y = fitness, color = plot_group)
) +
  
  geom_jitter(
    width = 0.25,
    height = 0,
    alpha = 0.6,
    size = 3
  ) +
  
  scale_color_manual(
    name = "Genotype",
    values = c(
      "chr11_and_SMF2" = "red", 
      "chr11_dup"      = "blue", 
      "SMF2_only"      = "darkgreen"
    ),
    labels = c(
      expression(paste("chr11 duplication & ", italic("SMF2"))),
      "chr11 duplication",
      expression(italic("SMF2"))
    )
  ) +
  labs(
    x = "Assay Environment",
    y = "Mean Fitness Change Per Cycle"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    
    # *** CHANGED: Legend position is now "right" ***
    legend.position = "right"
  )
ggsave("figures/fig_S5.pdf", p, width = 6, height = 6)
