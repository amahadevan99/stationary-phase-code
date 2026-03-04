###Figure S8###
library(ggplot2)
library(dplyr)
library(broom)
# ===============================
# Data Wrangling (Unchanged)
# ===============================
df <- read.csv("../R_data/adaptive_haploid_lineages.csv")
# Filter for ONLY the condition of interest
plot_df <- df %>%
    filter(evol_cond == "M05_6day") 
X_VAR <- "late_stationary_performance_mean_change_from_ancestor"
Y_VAR <- "earliest_stationary_phase_mean_change_from_ancestor"
# Ensure the columns and the grouping variable exist and are not empty
plot_df <- plot_df %>%
    filter(!is.na(.data[[X_VAR]]) & !is.na(.data[[Y_VAR]]) & !is.na(Contains_SMF2))
# Error columns
x_pos <- which(names(plot_df) == X_VAR)
y_pos <- which(names(plot_df) == Y_VAR)
X_ERR_VAR <- names(plot_df)[x_pos + 1]
Y_ERR_VAR <- names(plot_df)[y_pos + 1]
# Calculate statistics for each SMF2 group
stats_df <- plot_df %>%
    group_by(Contains_SMF2) %>% 
    summarise(
        n = n(),
        r2 = if(n >= 3) {
            fit <- lm(.data[[Y_VAR]] ~ .data[[X_VAR]])
            summary(fit)$r.squared
        } else NA_real_,
        pval = if(n >= 3) {
            fit <- lm(.data[[Y_VAR]] ~ .data[[X_VAR]])
            coef(summary(fit))[2,4]
        } else NA_real_,
        .groups = "drop"
    ) %>%
    filter(!is.na(r2), !is.na(pval)) %>%
    arrange(Contains_SMF2) %>% 
    mutate(
        label = paste0("R² = ", round(r2, 3), "\np = ", signif(pval, 2)),
        x = 0.0,
        y = 0.125 - (0.015 * (row_number() - 1)) 
    )
# --- MODIFIED: Using the "salmon" and "teal" default ggplot colors ---
smf2_colors <- c(
    `TRUE` = "#00BFC4",  # Default Teal/Cyan for SMF2 mutations
    `FALSE`  = "#F8766D"   # Default Salmon/Red for other points
)
# ============= Final Plot =============
p <- ggplot(plot_df, aes(
    x = .data[[X_VAR]],
    y = .data[[Y_VAR]],
    color = Contains_SMF2
)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.8) +
    geom_segment(aes(x=.data[[X_VAR]] - .data[[X_ERR_VAR]], xend=.data[[X_VAR]] + .data[[X_ERR_VAR]], y=.data[[Y_VAR]], yend=.data[[Y_VAR]]), linewidth=0.6) +
    geom_segment(aes(x=.data[[X_VAR]], xend=.data[[X_VAR]], y=.data[[Y_VAR]] - .data[[Y_ERR_VAR]], yend=.data[[Y_VAR]] + .data[[Y_ERR_VAR]]), linewidth=0.6) +
    geom_point(size = 3) +
    
    # Trendlines
    geom_smooth(method = "lm", se = FALSE, aes(group = Contains_SMF2), color = "black", linewidth = 2.0, show.legend = FALSE) +
    geom_smooth(method = "lm", se = FALSE, aes(group = Contains_SMF2, color = Contains_SMF2), linewidth = 1.5) +
    
    # Stats text
    geom_text(data = stats_df, aes(x = x, y = y, label = label, color = Contains_SMF2), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4, show.legend = FALSE) +
    
    # Color scale (now uses the new color vector)
    scale_color_manual(
        name = expression(paste(italic("SMF2"), " Mutation")), 
        values = smf2_colors,
        labels = c(`TRUE` = "Yes", `FALSE` = "No")
    ) +
    labs(
        title = "Gly/Eth 6 day",
        x = "Change in late stationary phase performance following Gly/Eth-limitation",
        y = "Change in earliest stationary phase performance following Gly/Eth-limitation"
    ) +
    theme_bw(base_size = 14) +
    theme(
        legend.position = "right",
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15)
    )
ggsave("figures/fig_S8.pdf", p, width = 8, height = 6)
