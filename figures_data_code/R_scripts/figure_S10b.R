###Figure S10b###
library(ggplot2)
library(dplyr)
library(broom)
df <- read.csv("../R_data/adaptive_haploid_no_SMF2.csv")
X_VAR <- "latest_stationary_phase_mean_change_from_ancestor"
Y_VAR <- "early_stationary_phase_mean_change_from_ancestor"
FACET_VAR <- "ancestor"
# Ancestor label mapping
label_map <- c(
    "Jason_ancestor"   = "This study",
    "Levy_ancestor"    = "Levy et al. 2015",
    "Yuping_ancestor"  = "Li et al. 2019"
)
df$ancestor_label <- label_map[df[[FACET_VAR]]]
df <- df %>% filter(!is.na(ancestor_label))
df$ancestor_label <- droplevels(factor(df$ancestor_label))
df <- df %>% filter(!is.na(.data[[X_VAR]]) & !is.na(.data[[Y_VAR]]))
# Error columns
x_pos <- which(names(df) == X_VAR)
y_pos <- which(names(df) == Y_VAR)
X_ERR_VAR <- names(df)[x_pos + 1]
Y_ERR_VAR <- names(df)[y_pos + 1]
# Statistics per ancestor group
stats_df <- df %>%
    group_by(ancestor_label) %>%
    summarise(
        n = n(),
        r2 = if(n >= 3) {
            fit <- lm(.data[[Y_VAR]] ~ .data[[X_VAR]])
            summary(fit)$r.squared
        } else NA_real_,
        pval = if(n >= 3) {
            fit <- lm(.data[[Y_VAR]] ~ .data[[X_VAR]])
            coef(summary(fit))[2, 4]
        } else NA_real_,
        .groups = "drop"
    ) %>%
    filter(!is.na(r2), !is.na(pval)) %>%
    arrange(desc(ancestor_label)) %>%
    mutate(
        label = paste0("R² = ", round(r2, 3), "\np = ", signif(pval, 2)),
        x = 0.05,
        y = 0.12 - (0.025 * (row_number() - 1))
    )
# Define consistent color mapping for all ancestors
ancestor_colors <- c(
    "Levy et al. 2015" = "#F8766D",
    "Li et al. 2019" = "#00BA38",
    "This study" = "#619CFF"
)
# ============== Plot ==============
p <- ggplot(df, aes(x = .data[[X_VAR]], y = .data[[Y_VAR]], color = ancestor_label, shape = ploidy_consensus)) +
    
    # Grey dashed lines at x=0 and y=0
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
    
    # Error crosses
    geom_segment(
        aes(x = .data[[X_VAR]] - .data[[X_ERR_VAR]], xend = .data[[X_VAR]] + .data[[X_ERR_VAR]], y = .data[[Y_VAR]], yend = .data[[Y_VAR]]),
        linewidth = 0.6
    ) +
    geom_segment(
        aes(x = .data[[X_VAR]], xend = .data[[X_VAR]], y = .data[[Y_VAR]] - .data[[Y_ERR_VAR]], yend = .data[[Y_VAR]] + .data[[Y_ERR_VAR]]),
        linewidth = 0.6
    ) +
    
    # Data points
    geom_point(size = 3, stroke = 1.1) +
    
    # Trendlines
    geom_smooth(
        method = "lm", 
        se = FALSE,
        aes(group = ancestor_label), 
        color = "black",        
        linewidth = 2.0,         
        show.legend = FALSE
    ) +
    geom_smooth(
        method = "lm", 
        se = FALSE,
        aes(group = ancestor_label), 
        linewidth = 1.5          
    ) +
    
    # Stats text
    geom_text(
        data = stats_df,
        aes(x = x, y = y, label = label, color = ancestor_label),
        inherit.aes = FALSE,
        hjust = 0, vjust = 1, size = 4, show.legend = FALSE
    ) +
    
    # Manual color and shape scales
    scale_color_manual(values = ancestor_colors, 
                       breaks = names(ancestor_colors)) +
    
    # --- MODIFIED: Added 'guide = "none"' to remove the ploidy legend ---
    scale_shape_manual(
        name = "Ploidy",
        values = c("haploid" = 16, "diploid" = 17, "undetermined" = 4), # circle, triangle, X
        labels = c("haploid" = "Haploid", "diploid" = "Diploid", "undetermined" = "Undetermined"),
        na.value = 4,
        guide = "none" # This line removes the shape legend
    ) +
    
    # Axis labels
    labs(
        color = "Ancestor",
        y = "Change in early stationary phase performance following Gly/Eth-limitation",
        x = "Change in latest stationary phase performance following Gly/Eth-limitation"
    ) +
    theme_bw(base_size = 14) +
    theme(
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15)
    )
ggsave("figures/fig_S10b.pdf", p, width = 8, height = 6)
