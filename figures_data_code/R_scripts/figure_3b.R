###figure 3b###
library(ggplot2)
library(dplyr)
library(broom)
library(grid) # For unit()
# ===============================
# Data Wrangling (Unchanged)
# ===============================
df <- read.csv("../R_data/adaptive_haploid_lineages.csv")
X_VAR <- "earliest_stationary_phase_mean_change_from_ancestor"
Y_VAR <- "glucose_stationary_phase_mean_change_from_ancestor"
FACET_VAR <- "ancestor"
# Ancestor label mapping
label_map <- c(
    "Jason_ancestor"   = "This study",
    "Levy_ancestor"    = "Levy et al. 2015",
    "Yuping_ancestor"  = "Li et al. 2019"
)
df$ancestor_label <- label_map[df[[FACET_VAR]]]
# Evolution condition code mapping & order
old_codes <- c("M3_1day", "M3_2day", "M3_5day", "M3_1_5day",
               "M05_2day", "M05_4day", "M05_6day", "M05_8day", "M05_10day")
new_names <- c(
    "Glucose 1 day", "Glucose 2 day", "Glucose 5 day", "Glucose 1/5 day",
    "Gly/Eth 2 day", "Gly/Eth 4 day", "Gly/Eth 6 day", "Gly/Eth 8 day", "Gly/Eth 10 day"
)
desired_order <- new_names
evol_map <- setNames(new_names, old_codes)
df$evol_cond_label <- evol_map[df$evol_cond]
df$evol_cond_label <- factor(df$evol_cond_label, levels = desired_order)
df$evol_cond_label <- droplevels(df$evol_cond_label)
df$ancestor_label <- droplevels(factor(df$ancestor_label,
                                       levels = c("This study", "Levy et al. 2015", "Li et al. 2019")))
df <- df %>% filter(!is.na(ancestor_label) & !is.na(evol_cond_label))
df <- df %>% filter(!is.na(.data[[X_VAR]]) & !is.na(.data[[Y_VAR]]))
# Filter panels with >=2 points
panel_counts <- df %>%
    group_by(evol_cond_label) %>%
    tally(name = "n") %>%
    filter(n >= 2)
df <- df %>% filter(evol_cond_label %in% panel_counts$evol_cond_label)
df$evol_cond_label <- droplevels(df$evol_cond_label)
ancestor_colors <- c(
    "Levy et al. 2015" = "#F8766D",   # RED
    "This study" = "#619CFF",         # BLUE
    "Li et al. 2019" = "#00BA38"      # GREEN
)
x_pos <- which(names(df) == X_VAR)
y_pos <- which(names(df) == Y_VAR)
X_ERR_VAR <- names(df)[x_pos + 1]
Y_ERR_VAR <- names(df)[y_pos + 1]
stats_x <- -0.15
stats_y <- 0.095
stats_df <- df %>%
    group_by(evol_cond_label) %>%
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
    mutate(
        label = paste0("R² = ", round(r2, 3), "\np = ", signif(pval, 2)),
        x = stats_x,
        y = stats_y
    )
stats_df <- stats_df %>% filter(evol_cond_label %in% levels(df$evol_cond_label))
# ============= Final Plot =============
p <- ggplot(df, aes(
    x = .data[[X_VAR]],
    y = .data[[Y_VAR]],
    color = ancestor_label
)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.8) +
    geom_segment(aes(x=.data[[X_VAR]] - .data[[X_ERR_VAR]], xend=.data[[X_VAR]] + .data[[X_ERR_VAR]], y=.data[[Y_VAR]], yend=.data[[Y_VAR]]), linewidth=0.6) +
    geom_segment(aes(x=.data[[X_VAR]], xend=.data[[X_VAR]], y=.data[[Y_VAR]] - .data[[Y_ERR_VAR]], yend=.data[[Y_VAR]] + .data[[Y_ERR_VAR]]), linewidth=0.6) +
    geom_point(size = 3) +
    geom_smooth(data= . %>% filter(evol_cond_label != "Gly/Eth 8 day"), method="lm", se=FALSE, aes(group=ancestor_label), color="black", linewidth=2.0, show.legend=FALSE) +
    geom_smooth(data= . %>% filter(evol_cond_label != "Gly/Eth 8 day"), method="lm", se=FALSE, aes(group=ancestor_label, color=ancestor_label), linewidth=1.5) +
    geom_text(data=stats_df, aes(x=x, y=y, label=label), inherit.aes=FALSE, hjust=0, vjust=1, size=5, color="black", show.legend=FALSE) +
    scale_color_manual(name = "Ancestor", values = ancestor_colors) +
    labs(
        color = "Ancestor",
        x = "Change in earliest stationary phase performance following Gly/Eth-limitation",
        y = "Change in stationary phase performance following Glucose-limitation"
    ) +
    theme_bw(base_size = 14) +
    facet_wrap(~ evol_cond_label, nrow = 3, ncol = 3) +

    # --- MODIFIED: Manual legend positioning and styling ---
    theme(
        legend.position = c(0.85, 0.2),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20)
    )
ggsave("figures/fig_3b.pdf", p, width = 12, height = 10)
