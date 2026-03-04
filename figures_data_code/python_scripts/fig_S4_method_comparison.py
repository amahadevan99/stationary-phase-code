import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from matplotlib.backends.backend_pdf import PdfPages

input_file = "../python_data/fitness_data_comparison.csv"

output_pdf = os.path.join("figures", "fig_S4.pdf")

df = pd.read_csv(input_file)

def match_col(substring, columns):
    for col in columns:
        if substring in col:
            return col
    raise KeyError(f"No column containing '{substring}' found. Available columns: {columns}")

columns = df.columns.tolist()
conditions = [
    {
        "label": "1 day",
        "cols": [
            match_col("1D_Fitness_Li2019_Li2019Neutrals", columns),
            match_col("1day_Grant_corrected", columns),
            match_col("M3_1_day_mean_Fitness_Per_Cycle_change_from_ancestor", columns)
        ]
    },
    {
        "label": "2 day",
        "cols": [
            match_col("2D_Fitness_Li2019_Li2019Neutrals", columns),
            match_col("2day_Grant_corrected", columns),
            match_col("M3_2_day_mean_Fitness_Per_Cycle_change_from_ancestor", columns)
        ]
    },
    {
        "label": "3 day",
        "cols": [
            match_col("3D_Fitness_Li2019_Li2019Neutrals", columns),
            match_col("3day_Grant_corrected", columns),
            match_col("M3_3_day_mean_Fitness_Per_Cycle_change_from_ancestor", columns)
        ]
    },
    {
        "label": "5 day",
        "cols": [
            match_col("5D_Fitness_Li2019_Li2019Neutrals", columns),
            match_col("5day_Grant_corrected", columns),
            match_col("M3_5_day_mean_Fitness_Per_Cycle_change_from_ancestor", columns)
        ]
    },
]

method_labels = ['Li et al. 2019', 'Kinsler et al. 2024', 'This study']

ancestor_groups = df['Ancestor'].unique()
def display_ancestor_label(x):
    if x == "Levy_ancestor":
        return "Levy et al. 2015"
    elif x == "Yuping_ancestor":
        return "Li et al. 2019"
    else:
        return x

colors = plt.cm.get_cmap('tab10', len(ancestor_groups))

# MathText bold with spaces
def make_bold(text):
    return r"$\bf{" + "}~{".join(text.split(" ")) + "}$"

with PdfPages(output_pdf) as pdf:
    fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(18, 15))
    plt.subplots_adjust(hspace=0.45, wspace=0.12)
    legend_handles = []
    legend_labels = []
    for cond_idx, cond in enumerate(conditions):
        cols = cond["cols"]
        pairs = [(0,1), (0,2), (1,2)]
        for pidx, (i,j) in enumerate(pairs):
            ax = axes[cond_idx, pidx]
            stats_lines = []
            handles = []
            for group_idx, ancestor in enumerate(ancestor_groups):
                group_df = df[df['Ancestor'] == ancestor]
                x = group_df[cols[i]]
                y = group_df[cols[j]]
                mask = x.notna() & y.notna()
                color = colors(group_idx)
                label = display_ancestor_label(ancestor)
                sc = ax.scatter(x[mask], y[mask], s=12, color=color, label=label, alpha=0.8)
                if cond_idx == 3 and pidx == 2:
                    legend_handles.append(sc)
                    legend_labels.append(label)
                r2 = np.nan
                n = mask.sum()
                if n > 1:
                    try:
                        coef = np.polyfit(x[mask], y[mask], 1)
                        y_pred = np.polyval(coef, x[mask])
                        fit_x = np.array([x[mask].min(), x[mask].max()])
                        fit_y = np.polyval(coef, fit_x)
                        ax.plot(fit_x, fit_y, color=color, lw=2, alpha=0.7)
                        r2 = r2_score(y[mask], y_pred)
                        r2 = max(r2, 0)
                    except Exception:
                        pass
                stats_lines.append(f"{label}: R²={r2:.3f}, n={n}")

            all_x = df[cols[i]]
            all_y = df[cols[j]]
            mask_all = all_x.notna() & all_y.notna()
            if mask_all.sum() > 0:
                minlim = min(all_x[mask_all].min(), all_y[mask_all].min())
                maxlim = max(all_x[mask_all].max(), all_y[mask_all].max())
            else:
                minlim, maxlim = 0, 1
            ax.plot([minlim, maxlim], [minlim, maxlim], 'r--', alpha=0.5)
            ax.set_xlim(minlim, maxlim)
            ax.set_ylim(minlim, maxlim)
            ax.set_xlabel(method_labels[i], fontsize=12)
            ax.set_ylabel(method_labels[j], fontsize=12)
            stats_text = "\n".join(stats_lines)

            # Top row: show bold and spaced label
            if cond_idx == 0:
                title_str = f'Measured in {make_bold(method_labels[i])} vs Measured in {make_bold(method_labels[j])}\n{stats_text}'
            else:
                title_str = stats_text
            ax.set_title(title_str, fontsize=12)

            # Only bottom right shows the legend
            if cond_idx == 3 and pidx == 2:
                ax.legend(legend_handles, legend_labels, title="Evolved in:", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=10)
            else:
                ax.legend().remove()
        axes[cond_idx, 0].annotate(
            cond['label'], xy=(-0.20, 0.5), xycoords="axes fraction",
            fontsize=13, ha='right', va='center', fontweight='bold', rotation=90,
            annotation_clip=False)
    fig.suptitle("Pairwise Comparison Across Methods by Condition and Ancestor", fontsize=17, fontweight='bold', y=1.03)
    plt.tight_layout(rect=[0.07, 0, 1, 0.99])
    pdf.savefig(fig)
    plt.close(fig)
print(f"PDF saved to {output_pdf}")
