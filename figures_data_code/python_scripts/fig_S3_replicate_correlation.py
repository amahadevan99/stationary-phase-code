import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from matplotlib.backends.backend_pdf import PdfPages
import re

# ---------- Configuration ----------
input_file = "../python_data/fitseq2_fitness_all_replicates.csv"
output_pdf_glucose = os.path.join("figures", "fig_S3_glucose.pdf")
output_pdf_glyeth = os.path.join("figures", "fig_S3_glyeth.pdf")

df = pd.read_csv(input_file)
fitness_cols = [col for col in df.columns if col.endswith("_Fitness_Per_Cycle")]
error_cols = [col for col in df.columns if col.endswith("_Error_Fitness")]

def exp_name(colname):
    parts = colname.split('_')
    return '_'.join(parts[:3])

def replicate_number(colname):
    parts = colname.split('_')
    for p in parts:
        if p.startswith('r'):
            return p
    return None

def days_from_expname(exp):
    try:
        return f"{int(exp.split('_')[1])} day"
    except:
        return ""

def day_number_from_expname(exp):
    try:
        return int(re.search(r'_([0-9]+)_day', exp).group(1))
    except:
        return 0

# Build experiment lookup tables
experiments = {}
errors_by_exp = {}
for col in fitness_cols:
    exp = exp_name(col)
    rep = replicate_number(col)
    experiments.setdefault(exp, {})[rep] = col
for col in error_cols:
    exp = exp_name(col)
    rep = replicate_number(col)
    errors_by_exp.setdefault(exp, {})[rep] = col

exp_list_all = sorted(experiments.keys())
exp_list_m3 = sorted([exp for exp in exp_list_all if exp.startswith("M3")], key=day_number_from_expname)
exp_list_m05 = sorted([exp for exp in exp_list_all if exp.startswith("M05")], key=day_number_from_expname)

replicate_pairs = [('r1', 'r2'), ('r1', 'r3'), ('r2', 'r3')]
rep_label = {'r1': "R1 fitness", 'r2': "R2 fitness", 'r3': "R3 fitness"}
pair_names = {('r1', 'r2'): "R1 vs R2", ('r1', 'r3'): "R1 vs R3", ('r2', 'r3'): "R2 vs R3"}

def make_panel(exp_list, output_pdf, panel_title):
    with PdfPages(output_pdf) as pdf:
        fig, axes = plt.subplots(nrows=len(exp_list), ncols=3, figsize=(15, max(3*len(exp_list), 6)), sharex=False, sharey=False)
        if len(exp_list) == 1:
            axes = np.expand_dims(axes, 0)  # ensure axes is always 2d

        plt.subplots_adjust(left=0.09, right=0.995, top=0.93, bottom=0.07, hspace=0.5, wspace=0)

        for exp_idx, exp in enumerate(exp_list):
            reps = experiments[exp]
            row_values = []
            for rep in ['r1', 'r2', 'r3']:
                if rep in reps:
                    col = reps[rep]
                    err_col = errors_by_exp[exp][rep]
                    vals = df.loc[(df[err_col] < 5) & df[col].notna(), col]
                    row_values.append(vals)
            if row_values and any(len(v) > 0 for v in row_values):
                all_vals = pd.concat(row_values)
                val_min = all_vals.min()
                val_max = all_vals.max()
                if val_min == val_max:
                    minval = val_min - 1
                    maxval = val_max + 1
                else:
                    minval = val_min
                    maxval = val_max
            else:
                minval = 0
                maxval = 1

            lims = (minval, maxval)
            locator = plt.MaxNLocator(nbins=6)
            ticks = locator.tick_values(minval, maxval)

            for pair_idx, (rep1, rep2) in enumerate(replicate_pairs):
                ax = axes[exp_idx, pair_idx]
                if rep1 in reps and rep2 in reps:
                    col1 = reps[rep1]
                    col2 = reps[rep2]
                    err1 = errors_by_exp[exp][rep1]
                    err2 = errors_by_exp[exp][rep2]
                    mask = (df[err1] < 5) & (df[err2] < 5) & df[col1].notna() & df[col2].notna()
                    x = df.loc[mask, col1]
                    y = df.loc[mask, col2]
                    xerr = df.loc[mask, err1]
                    yerr = df.loc[mask, err2]
                    if len(x) > 1 and len(y) > 1:
                        ax.errorbar(
                            x, y,
                            xerr=xerr, yerr=yerr,
                            fmt='none', ecolor='gray', alpha=0.4, zorder=1, capsize=2,
                            elinewidth=0.5, capthick=0.5,
                            rasterized=True
                        )
                        ax.scatter(x, y, s=10, zorder=2, rasterized=True)
                        if np.abs(x.max() - x.min()) > 0 and np.abs(y.max() - y.min()) > 0:
                            coef = np.polyfit(x, y, 1)
                            fit_x = np.array([minval, maxval])
                            fit_y = np.polyval(coef, fit_x)
                            ax.plot(fit_x, fit_y, color="blue", lw=2, alpha=0.7, label="Line of fit")
                        r2 = r2_score(x, y)
                        # Title logic:
                        if exp_idx == 0:
                            title_str = f"{pair_names[(rep1, rep2)]}\nR²={r2:.3f}, n={len(x)}"
                        else:
                            title_str = f"R²={r2:.3f}, n={len(x)}"
                        ax.set_title(title_str, fontsize=12)
                        ax.set_xlabel(rep_label[rep1], fontsize=11)
                        ax.set_ylabel(rep_label[rep2], fontsize=11)
                        ax.set_xlim(lims)
                        ax.set_ylim(lims)
                        ax.set_aspect('equal', adjustable='box')
                        ax.plot(lims, lims, 'r--', alpha=0.5, label="x=y")
                        ax.legend(loc="best", fontsize=9)
                    else:
                        if exp_idx == 0:
                            ax.set_title(f"{pair_names[(rep1, rep2)]}\nNot enough data", fontsize=12)
                        else:
                            ax.set_title(f"Not enough data", fontsize=12)
                        ax.set_xticks([])
                        ax.set_yticks([])
                else:
                    if exp_idx == 0:
                        ax.set_title(f"{pair_names[(rep1, rep2)]}\nReplicate(s) missing", fontsize=12)
                    else:
                        ax.set_title("Replicate(s) missing", fontsize=12)
                    ax.set_xticks([])
                    ax.set_yticks([])
                ax.set_xticks(ticks)
                ax.set_yticks(ticks)
                ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f"{val:.2g}"))
                ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f"{val:.2g}"))
            ax_left = axes[exp_idx, 0]
            ax_left.annotate(
                days_from_expname(exp),
                xy=(-0.32, 0.5), xycoords="axes fraction",
                fontsize=12, ha='right', va='center', fontweight='bold', rotation=90,
                annotation_clip=False
            )

        fig.suptitle(panel_title, fontsize=19, y=0.97, fontweight='bold')
        plt.tight_layout(rect=[0.10, 0, 1, 0.94])
        pdf.savefig(fig, dpi=300)
        plt.close(fig)
    print(f"Done! {panel_title} PDF saved to {output_pdf}")

make_panel(exp_list_m3, output_pdf_glucose, 'Glucose')
make_panel(exp_list_m05, output_pdf_glyeth, 'Gly/Eth')
