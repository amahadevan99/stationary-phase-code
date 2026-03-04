import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import math
import os  # Add this import

# --- Load data ---
data_path = '../python_data/chimeric_predictions_empty_combos.csv'
df = pd.read_csv(data_path)

# --- Find all *_count columns and matching *_prediction_my_method columns ---
count_cols = [col for col in df.columns if col.endswith('_count')]
pairs = []
for count_col in count_cols:
    pred_col = count_col.replace('_count', '_prediction_my_method')
    if pred_col in df.columns:
        pairs.append((count_col, pred_col))

if len(pairs) == 0:
    print("\nNo suitable count/prediction pairs found. Check your column names!")
    exit(1)
else:
    print(f"\nFound {len(pairs)} pairs for plotting:")
    for count_col, pred_col in pairs:
        print(f"  {count_col} — {pred_col}")

n_plots = len(pairs)
n_cols = 2
n_rows = math.ceil(n_plots / n_cols)

pdf_filename = os.path.join("figures", "fig_S15.pdf")

# --- Plot ---
fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows), squeeze=False)

for idx, (count_col, pred_col) in enumerate(pairs):
    ax = axes[idx // n_cols][idx % n_cols]
    x = df[pred_col]
    y = df[count_col]

    # Remove zeros and negatives for log plotting
    mask = (x > 0) & (y > 0)
    x = x[mask]
    y = y[mask]

    # Scatter
    sns.scatterplot(x=x, y=y, ax=ax, s=40)

    # Linear regression on log-transformed data
    log_x = np.log10(x)
    log_y = np.log10(y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_x, log_y)
    x_sorted = pd.Series(x).sort_values()
    log_x_sorted = np.log10(x_sorted)
    ax.plot(x_sorted, 10**(slope * log_x_sorted + intercept), color='blue', lw=2)

    # R^2 annotation
    ax.text(0.05, 0.92, f"$R^2 = {r_value**2:.3f}$", transform=ax.transAxes, fontsize=14, color='black')

    # Custom labels: FOS#-ROS# predicted/count
    # Extract FOS#-ROS# from the count_col name
    group_label = count_col.replace('_count','')
    group_label = group_label.replace('_','-')  # FOS2_ROS11 -> FOS2-ROS11

    ax.set_xlabel(f"{group_label} predicted", fontsize=11)
    ax.set_ylabel(f"{group_label} count", fontsize=11)
    ax.set_title(f"{group_label}: predicted vs count", fontsize=12)

    # Set both axes to log scale
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Add an x=y dashed line (1:1 line) over the log-log range
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    minval = max(xlims[0], ylims[0])
    maxval = min(xlims[1], ylims[1])
    x_to_plot = np.logspace(np.log10(minval), np.log10(maxval), 100)
    ax.plot(x_to_plot, x_to_plot, color='gray', linestyle='--', linewidth=1, label='x=y')

    # Optionally show legend just once, or none:
    ax.legend(['Regression', '$x=y$'], loc='lower right', fontsize=10)

# Remove empty subplots
for idx in range(n_plots, n_rows*n_cols):
    fig.delaxes(axes[idx // n_cols][idx % n_cols])

plt.tight_layout()
plt.savefig(pdf_filename)
plt.close(fig)
print(f"Saved PDF output to: {os.path.abspath(pdf_filename)}")
