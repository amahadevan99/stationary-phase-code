import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import os

file_path = '../python_data/adaptive_fitness_home_condition.csv'
df = pd.read_csv(file_path)

df["days between transfers"] = pd.to_numeric(df["days between transfers"], errors='coerce')
df = df.dropna(subset=["days between transfers", "fitness_change_in_evol_cond", "evolution media"])

media_list = ["glucose", "Gly/Eth"]
xtick_dict = {
    "glucose": [1, 2, 3, 4, 5],
    "Gly/Eth": [2, 4, 6, 8, 10]
}

fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

for ax, media in zip(axes, media_list):
    sub = df[df["evolution media"] == media]
    x = sub["days between transfers"]
    y = sub["fitness_change_in_evol_cond"]

    jitter = 0.3  # Increased from 0.15 to 0.3 for more spread
    x_jittered = x + (jitter * (np.random.rand(len(x)) - 0.5))
    ax.scatter(x_jittered, y, alpha=0.7, color="black")

    # Add horizontal lines for means at each x value
    unique_x = sorted(x.unique())
    for x_val in unique_x:
        y_vals = y[x == x_val]
        if len(y_vals) > 0:
            mean_y = y_vals.mean()
            # Add horizontal line at the mean, spanning a small width
            line_width = 0.4
            ax.hlines(y=mean_y, xmin=x_val - line_width, xmax=x_val + line_width,
                     color='red', linewidth=2.5)

    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    x_fit = np.array(xtick_dict[media])
    y_fit = slope * x_fit + intercept
    ax.plot(x_fit, y_fit, color='blue')

    stat_txt = f"R² = {r_value**2:.3f}\nP = {p_value:.2g}"
    y_max = y.max()
    ax.text(x_fit[0], y_max * 0.97, stat_txt, va='top', ha='left', fontsize=16, color='black')

    ax.set_xticks(xtick_dict[media])
    ax.tick_params(axis='both', labelsize=14)

    # Capitalize "glucose"
    panel_title = media.capitalize() if media == "glucose" else media
    ax.set_title(panel_title, fontsize=18)

plt.tight_layout()

# Put axis labels even further away from the figure
fig.text(0.5, -0.03, "days between transfers", ha="center", fontsize=18)    # x-label (lower than before)
fig.text(-0.01, 0.5, "fitness effect per cycle in evolution condition", va="center", rotation="vertical", fontsize=18)  # y-label (further left than before)

outfile = os.path.join("figures", "fig_2b.pdf")

plt.savefig(outfile, bbox_inches="tight", format="pdf")
# plt.show()

print(f"PDF saved to: {os.path.abspath(outfile)}")
