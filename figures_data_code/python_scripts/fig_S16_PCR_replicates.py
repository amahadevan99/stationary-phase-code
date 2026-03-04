import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from collections import defaultdict
from matplotlib.backends.backend_pdf import PdfPages
import re

# === CONFIGURATION ===
input_sets = {
    'Raw counts': {
        '1': '../python_data/lane1_raw_counts.csv',
        '2': '../python_data/lane2_raw_counts.csv'
    },
    'Chimera corrected': {
        '1': '../python_data/lane1_chimera_corrected.csv',
        '2': '../python_data/lane2_chimera_corrected.csv'
    }
}
all_samples_file = '../python_data/all_samples.csv'
excluded_samples = {
    "M05_6_day_r3_T2",
    "M05_6_day_r3_T3",
    "M3_3_day_r3_T2",
    "M3_3_day_r3_T3"
}

data_dfs_sets = {
    setname: {lanename: pd.read_csv(path, dtype=float)
              for lanename, path in files.items()}
    for setname, files in input_sets.items()
}
all_samples_df = pd.read_csv(all_samples_file, dtype=str)

def get_column(df, ros, fos):
    ros_num = str(ros).replace('ROS', '').strip()
    fos_num = str(fos).replace('FOS', '').strip()
    base_col = f'FOS{fos_num}_ROS{ros_num}'
    if base_col in df.columns:
        return base_col
    elif (base_col + " Corrected") in df.columns:
        return base_col + " Corrected"
    else:
        return None

def beautify_name(name):
    name = name.replace("M05", "Gly/Eth").replace("M3", "Glucose")
    name = name.replace("_", " ")
    return name

def extract_group_and_day(pretty_sample):
    pretty_sample = pretty_sample.lower()
    group_val = 0 if 'gly' in pretty_sample else 1  # Gly/Eth first
    # Extract day number (assume space before and after) or fallback
    day_found = re.search(r'(\d+)\s*day', pretty_sample)
    if day_found:
        day_val = int(day_found.group(1))
    else:
        day_val = 999
    return group_val, day_val, pretty_sample

sample_groups = defaultdict(list)
for idx, row in all_samples_df.iterrows():
    sample_groups[row['sample']].append(row)

plot_pairs_with_names = []
for sample_name, rows in sample_groups.items():
    if sample_name in excluded_samples:
        continue
    if len(rows) < 2:
        continue
    pretty_sample = beautify_name(sample_name)
    for i in range(len(rows)):
        for j in range(i+1, len(rows)):
            row1, row2 = rows[i], rows[j]
            lane1, ros1, fos1 = row1['Lane'], row1['ROS'], row1['FOS']
            lane2, ros2, fos2 = row2['Lane'], row2['ROS'], row2['FOS']
            plots = []
            for setname, _ in input_sets.items():
                df1 = data_dfs_sets[setname].get(str(lane1))
                df2 = data_dfs_sets[setname].get(str(lane2))
                if df1 is None or df2 is None:
                    continue
                col1 = get_column(df1, ros1, fos1)
                col2 = get_column(df2, ros2, fos2)
                if not col1 or not col2:
                    continue
                x = df1[col1]
                y = df2[col2]
                df_xy = pd.DataFrame({'x': x, 'y': y}).dropna()
                df_xy = df_xy[(df_xy['x'] > 0) & (df_xy['y'] > 0)]
                if len(df_xy) == 0:
                    continue
                df_xy['log_x'] = np.log10(df_xy['x'])
                df_xy['log_y'] = np.log10(df_xy['y'])
                r = df_xy['log_x'].corr(df_xy['log_y'])
                r2 = r ** 2 if pd.notnull(r) else float('nan')
                pretty_col1 = col1.replace("_", " ").replace("Corrected", "").replace("(count)", "").replace("count", "").strip()
                pretty_col2 = col2.replace("_", " ").replace("Corrected", "").replace("(count)", "").replace("count", "").strip()
                label_lane1 = f"lane {str(lane1)}" if "lane" not in str(lane1) else str(lane1).replace("lane", "lane ")
                label_lane2 = f"lane {str(lane2)}" if "lane" not in str(lane2) else str(lane2).replace("lane", "lane ")
                plots.append({
                    'log_x': df_xy['log_x'],
                    'log_y': df_xy['log_y'],
                    'xlabel': f"{pretty_col1} ({label_lane1}) log10",
                    'ylabel': f"{pretty_col2} ({label_lane2}) log10",
                    'title': f"R$^2$={r2:.4f}",
                    'setname': setname,
                    'sample': pretty_sample
                })
            if len(plots) == 2:
                plot_pairs_with_names.append((pretty_sample, plots))

# === SORT by Gly/Eth then Glucose, and 8 day before 10 day ===
plot_pairs_with_names.sort(key=lambda pair: extract_group_and_day(pair[0]))

def save_plot_grid_pairs(plot_pairs_with_names, outpdf, n_rows=5):
    pairs_per_row = 2
    plots_per_row = pairs_per_row * 2
    num_pairs = len(plot_pairs_with_names)
    num_rows = n_rows

    fig, axes = plt.subplots(num_rows, plots_per_row, figsize=(plots_per_row*6, num_rows*6), squeeze=False)
    plt.subplots_adjust(left=0.07, right=0.98, top=0.88, bottom=0.08, wspace=0.35, hspace=2.0)

    pair_idx = 0
    for row in range(num_rows):
        for pair_col in range(pairs_per_row):
            if pair_idx >= num_pairs:
                for ax in axes[row][(pair_col*2):(pair_col*2+2)]:
                    ax.set_visible(False)
                continue
            sample_label, plots = plot_pairs_with_names[pair_idx]
            for plot_num in range(2):
                ax = axes[row][pair_col*2 + plot_num]
                plot = plots[plot_num]
                ax.scatter(plot['log_x'], plot['log_y'], alpha=0.5)
                ax.set_xlabel(plot['xlabel'], fontsize=15)
                ax.set_ylabel(plot['ylabel'], fontsize=15)
                ax.set_title(plot['title'], fontsize=19)
                allmin = min(plot['log_x'].min(), plot['log_y'].min())
                allmax = max(plot['log_x'].max(), plot['log_y'].max())
                ax.plot([allmin, allmax], [allmin, allmax], ls='--', color='gray')
                ax.tick_params(labelsize=13)
            pair_idx += 1

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])

    fig.suptitle("PCR replicates", fontsize=28)

    column_label_y = 0.95
    column_labels = ["Raw counts", "Chimera corrected"] * pairs_per_row
    for col_idx in range(plots_per_row):
        ax = axes[0][col_idx]
        bbox = ax.get_position()
        x = (bbox.x0 + bbox.x1) / 2
        fig.text(x, column_label_y, column_labels[col_idx], ha="center", va="bottom", fontsize=21, fontweight="bold", color="black")

    sample_label_y_offset = 0.0005
    pair_idx = 0
    for row in range(num_rows):
        for pair_col in range(pairs_per_row):
            if pair_idx >= num_pairs:
                continue
            sample_label = plot_pairs_with_names[pair_idx][0]
            ax_left = axes[row][pair_col*2]
            ax_right = axes[row][pair_col*2+1]
            bbox_left = ax_left.get_position()
            bbox_right = ax_right.get_position()
            x_label = (bbox_left.x0 + bbox_right.x1) / 2
            y_label = bbox_left.y1 + sample_label_y_offset
            fig.text(x_label, y_label, sample_label, ha="center", va="bottom", fontsize=25, fontweight="bold", color="navy")
            pair_idx += 1

    with PdfPages(outpdf) as pdf:
        pdf.savefig(fig)
        plt.close(fig)
    print(f"Saved {min(num_rows*pairs_per_row*2, len(plot_pairs_with_names)*2)} plots to: {outpdf}")

n_rows_per_pdf = 5
n_pairs_per_pdf = n_rows_per_pdf * 2
num_total_pairs = len(plot_pairs_with_names)
num_pdfs = math.ceil(num_total_pairs / n_pairs_per_pdf)

for idx in range(num_pdfs):
    start = idx * n_pairs_per_pdf
    end = start + n_pairs_per_pdf
    pdf_name = f'./figures/fig_S16_part{idx+1}.pdf'
    print(f"Making PDF #{idx+1}: pairs {start} to {end-1}")
    save_plot_grid_pairs(plot_pairs_with_names[start:end], pdf_name, n_rows=n_rows_per_pdf)
