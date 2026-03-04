import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.colors as mcolors
import random
import re
from matplotlib.backends.backend_pdf import PdfPages

# --- Input Folders ---
traj_folder = '../fitseq2_input'
t0_folder = '../T0_counts_aggregated'
fitness_folder = '../fitseq2_output'

# --- Output Folder ---
# Create figures folder in the current working directory
output_folder = 'figures'
os.makedirs(output_folder, exist_ok=True)

# --- Constants ---
EXTINCT_Y = 1e-9
NOT_DETECTED_Y = 1e-8
MAX_FREQ = 1.2

# --- Helper Functions ---
def assign_categories(original, extinct_y, not_detected_y):
    num_tp = len(original)
    categories = np.zeros(num_tp)
    extinct_flag = False
    for i in range(num_tp):
        if extinct_flag:
            categories[i] = extinct_y
        elif original[i] == 0:
            if i == num_tp - 1:
                categories[i] = extinct_y
                extinct_flag = True
            elif np.all(original[i+1:] == 0):
                categories[i:] = extinct_y
                extinct_flag = True
                break
            else:
                categories[i] = not_detected_y
        else:
            categories[i] = original[i]
    return categories

def find_fitness_file(run_id, fitness_folder):
    target_suffix = '_FitSeq2_Result.csv'
    for fname in os.listdir(fitness_folder):
        if run_id in fname and fname.endswith(target_suffix):
            return os.path.join(fitness_folder, fname)
    return None

def parse_traj_filename(fname):
    m = re.match(r'(M3|M05)_(\d+)_day_r(\d+)\.csv', fname)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return None, None, None

# --- Main Logic ---
traj_files = [f for f in os.listdir(traj_folder) if re.match(r'(M3|M05)_\d+_day_r\d+\.csv', f)]

# Panel layouts
m3_days = [1,2,3,5]
m05_days = [2,4,6,8,10]
reps = [1,2,3]
m3_grid = {(d, r): None for d in m3_days for r in reps}
m05_grid = {(d, r): None for d in m05_days for r in reps}

for fname in traj_files:
    exp_type, day, rep = parse_traj_filename(fname)
    if exp_type == "M3" and day in m3_days and rep in reps:
        m3_grid[(day, rep)] = fname
    elif exp_type == "M05" and day in m05_days and rep in reps:
        m05_grid[(day, rep)] = fname

def plot_panel(grid, days, reps, pdf_path, exp_prefix, top_label):
    nrows = len(days)
    ncols = len(reps)
    row_labels = [f"{day} day" if day == 1 else f"{day} days" for day in days]

    fig_width = 16
    fig_height = 3.2*nrows + 1.2

    orig_left = 0.17
    orig_right = 0.93
    orig_grid_frac = orig_right - orig_left
    orig_grid_inches = (4.5 * ncols + 1.5) * orig_grid_frac
    new_left = orig_left
    new_grid_frac = orig_grid_inches / fig_width
    new_right = new_left + new_grid_frac

    with PdfPages(pdf_path) as pdf:
        fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), squeeze=False)
        plt.subplots_adjust(left=new_left, right=new_right, top=0.87, bottom=0.13, hspace=0.7, wspace=0.6)

        fitness_min, fitness_max = None, None
        for day in days:
            for rep in reps:
                fname = grid.get((day, rep), None)
                if fname is not None:
                    run_id = fname.replace('.csv', '')
                    fitness_path = find_fitness_file(run_id, fitness_folder)
                    try:
                        fitness_df_for_range = pd.read_csv(fitness_path)
                        if 'Error_Fitness' in fitness_df_for_range.columns:
                           fitness_df_for_range = fitness_df_for_range[fitness_df_for_range['Error_Fitness'] < 5]
                        if 'Fitness_Per_Cycle' not in fitness_df_for_range.columns: fitness_df_for_range.columns = fitness_df_for_range.columns.str.strip()
                        if 'Fitness_Per_Cycle' in fitness_df_for_range.columns:
                            vals = fitness_df_for_range['Fitness_Per_Cycle'].dropna().values
                            if len(vals):
                                this_min, this_max = np.min(vals), np.max(vals)
                                if fitness_min is None or this_min < fitness_min: fitness_min = this_min
                                if fitness_max is None or this_max > fitness_max: fitness_max = this_max
                    except Exception:
                        continue
        if fitness_min is None or fitness_max is None: fitness_min, fitness_max = 0, 1
        cmap, norm = plt.get_cmap('viridis'), mcolors.Normalize(vmin=fitness_min, vmax=fitness_max)

        for irow, day in enumerate(days):
            for icol, rep in enumerate(reps):
                ax = axes[irow, icol]
                fname = grid.get((day, rep), None)
                if fname is not None:
                    file_path = os.path.join(traj_folder, fname)
                    run_id = fname.replace('.csv', '')
                    fitness_path = find_fitness_file(run_id, fitness_folder)
                    t0_file = None
                    exp_type, day_num, _ = parse_traj_filename(fname)
                    if exp_type == "M3":
                        if day_num == 1: t0_file = 'FOS3_ROS1_agg.csv'
                        elif day_num in [2, 3]: t0_file = 'FOS5_ROS9_agg.csv'
                        elif day_num == 5: t0_file = 'FOS9_ROS7_agg.csv'
                    elif exp_type == "M05": t0_file = 'FOS9_ROS7_agg.csv'
                    else: t0_file = 'FOS9_ROS7_agg.csv'
                    t0_path = os.path.join(t0_folder, t0_file)

                    try:
                        traj_df, t0_df, fitness_df = pd.read_csv(file_path, header=None), pd.read_csv(t0_path, header=None), pd.read_csv(fitness_path)
                    except Exception as e:
                        ax.text(0.5, 0.5, f"File error!\n{e}", ha="center", va="center", fontsize=8); ax.axis('off'); continue

                    error_threshold, error_column_name = 5, 'Error_Fitness'
                    if error_column_name in fitness_df.columns:
                        fitness_df = fitness_df[fitness_df[error_column_name] < error_threshold].copy()
                    else:
                        print(f"Warning: '{error_column_name}' column not found. Cannot filter by error for {fname}.")

                    if fitness_df.empty: ax.text(0.5, 0.5, "No variants left\nafter error filter", ha="center", va="center", fontsize=8); ax.axis('off'); continue

                    combined_df = pd.concat([t0_df, traj_df], axis=1)
                    freq_df_full = combined_df.div(combined_df.sum(axis=0), axis=1)
                    common_indices = fitness_df.index.intersection(freq_df_full.index)
                    fitness_df, freq_df, combined_df = fitness_df.loc[common_indices], freq_df_full.loc[common_indices], combined_df.loc[common_indices]

                    if freq_df.empty: ax.text(0.5, 0.5, "No matching variants", ha="center", va="center", fontsize=8); ax.axis('off'); continue

                    n_cols_orig = combined_df.shape[1]
                    COL_GENERATIONS = [i * 8 for i in range(n_cols_orig)]
                    mask_empty = (combined_df.isna() | (combined_df == 0)).all(axis=0)

                    combined_df = combined_df.loc[:, ~mask_empty]
                    freq_df = freq_df.loc[:, ~mask_empty]

                    indices_kept = [i for i, keep in enumerate(~mask_empty) if keep]
                    time_labels = [COL_GENERATIONS[i] for i in indices_kept]

                    if 'Fitness_Per_Cycle' not in fitness_df.columns: fitness_df.columns = fitness_df.columns.str.strip()
                    if 'Fitness_Per_Cycle' not in fitness_df.columns: ax.text(0.5, 0.5, "No Fitness_Per_Cycle", ha="center", va="center", fontsize=8); ax.axis('off'); continue

                    fitness_vals = fitness_df['Fitness_Per_Cycle']

                    top_idx = fitness_vals.nlargest(100).index.tolist()
                    remaining_idx = freq_df.index.difference(top_idx).tolist()
                    rand_idx = random.sample(remaining_idx, min(100, len(remaining_idx)))
                    plot_idx = top_idx + rand_idx
                    random.shuffle(plot_idx)

                    for idx in plot_idx:
                        trajectory = freq_df.loc[idx].values.astype(float)
                        values = assign_categories(trajectory, EXTINCT_Y, NOT_DETECTED_Y)
                        values = np.clip(values, 0, MAX_FREQ)
                        fitness_value = fitness_vals.loc[idx]
                        color = cmap(norm(fitness_value))
                        ax.plot(time_labels, values, linewidth=1, alpha=0.8, color=color)

                    ax.set_yscale('log'); ax.set_ylim(bottom=1e-10, top=MAX_FREQ); ax.set_xticks(time_labels)
                    yticks, ylabels = [EXTINCT_Y, NOT_DETECTED_Y, 1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1.0], ['extinct', 'not detected', '1e-7','1e-6','1e-5','1e-4','0.001','0.01','0.1','1']
                    ax.set_yticks(yticks); ax.set_yticklabels(ylabels)

                    col_sums, nonzero_counts = combined_df.sum(axis=0).values.astype(int), (combined_df > 0).sum(axis=0).values.astype(int)
                    n_tp, t0_x, counts_start, counts_end = len(time_labels), -0.25, 0.18, 1.08
                    for j, x in enumerate(time_labels):
                        if j == 0:
                            x_pos = t0_x
                            ax.text(x_pos, -0.29, f'Sum: {col_sums[j]:,}', ha='left', va='top', fontsize=7, transform=ax.transAxes)
                            ax.text(x_pos, -0.39, f'Active: {nonzero_counts[j]:,}', ha='left', va='top', fontsize=7, transform=ax.transAxes)
                        else:
                            frac = (j-1) / (n_tp-2) if n_tp > 2 else 0
                            x_pos = counts_start + (counts_end - counts_start) * frac
                            ax.text(x_pos, -0.29, f'{col_sums[j]:,}', ha='center', va='top', fontsize=7, transform=ax.transAxes)
                            ax.text(x_pos, -0.39, f'{nonzero_counts[j]:,}', ha='center', va='top', fontsize=7, transform=ax.transAxes)
                else:
                    ax.axis('off')

        for icol, rep in enumerate(reps): axes[0, icol].set_title(f"Replicate {rep}", fontsize=14, pad=18)
        for ax_row in axes:
            for ax in ax_row: ax.set_xlabel(''); ax.set_ylabel('')
        for irow, label in enumerate(row_labels):
            fig.text(0.09, axes[irow,0].get_position().get_points().mean(axis=0)[1], label, fontsize=14, va='center', ha='left', rotation='vertical', bbox=dict(boxstyle="round,pad=0.3", fc="w", alpha=0.7))
        fig.text(0.54, 0.06, 'Generation', fontsize=14, ha='center')
        fig.text(0.06, 0.5, 'Frequency', fontsize=14, va='center', rotation='vertical')
        cbar_ax = fig.add_axes([new_right + 0.033, 0.15, 0.018, 0.7])

        # --- FIX IS HERE: Broke the single line into two ---
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = fig.colorbar(sm, cax=cbar_ax)

        cbar.set_label('Fitness per cycle')
        fig.text(0.5, 0.93, top_label, fontsize=18, ha='center', va='center', fontweight='bold')

        pdf.savefig()
        plt.close()
        print(f"Saved panel PDF: {pdf_path}")

# --- Execute Plotting ---
m3_pdf = os.path.join(output_folder, 'fig_S2_glucose.pdf')
plot_panel(m3_grid, m3_days, reps, m3_pdf, 'M3', top_label="Glucose")

m05_pdf = os.path.join(output_folder, 'fig_S2_glyeth.pdf')
plot_panel(m05_grid, m05_days, reps, m05_pdf, 'M05', top_label="Gly/Eth")

print(f"✅ Saved all panel PDFs to: {output_folder}")
