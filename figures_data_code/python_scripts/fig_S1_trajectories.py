import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
import os
import re

counts_folder = '../fitmut_input/filtered'
fitness_folder = '../fitmut_output_filtered'

EXTINCT_Y = 1e-9
NOT_DETECTED_Y = 1e-8
MAX_FREQ = 1.2

def assign_categories(raw, freq, extinct_y, not_detected_y):
    num_tp = len(raw)
    out = np.zeros(num_tp)
    extinct_flag = False
    for i in range(num_tp):
        if extinct_flag:
            out[i] = extinct_y
        elif raw[i] <= 1:
            if i == num_tp - 1:
                out[i] = extinct_y
                extinct_flag = True
            elif np.all(raw[i+1:] <= 1):
                out[i:] = extinct_y
                extinct_flag = True
                break
            else:
                out[i] = not_detected_y
        else:
            out[i] = freq[i]
    return out

def extract_info_from_filename(counts_file):
    m = re.match(r"(.+?)_([0-9]+)_day_r([0-9a-zA-Z]+)_counts\.csv", counts_file)
    if m:
        carbon = m.group(1).replace('_', '/')
        days = int(m.group(2))
        rep = m.group(3)
        return carbon, days, rep
    m = re.match(r"([0-9]+)_day_r([0-9a-zA-Z]+)_counts\.csv", counts_file)
    if m:
        carbon = ""
        days = int(m.group(1))
        rep = m.group(2)
        return carbon, days, rep
    return None, None, None

# ----------- STAR ANNOTATION LOGIC -----------
def T_label_to_gen(tlabel):
    if tlabel.startswith("T") and tlabel[1:].isdigit():
        return int(tlabel[1:]) * 8
    return None

star_map = {
    "10 day #1": ["T2", "T5", "T3", "T3", "T4"],
    "10 day #2": ["T2", "T3", "T4", "T5"],
    "2 day #2": ["T13", "T14", "T12"],
    "4 day #1": ["T4", "T5", "T6", "T7", "T8"],
    "6 day #1": ["T4", "T5", "T6"],
    "6 day #3": ["T4", "T5", "T6"],
    "8 day #2": ["T8", "T9", "T3", "T4", "T5"],
    "8 day #3b": ["T8", "T9"]
}
star_map_gens = {k: [T_label_to_gen(t) for t in v if T_label_to_gen(t) is not None] for k, v in star_map.items()}

def title_to_key(title):
    # Example: "10 days replicate 1" -> "10 day #1"
    match = re.match(r'(\d+)\s+days?\s+replicate\s*([0-9a-zA-Z]+)', title)
    if match:
        return f"{match.group(1)} day #{match.group(2)}"
    match = re.match(r'(\d+)\s+day#([0-9a-zA-Z]+)', title)
    if match:
        return f"{match.group(1)} day#{match.group(2)}"
    # fallback: strip whitespace
    return title.strip()

counts_pat = re.compile(r"_day_r[0-9a-zA-Z]+_counts\.csv$")
counts_files = [f for f in os.listdir(counts_folder) if counts_pat.search(f)]

plotdata = []
titles = []
order_keys = []

for counts_file in counts_files:
    base = counts_file.replace("_counts.csv", "")
    mutseq_file = f"{base}_MutSeq_Result.csv"
    counts_path = os.path.join(counts_folder, counts_file)
    mutseq_path = os.path.join(fitness_folder, mutseq_file)
    if not os.path.exists(mutseq_path):
        print(f"⚠️ {mutseq_file} not found, skipping {counts_file}")
        continue

    carbon, days, rep = extract_info_from_filename(counts_file)
    title_parts = []
    if carbon and carbon.strip():
        title_parts.append(carbon.strip())
    title_parts.append(f"{days} days")
    title_parts.append(f"replicate {rep}")
    plot_title = " ".join(title_parts)

    traj_df = pd.read_csv(counts_path, header=None)
    n_timepoints = traj_df.shape[1]
    generation_labels = [i * 8 for i in range(n_timepoints)]

    fitness_df = pd.read_csv(mutseq_path)
    fitness = fitness_df['Fitness'] * 8  # Per cycle

    assert len(traj_df) == len(fitness), f"Trajectory and fitness length mismatch for {cell_file}"

    total_counts = traj_df.sum(axis=1)
    top_idx = total_counts.nlargest(2000).index
    traj_df_top = traj_df.loc[top_idx]
    fitness_top = fitness.loc[top_idx].reset_index(drop=True)
    freq_df_top = traj_df_top.div(traj_df_top.sum(axis=0), axis=1)

    rep_digits = re.findall(r'\d+', str(rep))
    rep_key = int(rep_digits[0]) if rep_digits else 0
    order_keys.append((days, rep_key, str(rep)))
    plotdata.append((traj_df_top, freq_df_top, fitness_top, generation_labels))
    titles.append(plot_title)

sorting = sorted(zip(order_keys, plotdata, titles))
plotdata = [p for _, p, t in sorting]
titles = [t for _, p, t in sorting]

n_plots = len(plotdata)
grid_n = int(np.ceil(np.sqrt(n_plots)))

all_fitness_vals = pd.concat([d[2] for d in plotdata], ignore_index=True)
fitness_min = all_fitness_vals.min()
fitness_max = all_fitness_vals.max()

cmap = plt.get_cmap('viridis')
norm = mcolors.Normalize(vmin=fitness_min, vmax=fitness_max)
# timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_pdf_path = os.path.join('./figures/', 'fig_S1.pdf')

with PdfPages(output_pdf_path) as pdf:
    fig, axes = plt.subplots(grid_n, grid_n, figsize=(5*grid_n+1.5, 4*grid_n), sharey=True)  # no sharex!
    axes = axes.flatten()
    for i in range(grid_n * grid_n):
        ax = axes[i]
        if i < len(plotdata):
            traj_df_top, freq_df_top, fitness_top, generation_labels = plotdata[i]
            title = titles[i]
            for idx in freq_df_top.index:
                raw_counts = traj_df_top.loc[idx].values
                freq = freq_df_top.loc[idx].values
                fitness_value = fitness_top.iloc[freq_df_top.index.get_loc(idx)]
                color = cmap(norm(fitness_value))
                values = assign_categories(raw_counts, freq, EXTINCT_Y, NOT_DETECTED_Y)
                values = np.clip(values, 0, MAX_FREQ)
                ax.plot(generation_labels, values, color=color, alpha=0.8, linewidth=1)
            ax.set_title(title, fontsize=13)
            ax.set_yscale('log')
            ax.set_ylim(bottom=EXTINCT_Y/10, top=MAX_FREQ)
            # x-limits and ticks are per-plot
            x_min = min(generation_labels)
            x_max = max(generation_labels)
            ax.set_xlim(x_min, x_max)
            # Show ticks at every 20 between min/max
            ax.set_xticks([x for x in range(x_min, x_max+1, 20) if x_min <= x <= x_max])
            ax.axhline(EXTINCT_Y, color='gray', linestyle=':', linewidth=1.2, label='extinct')
            ax.axhline(NOT_DETECTED_Y, color='black', linestyle=':', linewidth=1.2, label='not detected')
            ax.tick_params(axis='x', which='both', labelbottom=True)
            ax.tick_params(axis='y', which='both', labelleft=True)
            ax.set_yticks(
                [EXTINCT_Y, NOT_DETECTED_Y, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0])
            ax.set_yticklabels(
                ['extinct', 'not detected', '1e-7', '1e-6', '1e-5', '1e-4',
                 '0.001', '0.01', '0.1', '1'])

            # ======== ADD STAR MARKERS (FOR THIS PANEL, FROM MAP) =========
            key = title_to_key(title)
            star_xs = star_map_gens.get(key, [])
            if star_xs:
                star_y = ax.get_ylim()[0] * 2
                ax.plot(star_xs, [star_y] * len(star_xs), marker='*', color='red', markersize=14,
                        linestyle='None', zorder=10)
        else:
            ax.axis('off')

    fig.text(0.5, 0.08, 'Generation', ha='center', va='center', fontsize=15)
    fig.text(0.04, 0.5, 'Frequency', ha='center', va='center', fontsize=15, rotation='vertical')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.017, 0.7])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Fitness per cycle', fontsize=13)

    fig.tight_layout(rect=[0.07, 0.09, 0.90, 0.98])
    pdf.savefig(fig)
    plt.close(fig)

print(f"✅ PDF saved to {output_pdf_path}")
