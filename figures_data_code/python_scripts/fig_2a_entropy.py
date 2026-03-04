import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

exp_name_list = ['2_day_r2', '4_day_r1', '4_day_r2', '6_day_r1', '6_day_r3',
                 '8_day_r2', '8_day_r3b', '10_day_r1', '10_day_r2']

bc_dict = {}
xrange_dict = {}

for exp_name in exp_name_list:
    bc_dict[exp_name] = pd.read_csv(
        f'../fitmut_input/filtered/{exp_name}_counts.csv', header=None).values
    xrange_dict[exp_name] = pd.read_csv(
        f'../fitmut_input/params0/{exp_name}_timepoints.csv', header=None)[0].values

cmap = mpl.cm.get_cmap('plasma')
norm = mpl.colors.Normalize(vmin=2, vmax=10.5)

day_labels_seen = set()

fig, ax = plt.subplots(figsize=(8, 6))

for exp_name in exp_name_list:
    if exp_name == '4_day_r2':
        continue
    days = int(exp_name.split('_')[0])
    rep = int(exp_name.split('_')[-1].strip('r').strip('b')) - 1

    bc_counts = bc_dict[exp_name].T  # shape: (timepoints, barcodes)
    total_reads = bc_counts.sum(axis=1)
    bc_freq = bc_counts / total_reads[:, np.newaxis]

    log_freq = np.ma.log(bc_freq).filled(0)
    shannon_entropy = np.sum(-bc_freq * log_freq, axis=1)

    xrange = xrange_dict[exp_name]

    label = f'{days}' if days not in day_labels_seen else None
    day_labels_seen.add(days)

    ax.plot(xrange, shannon_entropy, label=label,
            color=cmap(norm(days)), linewidth=4,
            marker='o', markersize=7)

ax.set_xlabel('Generation', fontsize=16)
ax.set_ylabel('Shannon diversity', fontsize=16)
ax.tick_params(axis='both', labelsize=13)
ax.legend(title='Days between\ntransfers', fontsize=12, title_fontsize=13)

plt.tight_layout()
outfile = os.path.join('figures', 'fig_2a.pdf')
plt.savefig(outfile, bbox_inches='tight', format='pdf')
plt.close()
print(f"PDF saved to: {os.path.abspath(outfile)}")
