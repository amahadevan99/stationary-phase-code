import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os 

infile = '../python_data/evo_transfer_record.csv'
df = pd.read_csv(infile)

def explode_replicates(row):
    reps = str(row['replicate #']).split(',')
    reps = [r.strip() for r in reps]
    return [dict(row, **{'replicate #': rep}) for rep in reps]

rows = []
for _, row in df.iterrows():
    if ',' in str(row['replicate #']):
        rows.extend(explode_replicates(row))
    else:
        row_dict = row.to_dict()
        row_dict['replicate #'] = str(row_dict['replicate #']).strip()
        rows.append(row_dict)
df_expanded = pd.DataFrame(rows)

df_expanded['% diploid'] = pd.to_numeric(df_expanded['% diploid'], errors='coerce')
df_expanded['timepoint'] = pd.to_numeric(df_expanded['timepoint'], errors='coerce')
df_expanded['replicate #'] = df_expanded['replicate #'].astype(str)
df_expanded['Generations'] = df_expanded['timepoint'] * 8

panel_order = ['2 day', '4 day', '6 day', '8 day', '10 day']
existing_panels = [p for p in panel_order if p in df_expanded['evolution condition'].unique()]

num_panels = len(existing_panels)
nrows = 2
ncols = 3

fig, axes = plt.subplots(nrows, ncols, figsize=(16, 8), sharey=True)
axes = axes.flatten()

for ax, cond in zip(axes, existing_panels):
    cond_df = df_expanded[df_expanded['evolution condition'] == cond]
    cond_df = cond_df.dropna(subset=['% diploid', 'Generations'])

    # Exclude replicate #1 in '8 day' plot
    if cond == '8 day':
        cond_df = cond_df[cond_df['replicate #'] != '1']

    for rep, rep_df in cond_df.groupby('replicate #'):
        rep_df_sorted = rep_df.sort_values('Generations')
        ax.plot(rep_df_sorted['Generations'], rep_df_sorted['% diploid'],
                marker='o', label=f'Replicate {rep}', linestyle='-')
    ax.set_title(f'Gly/Eth {cond}', fontsize=14)
    ax.set_xlabel("Generations", fontsize=12)
    ax.tick_params(axis='both', labelsize=11)
    ax.legend(fontsize=10)
    ax.grid(True)

# Set y-label only on left plots (first col)
for idx, ax in enumerate(axes):
    if idx % ncols == 0:
        ax.set_ylabel("% diploid", fontsize=12)

# Hide unused subplot(s)
if num_panels < len(axes):
    for ax in axes[num_panels:]:
        ax.set_visible(False)

plt.tight_layout()

# Create figures folder in the current working directory
pdf_outfile = os.path.join("figures", "fig_S14.pdf")

with PdfPages(pdf_outfile) as pdf:
    pdf.savefig(fig)
plt.close(fig)
print(f"Panel PDF saved to: {os.path.abspath(pdf_outfile)}")
