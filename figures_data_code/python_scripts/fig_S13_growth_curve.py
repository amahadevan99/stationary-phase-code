import pandas as pd
import matplotlib.pyplot as plt
import os

df = pd.read_csv('../python_data/growth_curve_data.csv')
df.columns = df.columns.str.strip()

fig, ax = plt.subplots(figsize=(8, 6))

for curve, grp in df.groupby('growth curve'):
    sub = grp.dropna(subset=['density from plates'])
    ax.scatter(sub['elapsed time'], sub['density from plates'],
               label=f'Curve {curve}', s=40,color='black')
ax.axvline(48,linestyle='--',color='black')
ax.set_xlabel('time (hours)', fontsize=16)
ax.set_ylabel('CFU/mL', fontsize=16)
ax.tick_params(axis='both', labelsize=13)
# ax.legend(fontsize=13)
ax.set_yscale('log')

plt.tight_layout()
outfile = os.path.join('figures', 'fig_S13.pdf')
plt.savefig(outfile, bbox_inches='tight', format='pdf')
plt.close()
print(f"PDF saved to: {os.path.abspath(outfile)}")
