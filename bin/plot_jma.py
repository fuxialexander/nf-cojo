# %%
import pandas as pd
import seaborn as sns 
from sys import argv
import numpy as np
import matplotlib.pyplot as plt
# %%
# argv[1] = "../results/test_gwas_merged.jma.cojo"
d = pd.read_csv(argv[1], sep='\t', index_col=1)

# %%
sns.distplot(d.LD_r**2)

# %%
fig, axs = plt.subplots(1, 2, figsize=(6,3))
sns.scatterplot(-np.log(d.p), -np.log(d.pJ), pd.Series(d.LD_r**2, name='LD R2'), data = d, s=20, ax=axs[0])
lims = [
    np.min([axs[0].get_xlim(), axs[0].get_ylim()]),  # min of both axes
    np.max([axs[0].get_xlim(), axs[0].get_ylim()]),  # max of both axes
]

# now plot both limits against eachother
axs[0].plot(lims, lims, 'k-', alpha=0.75, zorder=0)
axs[0].set_aspect('equal')
axs[0].set_xlim(lims)
axs[0].set_ylim(lims)
axs[0].set_xlabel("Original -log P")
axs[0].set_ylabel("Joint model -log P")

sns.distplot(d.LD_r**2, ax=axs[1])
axs[1].set_xlabel("LD R2")
axs[1].set_aspect(1./axs[1].get_data_ratio())
fig.savefig(argv[1]+'.png')
# %%
