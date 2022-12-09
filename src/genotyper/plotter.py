from typing import Optional, Tuple

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection


def plot_complex_repeats(df, cols, alleles, img_path: str):
    fig, axes = plt.subplots(nrows=len(cols), ncols=1, figsize=(8, 6*len(cols)))
    for idx, col in enumerate(cols):
        ex_df = df
        val1 = alleles[0][idx]
        val2 = alleles[2][idx]
        name = 'repeat numbers: '+str(val1)+','+str(val2)
        ex_df[name] = ''
        axes[idx] = sns.violinplot(data=ex_df, x=name, y=col, hue='allele', orient='vertical',
                                   split=False, scale='count', whis=np.inf, inner=None, ax=axes[idx])
        plt.setp(axes[idx].collections, alpha=.3)
        first = [r for r in axes[idx].get_children() if type(r) == PolyCollection]
        c1 = first[0].get_facecolor()[0]
        c2 = first[1].get_facecolor()[0]
        if val1 != '-':
            axes[idx].axhline(y=val1, color=c1, linestyle='--')
        if val2 != '-':
            axes[idx].axhline(y=val2, color=c2, linestyle='--')
        axes[idx] = sns.stripplot(data=ex_df, x=name, y=col, hue='allele', orient='vertical',
                                  dodge=True, size=6, alpha=0.8, jitter=0.3, ax=axes[idx])
        axes[idx].get_legend().remove()

    plt.savefig(img_path, bbox_inches='tight', format='svg')
    plt.close()


def plot_clustering_preds(img_path: str, vals, bvals, alleles: Tuple[int, int], alleles_bc: Optional[Tuple[int, int]]):
    if bvals is not None:
        _, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 7), sharey=True, gridspec_kw={'width_ratios': [1, 1]})
        axes[0] = sns.violinplot(data=vals, orient='vertical', split=True, scale='count',
                                 whis=np.inf, inner=None, color='.8', width=0.8, scale_hue=False, ax=axes[0])
        axes[0] = sns.stripplot(data=vals, orient='vertical', size=6, alpha=0.6,
                                jitter=0.3, palette='bright', ax=axes[0])
        axes[0].title.set_text('WarpSTR')
        axes[0].set(xlabel=f'Alleles of length {alleles[0]},{alleles[1]}')
        axes[0].set(ylabel='Filtered predictions')
        axes[1] = sns.violinplot(data=bvals, orient='vertical', split=True, scale='count',
                                 whis=np.inf, inner=None, color='.8', width=0.8, scale_hue=False, ax=axes[1])
        axes[1] = sns.stripplot(data=bvals, orient='vertical', size=6, alpha=0.6,
                                jitter=0.3, palette='bright', ax=axes[1])
        axes[1].title.set_text('Basecalled sequences')
        axes[1].set(xlabel=f'Alleles of length {alleles_bc[0]},{alleles_bc[1]}')
        axes[1].set(ylabel='Filtered predictions')
        plt.savefig(img_path, bbox_inches='tight', format='svg')
        plt.close()
        return

    fig = plt.figure(figsize=(16, 7))
    fig.suptitle('WarpSTR', fontsize=20)
    plt.xlabel(f'Alleles of length {alleles[0]},{alleles[1]}', fontsize=18)
    plt.ylabel('Filtered predictions', fontsize=16)
    sns.violinplot(data=vals, orient='vertical', split=True, scale='count',
                   whis=np.inf, inner=None, color='.8', width=0.8, scale_hue=False)
    sns.stripplot(data=vals, orient='vertical', size=6, alpha=0.6, jitter=0.3, palette='bright')
    plt.savefig(img_path, bbox_inches='tight', format='svg')
    plt.close()
