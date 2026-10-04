#!/usr/bin/env python3
"""Checkpoint figures only; novel-treatment figures are gated on the freeze."""
import csv
from pathlib import Path
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
RESULTS=ROOT/'empirical/neoaves_chr4/results'
FIGURES=ROOT/'empirical/neoaves_chr4/figures'

def main():
    with (RESULTS/'stage4d_figure4_reproduction.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    fig,axes=plt.subplots(1,2,figsize=(9,3.8),sharey=True)
    for ax,treatment,title in zip(axes,('T0','T_PNAS'),('All loci','Published outlier loci removed')):
        r=[x for x in rows if x['treatment']==treatment]
        x=range(len(r)); means=[float(v['mean']) for v in r]; medians=[float(v['median']) for v in r]; se=[float(v['se']) for v in r]
        ax.errorbar(x,means,yerr=se,fmt='o',ms=2.8,lw=.7,capsize=1.5,color='#2040DE',label='mean ± SE')
        ax.scatter(x,medians,s=9,color='#20A060',label='median',zorder=3)
        ax.axhline(0,color='#B02020',lw=.8)
        ax.set_title(title,fontsize=10);ax.set_xlim(-1,len(r));ax.set_xticks([]);ax.grid(axis='x',alpha=.15)
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Δ quartet score (S2024 − J2014)')
    axes[0].legend(frameon=False,fontsize=8,loc='lower left')
    fig.supxlabel('Original Figure-4 taxon-removal order',fontsize=9)
    fig.tight_layout()
    for ext in ('pdf','png'):
        fig.savefig(FIGURES/f'stage4d_figure_A_original_figure4_reproduction.{ext}',dpi=300,bbox_inches='tight')
    plt.close(fig)

if __name__=='__main__':main()
