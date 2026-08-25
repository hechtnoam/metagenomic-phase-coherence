#!/usr/bin/env python3
"""Generate the final six-state EM main-text figure.

Panel A: empirical A2424 3x4 base-probability matrices at L=59,60,61.
Panel B: synthetic six-state assignment accuracy for positive vs IID-null data.

Run from repository root.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import t

BASES=("A","T","G","C")
COLORS={"A":"tab:blue","T":"tab:orange","G":"tab:green","C":"tab:red"}

def load_base_probs(path):
    df=pd.read_csv(path,sep="\t")
    req={"length","position_class","base","prob"}
    miss=req-set(df.columns)
    if miss: raise ValueError(f"{path}: missing {sorted(miss)}")
    return df

def load_validation(path):
    df=pd.read_csv(path,sep="\t")
    req={"condition","replicate","state_accuracy_after_symmetry"}
    miss=req-set(df.columns)
    if miss: raise ValueError(f"{path}: missing {sorted(miss)}")
    return df

def ci95(x):
    x=np.asarray(x,dtype=float)
    m=float(np.mean(x))
    if len(x)<2: return m,0.0
    h=float(t.ppf(0.975,len(x)-1)*np.std(x,ddof=1)/np.sqrt(len(x)))
    return m,h

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--base-probs",type=Path,
        default=Path("results/em/A2424_original/base_probs_by_length.tsv"))
    ap.add_argument("--validation",type=Path,
        default=Path("results/em_synthetic_validation/validation_summary.tsv"))
    ap.add_argument("--output-prefix",type=Path,
        default=Path("figures/fig_six_state_em"))
    args=ap.parse_args()

    bp=load_base_probs(args.base_probs)
    val=load_validation(args.validation)

    lengths=[59,60,61]
    for L in lengths:
        if L not in set(bp.length.astype(int)):
            raise ValueError(f"Length {L} absent from {args.base_probs}")

    fig=plt.figure(figsize=(12,7.6))
    gs=fig.add_gridspec(2,6,height_ratios=[1.0,1.15],hspace=0.48,wspace=0.48)
    axesA=[
        fig.add_subplot(gs[0,0:2]),
        fig.add_subplot(gs[0,2:4]),
        fig.add_subplot(gs[0,4:6]),
    ]
    axB=fig.add_subplot(gs[1,1:5])

    width=0.19
    x=np.arange(3)
    offsets=np.linspace(-1.5*width,1.5*width,4)
    for ax,L in zip(axesA,lengths):
        sub=bp[bp.length.astype(int)==L]
        for j,b in enumerate(BASES):
            vals=[]
            for pc in (1,2,3):
                q=sub[(sub.position_class.astype(int)==pc)&(sub.base==b)]
                if len(q)!=1: raise ValueError(f"Expected one row for L={L}, class={pc}, base={b}")
                vals.append(float(q.prob.iloc[0]))
            ax.bar(x+offsets[j],vals,width=width,label=b,color=COLORS[b])
        ax.set_title(f"L={L}",fontsize=11)
        ax.set_xticks(x,["1","2","3"])
        ax.set_xlabel("Position class",fontsize=9)
        ax.set_ylim(0,0.45)
        ax.tick_params(labelsize=8)
    axesA[0].set_ylabel("Estimated base probability",fontsize=9)

    handles=[plt.Rectangle((0,0),1,1,color=COLORS[b]) for b in BASES]
    fig.legend(handles, BASES, loc="upper center", bbox_to_anchor=(0.5, 0.985),
               ncol=4, frameon=False, fontsize=9)
    axesA[0].text(-0.18,1.12,"A",transform=axesA[0].transAxes,
                  fontsize=15,fontweight="bold",va="top")

    conds=[("positive","Known six-state\nstructure"),("iid_null","IID null")]
    xpos=[0,1]
    for xi,(cond,label) in zip(xpos,conds):
        y=val.loc[val.condition==cond,"state_accuracy_after_symmetry"].to_numpy(float)
        # deterministic symmetric jitter: preserves visibility without randomness.
        jitter=np.linspace(-0.065,0.065,len(y))
        axB.scatter(np.full(len(y), xi) + jitter, y, s=22, alpha=0.75, zorder=3)
        m, h = ci95(y)
        axB.errorbar(
            xi, m, yerr=h, fmt="o", markersize=7, capsize=5,
            linewidth=1.5, color="black", markerfacecolor="black",
            markeredgecolor="black", zorder=5
        )
    # Chance accuracy for the balanced six-state synthetic validation:
    # three phase offsets x two orientations = six equally represented states.
    chance = 1 / 6
    axB.axhline(chance, linestyle="--", linewidth=1.2, color="black")
    axB.text(
        1.40, chance + 0.018, "Chance (1/6 = 16.7%)",
        ha="right", va="bottom", fontsize=9
    )
    axB.set_xticks(xpos,[c[1] for c in conds])
    axB.set_ylabel("Six-state assignment accuracy",fontsize=10)
    axB.set_ylim(0,1.0)
    axB.set_xlim(-0.45,1.45)
    axB.tick_params(labelsize=9)
    axB.text(-0.12,1.08,"B",transform=axB.transAxes,
             fontsize=15,fontweight="bold",va="top")

    fig.subplots_adjust(top=0.92, bottom=0.09, left=0.08, right=0.985)

    out=args.output_prefix
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out.with_suffix(".png"),dpi=600,bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"),bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out.with_suffix('.png')}")
    print(f"Wrote {out.with_suffix('.pdf')}")

if __name__=="__main__":
    main()
