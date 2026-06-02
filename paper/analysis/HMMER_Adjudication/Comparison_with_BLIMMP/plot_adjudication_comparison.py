#!/usr/bin/env python3
"""
plot_adjudication_comparison.py

Compare hit_conf values between preprocess_domtblout.py output
and BLIMMP _dk.csv output for matched KOs.

Usage:
    python plot_adjudication_comparison.py \
        --preprocess test_SAMD00000344_preprocess.csv \
        --blimmp new_SAMD00000344__BLIMMP_dk.csv \
        --out adjudication_comparison.png \
        --sample SAMD00000344
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def load_and_merge(preprocess_path, blimmp_path):
    pre = pd.read_csv(preprocess_path)
    bli = pd.read_csv(blimmp_path)

    # Keep only BLIMMP KOs with actual hits
    bli_hits = bli[bli['KO_annotation_conf'] > 0]

    # Merge on KO id
    merged = pre[['KO id', 'score', 'hit_conf']].merge(
        bli_hits[['KO id', 'score', 'KO_annotation_conf', 'kofam_score_threshold', 'flag_is_below_kofam_threshold']],
        on='KO id',
        suffixes=('_pre', '_bli')
    )
    return merged

def plot_comparison(merged, out_path, sample_name=""):
    flag = merged['flag_is_below_kofam_threshold'].astype(str).str.lower() == 'true'
    above = merged[~flag]
    below = merged[flag]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # --- Left panel: score scatter ---
    ax = axes[0]
    ax.scatter(above['score_pre'], above['score_bli'],
               color='#3266ad', alpha=0.6, s=30, label='above KOfam threshold', zorder=3)
    ax.scatter(below['score_pre'], below['score_bli'],
               color='#D85A30', alpha=0.75, s=40, marker='^', label='below KOfam threshold', zorder=4)

    lim = max(merged['score_pre'].max(), merged['score_bli'].max()) * 1.05
    ax.plot([0, lim], [0, lim], 'k--', linewidth=0.8, alpha=0.4, label='y = x')
    ax.set_xlabel('preprocess score', fontsize=12)
    ax.set_ylabel('BLIMMP score', fontsize=12)
    ax.set_title('HMM scores (should be identical)', fontsize=12)
    ax.legend(fontsize=9)

    corr = merged['score_pre'].corr(merged['score_bli'])
    ax.text(0.05, 0.95, f'r = {corr:.4f}', transform=ax.transAxes,
            fontsize=10, va='top', color='#333')

    # --- Right panel: hit_conf scatter ---
    ax = axes[1]
    ax.scatter(above['hit_conf'], above['KO_annotation_conf'],
               color='#3266ad', alpha=0.6, s=30, label='above KOfam threshold', zorder=3)
    ax.scatter(below['hit_conf'], below['KO_annotation_conf'],
               color='#D85A30', alpha=0.75, s=40, marker='^', label='below KOfam threshold', zorder=4)

    ax.plot([0, 1], [0, 1], 'k--', linewidth=0.8, alpha=0.4, label='y = x')
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel('preprocess hit_conf', fontsize=12)
    ax.set_ylabel('BLIMMP KO_annotation_conf', fontsize=12)
    ax.set_title('Annotation confidence', fontsize=12)
    ax.legend(fontsize=9)

    n_diff = (abs(merged['hit_conf'] - merged['KO_annotation_conf']) > 0.01).sum()
    ax.text(0.05, 0.95,
            f'n matched = {len(merged)}\nn differing = {n_diff} (below KOfam threshold)',
            transform=ax.transAxes, fontsize=9, va='top', color='#333')

    title = f'Adjudication comparison: preprocess vs BLIMMP'
    if sample_name:
        title += f'\n{sample_name}'
    fig.suptitle(title, fontsize=13, fontweight='500', y=1.01)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Saved to {out_path}")
    plt.close()

def main():
    ap = argparse.ArgumentParser(description="Compare preprocess_domtblout vs BLIMMP hit_conf")
    ap.add_argument('--preprocess', required=True, help='Path to preprocess_domtblout output CSV')
    ap.add_argument('--blimmp', required=True, help='Path to BLIMMP _dk.csv output')
    ap.add_argument('--out', default='adjudication_comparison.png', help='Output figure path')
    ap.add_argument('--sample', default='', help='Sample name for plot title')
    args = ap.parse_args()

    merged = load_and_merge(args.preprocess, args.blimmp)
    print(f"Matched KOs: {len(merged)}")
    print(f"Score correlation: {merged['score_pre'].corr(merged['score_bli']):.4f}")
    print(f"Conf differences (>0.01): {(abs(merged['hit_conf'] - merged['KO_annotation_conf']) > 0.01).sum()}")

    plot_comparison(merged, args.out, sample_name=args.sample)

if __name__ == '__main__':
    main()
