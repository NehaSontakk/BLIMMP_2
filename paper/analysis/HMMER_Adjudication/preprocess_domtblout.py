#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Process a HMMER .domtblout file for ATB dataset

Pipeline:
  1. Parse domtblout columns (auto-detect KO column)
  2. Use full_score / full_Evalue (not domain-level)
  3. Compute strand, hmm_len, per-hit HMM coverage
  4. Compute HMM union coverage per (strand, target, KO) group via Numba
  5. Assign strand-aware overlap groups via Numba clustering
  6. Compute per-position softmax within each overlap group
  7. Compute noise-floor hit confidence and pick one winner per overlap group
  8. Write results to CSV

Usage:
  python process_domtblout.py path/to/file.domtblout \
      --out grouped_hits_winners.csv \
      --e-threshold 1e-4
"""

import sys
import re
import argparse
import pandas as pd
import numpy as np
from numba import njit


KO_RE = re.compile(r'^K\d{5}$')


#  Numba kernels 

@njit
def _hmm_union_len_per_group(gids, starts, ends, n_groups):
    """Compute the union of HMM coordinate intervals per group id."""
    covered = np.zeros(n_groups, dtype=np.int64)

    cur_gid = gids[0]
    cur_s   = starts[0]
    cur_e   = ends[0]

    for i in range(1, len(starts)):
        g = gids[i]
        s = starts[i]
        e = ends[i]

        if g != cur_gid:
            covered[cur_gid] += (cur_e - cur_s + 1)
            cur_gid = g
            cur_s   = s
            cur_e   = e
        else:
            if s <= cur_e:
                if e > cur_e:
                    cur_e = e
            else:
                covered[cur_gid] += (cur_e - cur_s + 1)
                cur_s = s
                cur_e = e

    covered[cur_gid] += (cur_e - cur_s + 1)
    return covered


@njit
def _assign_groups_numba(starts, ends, frac_thresh):
    """Assign overlap group ids by fractional overlap against the shorter segment."""
    n = len(starts)
    g_st = np.empty(n, dtype=np.float64)
    g_en = np.empty(n, dtype=np.float64)
    gcount = 0
    grp_ids = np.empty(n, dtype=np.int32)

    for i in range(n):
        s = starts[i]
        e = ends[i]
        assigned = False
        for g in range(gcount - 1, -1, -1):
            gs = g_st[g]
            ge = g_en[g]
            if s > ge:
                break
            overlap = min(e, ge) - max(s, gs)
            if overlap <= 0.0:
                continue
            short_len = (e - s) if (e - s) < (ge - gs) else (ge - gs)
            if (overlap / short_len) >= frac_thresh:
                grp_ids[i] = g + 1
                if s < gs:
                    g_st[g] = s
                if e > ge:
                    g_en[g] = e
                assigned = True
                break
        if not assigned:
            g_st[gcount] = s
            g_en[gcount] = e
            gcount += 1
            grp_ids[i] = gcount

    return grp_ids


# Parsing 

def process_domtblout(path: str) -> pd.DataFrame:
    """Read a .domtblout and return a normalized DataFrame (BLIMMP-aligned)."""
    cols = [
        'target name', 'target_accession', 'tlen', 'query_name',
        'query_accession', 'qlen', 'full_Evalue', 'full_score',
        'full_bias', 'n_domains', 'of_domains', 'c_Evalue',
        'i_Evalue', 'i_score', 'i_bias', 'hmm from', 'hmm to',
        'ali from', 'ali to', 'env from', 'env to', 'acc'
    ]

    usecols = [
        'target name', 'query_name',
        'hmm from', 'hmm to', 'tlen',
        'ali from', 'ali to', 'qlen',
        'full_score', 'full_Evalue',
        'i_score', 'i_Evalue'
    ]

    df = pd.read_csv(
        path,
        comment='#',
        header=None,
        names=cols,
        usecols=list(range(22)),
        sep=r"\s+",
        engine='c',
        low_memory=False,
        memory_map=True
    )

    df = df[usecols].copy()

    # Auto-detect which column contains KO IDs (KNNNNN)
    q_matches = df['query_name'].astype(str).str.fullmatch(KO_RE.pattern, na=False)
    t_matches = df['target name'].astype(str).str.fullmatch(KO_RE.pattern, na=False)

    fq = float(q_matches.mean()) if len(df) else 0.0
    ft = float(t_matches.mean()) if len(df) else 0.0
    print(f"Detecting KO column ... "
          f"(query_name match frac={fq:.3f}, target name match frac={ft:.3f})",
          file=sys.stderr)

    if (fq > 0 or ft > 0) and (fq >= ft):
        df = df.rename(columns={'query_name': 'KO id', 'target name': 'target name'})
        df['hmm_len'] = df['qlen']
    elif ft > 0:
        df = df.rename(columns={'target name': 'KO id', 'query_name': 'target name'})
        df['hmm_len'] = df['tlen']
    else:
        raise ValueError("Could not detect KO ID column (no values like KNNNNN).")

    # Use full-sequence scores (aligned with BLIMMP)
    df = df.rename(columns={'full_score': 'score', 'full_Evalue': 'E-value'})

    # Compute spans and strand, filter zero-length
    df['ali_span'] = (df['ali to'] - df['ali from']).abs()
    df['hmm_span'] = (df['hmm to'] - df['hmm from']).abs()
    df['strand'] = np.where(df['ali to'] >= df['ali from'], '+', '-')

    before = len(df)
    df = df[(df['ali_span'] > 0) & (df['hmm_span'] > 0)].copy()
    print(f"Filtered zero-length hits: kept {len(df)}/{before}", file=sys.stderr)

    df['per_hit_hmm_coverage'] = (df['hmm_span'] + 1) / df['hmm_len']

    df.drop(columns=['ali_span', 'hmm_span'], inplace=True)

    return df


# HMM union coverage 

def compute_hmm_union_coverage(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute HMM union coverage per (strand, target name, KO id) group.
    Adds: group_id, hmm_covered_len, hmm_coverage_fraction.
    """
    df = df.copy()

    df['group_id'] = df.groupby(
        ['strand', 'target name', 'KO id'],
        sort=False
    ).ngroup()

    n_groups = df['group_id'].max() + 1

    starts = df[['hmm from', 'hmm to']].min(axis=1).to_numpy(np.int64)
    ends   = df[['hmm from', 'hmm to']].max(axis=1).to_numpy(np.int64)
    gids   = df['group_id'].to_numpy(np.int64)

    order = np.lexsort((starts, gids))
    starts = starts[order]
    ends   = ends[order]
    gids   = gids[order]

    covered = _hmm_union_len_per_group(gids, starts, ends, n_groups)

    hmm_len_per_group = (
        df.groupby('group_id', sort=False)['hmm_len']
        .first()
        .to_numpy()
    )

    coverage_per_group = covered / hmm_len_per_group

    gid_full = df['group_id'].to_numpy()
    df['hmm_covered_len']       = covered[gid_full]
    df['hmm_coverage_fraction'] = coverage_per_group[gid_full]

    return df


# Overlap grouping 

def cluster_strand(df: pd.DataFrame,
                   from_col: str = "ali from",
                   to_col: str = "ali to",
                   frac_thresh: float = 0.6) -> pd.DataFrame:
    """Single-strand overlap clustering via Numba."""
    starts = df[[from_col, to_col]].min(axis=1).to_numpy(dtype=np.float64)
    ends   = df[[from_col, to_col]].max(axis=1).to_numpy(dtype=np.float64)
    orig_idx = df.index.to_numpy()

    #order = np.argsort(starts, kind="mergesort")
    ends_arr = df[[from_col, to_col]].max(axis=1).to_numpy(dtype=np.float64)
    scores = df['score'].to_numpy(dtype=np.float64)
    ko_rank = df['KO id'].rank(method='dense').to_numpy(dtype=np.float64)
    order = np.lexsort((ko_rank, -scores, ends_arr, starts))
    starts = starts[order]
    ends   = ends[order]
    orig_idx = orig_idx[order]

    grp_ids = _assign_groups_numba(starts, ends, float(frac_thresh))

    out = pd.Series(grp_ids, index=orig_idx).sort_index()
    return out.to_frame("grp_id")


def assign_overlap_groups(df_hits: pd.DataFrame) -> pd.DataFrame:
    """Add overlap_group labels, clustered per (target name, strand)."""
    if df_hits.empty:
        return df_hits.assign(overlap_group=pd.Series(dtype=str))

    df = df_hits.copy()

    out = []
    for (tgt, strand), sub in df.groupby(["target name", "strand"], sort=False):
        clustered = cluster_strand(sub)
        sub = sub.join(clustered, how="left")
        out.append(sub)

    result = pd.concat(out).sort_index()
    result["overlap_group"] = (
        result["target name"].astype(str) + "_" +
        result["grp_id"].astype(str) + "_" +
        result["strand"].astype(str)
    )
    return result


# Softmax + noise-floor winner selection 

def compute_softmax(df: pd.DataFrame) -> pd.DataFrame:
    """Per-position softmax within each overlap group (BLIMMP-aligned)."""
    df = df.copy().astype({'score': float})
    df['log_per_hit_weight'] = df['score'] * np.log(2.0)
    df['group_log_sum'] = (
        df.groupby('overlap_group')['log_per_hit_weight']
        .transform(lambda x: np.logaddexp.reduce(x.values))
    )
    with np.errstate(divide='ignore', invalid='ignore'):
        df['overlap_relative_position_confidence'] = (
            np.exp(df['log_per_hit_weight'] - df['group_log_sum']).fillna(0.0)
        )
    return df


def compute_winners(df: pd.DataFrame, e_threshold: float = 1e-4) -> pd.DataFrame:
    """
    Noise-floor hit confidence per overlap group, then pick one winner per group.
    Returns a reduced DataFrame with one row per overlap_group.
    """
    df = df.copy()

    noise_logw = -np.log(e_threshold)
    df['log_per_hit_weight'] = df['score'] * np.log(2.0)

    df['group_log_sum'] = (
        df.groupby('overlap_group')['log_per_hit_weight']
        .transform(lambda x: np.logaddexp.reduce(x.values))
    )

    df['total_log_weight'] = np.logaddexp(df['group_log_sum'], noise_logw)

    df['hit_conf'] = np.exp(df['log_per_hit_weight'] - df['total_log_weight'])

    # Debug NaNs
    nan_rows = df[df['hit_conf'].isna()]
    if not nan_rows.empty:
        print(f"Warning: {len(nan_rows)} rows with NaN hit_conf", file=sys.stderr)

    # Pick max-confidence row per overlap_group
    idx = (
        df.groupby('overlap_group')['hit_conf']
        .idxmax()
        .dropna()
        .astype(int)
    )

    winners = df.loc[idx].reset_index(drop=True)
    return winners
def filter_by_kofam_threshold(df: pd.DataFrame, ko_list_path: str) -> pd.DataFrame:
    """Drop winners whose score is below the KOfam threshold for their KO."""
    ko = pd.read_csv(ko_list_path, sep='\t', comment='#')
    ko.columns = ko.columns.str.strip()
    ko['knum']       = ko['knum'].astype(str).str.strip()
    ko['score_type'] = ko['score_type'].astype(str).str.strip().str.lower()
    ko['threshold']  = pd.to_numeric(ko['threshold'], errors='coerce')
    ko = ko[['knum', 'threshold', 'score_type']]

    original_columns = df.columns.tolist()  # save before merge

    merged = df.merge(ko, left_on='KO id', right_on='knum', how='left')

    no_thresh   = merged['threshold'].isna()
    full_pass   = (merged['score_type'] == 'full')   & (merged['score']   >= merged['threshold'])
    domain_pass = (merged['score_type'] == 'domain') & (merged['i_score'] >= merged['threshold'])
    keep = no_thresh | full_pass | domain_pass

    # Print proof of dropped rows
    dropped = merged[~keep][['KO id', 'score_type', 'score', 'i_score', 'threshold']].copy()
    dropped['score_used'] = dropped.apply(
        lambda r: r['score'] if r['score_type'] == 'full' else r['i_score'], axis=1
    )
    print("\nDropped rows (score < threshold):", file=sys.stderr)
    print(dropped[['KO id', 'score_type', 'score_used', 'threshold']].to_string(index=False), file=sys.stderr)
    print(f"\nThreshold filtering: {(~keep).sum()} rows dropped, {keep.sum()} kept", file=sys.stderr)

    # Restore exact original columns, no extras
    return merged.loc[keep, original_columns].reset_index(drop=True)
#  Main 

def main():
    ap = argparse.ArgumentParser(
        description="Process a HMMER .domtblout into one winner per overlap group "
                    "(BLIMMP-aligned pipeline)."
    )
    ap.add_argument("domtblout", help="Path to .domtblout file")
    ap.add_argument("--out", default="grouped_hits_winners.csv",
                    help="Output CSV (default: grouped_hits_winners.csv)")
    ap.add_argument("--e-threshold", type=float, default=1e-4,
                    help="Noise E-value threshold (default: 1e-4)")
    ap.add_argument("--ko-list", default=None,
                help="Path to KOfam ko_list file for threshold filtering")
    args = ap.parse_args()

    # Warm up Numba JIT
    _ = _assign_groups_numba(np.array([0., 1.]), np.array([1., 2.]), 0.6)
    _ = _hmm_union_len_per_group(
        np.array([0, 0], dtype=np.int64),
        np.array([0, 1], dtype=np.int64),
        np.array([1, 2], dtype=np.int64),
        1
    )

    # 1. Parse
    print("Processing domtblout file...", file=sys.stderr)
    hits = process_domtblout(args.domtblout)

    # 2. HMM union coverage
    print("Computing HMM union coverage...", file=sys.stderr)
    hits = compute_hmm_union_coverage(hits)
    hits = hits.sort_values(['target name', 'strand', 'ali from', 'ali to', 'score', 'KO id'],
                        ascending=[True, True, True, True, False, True]).reset_index(drop=True)

    # 3. Overlap groups
    print("Assigning overlap groups...", file=sys.stderr)
    grouped = assign_overlap_groups(hits)
    grouped = grouped.drop(columns=["grp_id"], errors="ignore")

    # 3b. Genome-wide KO dedup (keep top-scoring row per KO, matching BLIMMP)
    before_dedup = len(grouped)
    grouped = (grouped.sort_values('score', ascending=False)
               .drop_duplicates(subset='KO id', keep='first')
               .reset_index(drop=True))
    print(f"KO dedup: {before_dedup} → {len(grouped)} rows", file=sys.stderr)

    # 4. Softmax within groups
    print("Computing per-position softmax...", file=sys.stderr)
    grouped = compute_softmax(grouped)

    # 5. Winner per overlap group
    print("Selecting winners per overlap group...", file=sys.stderr)
    winners = compute_winners(grouped, e_threshold=args.e_threshold)

    print(f"Overlap groups: {grouped['overlap_group'].nunique()}", file=sys.stderr)
    print(f"Winners output: {len(winners)} rows", file=sys.stderr)

    # 6. Filter by KOfam threshold and save
    filtered = filter_by_kofam_threshold(winners, args.ko_list)

    filtered.to_csv(args.out, index=False)
    print(f"Written to {args.out}", file=sys.stderr)

    # Preview
    with pd.option_context('display.max_columns', None, 'display.width', 200):
        print("\nPreview (top 5 filtered winners):", file=sys.stderr)
        print(filtered.head(5), file=sys.stderr)

if __name__ == "__main__":
    main()
