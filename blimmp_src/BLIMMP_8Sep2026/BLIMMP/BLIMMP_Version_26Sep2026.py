#BLIMMP Source Code
# Development Version August 18
#Author: Neha S

import os
import glob
import re
import sys
import math
import ast
import json
import argparse
import pandas as pd
import numpy as np
import copy
from math import exp
from pathlib import Path
import operator as op
from dataclasses import dataclass
from numba import njit
from scipy.stats import beta
import warnings
import logging
from typing import Dict, Tuple, Set, Union, List, Any, Optional
    


## CONSTANTS
KINGDOM = {"bacillati", "fusobacteriati", "mycoplasmatota", "pseudomonadati", "thermotogati"}
PHYLUM  = {"bacillota", "acidobacteriota", "actinomycetota", "campylobacterota", "cyanobacteriota",
           "deinococcota", "fcb_group", "mycoplasmatota", "myxococcota", "pseudomonadota",
           "pvc_group", "spirochaetota", "thermodesulfobacteriota", "thermotogota"}
KO_RE = re.compile(r'^K\d{5}$')

## Configurations
@dataclass(frozen=True)
class Paths:
    counts_dir: Path
    onehop_dir: Path
    twohop_dir: Path
    module_neighbor_dir: Path
    module_eq_json: Path
    module_json_dir: Path
    kofam_ko_list_path: Path
    module_frequencies: Path
    module_reaction_dir: Path
    module_descriptions_path: Path
    ko_reaction_path: Path
    cooccurrence_lookup_dir: Optional[Path] = None


@dataclass(frozen=True)
class RunConfig:
    input_file: str
    fmt: str                   # 'tbl' or 'domtblout'
    sigma: float               # 0 to 1
    taxonomy: str              # bacteria / phylum / kingdom tag
    output_prefix: str
    verbose: bool = False
    logfile_path: Optional[str] = None
    no_substitutes: bool = False


@njit
def _hmm_union_len_per_group_py(gids, starts, ends, n_groups):
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


class HMMParsers:
    # Parse a HMMER domtblout file into a DataFrame of per-domain hits.
    # We auto-detect which column holds KO ids (query vs. target) and compute per-hit HMM coverage so downstream steps can decide how well each KO was found.
    @staticmethod
    def process_domtblout(path):
        print("Parsing HMM search results from: %s", path)
        # All columns defined by the domtblout format; we read the full set and then subset to only what BLIMMP actually uses below.
        cols = ['target name', 'target_accession', 'tlen', 'query_name',
            'query_accession', 'qlen', 'full_Evalue', 'full_score',
            'full_bias', 'n_domains', 'of_domains', 'c_Evalue',
            'i_Evalue', 'i_score', 'i_bias', 'hmm from', 'hmm to',
            'ali from', 'ali to', 'env from', 'env to', 'acc'
        ]

        usecols=[
            'target name','query_name',
            'hmm from','hmm to','tlen',
            'ali from','ali to','qlen',
            'full_score','full_Evalue',
            'i_score','i_Evalue'
        ]

        df = pd.read_csv(
            path,
            comment='#',
            header=None,
            names=cols,
            usecols=list(range(22)),
            sep= r"\s+",  
            engine='c',
            low_memory=False,
            memory_map=True
        )
        
        df = df[usecols].copy()

        # Auto-detect which column holds KEGG Ortholog ids (K#####).
        # HMMER output differs by tool: some put the KO in query_name (kofamscan),
        # others put it in target name. We check which column has more KO-shaped values.
        pat = KO_RE
        #Does the pattern match query or target
        q_matches = df['query_name'].astype(str).str.fullmatch(pat.pattern, na=False)
        t_matches = df['target name'].astype(str).str.fullmatch(pat.pattern, na=False)

        #Mean of number of matches
        fq = q_matches.mean()
        ft = t_matches.mean() 

        if (fq > 0 or ft > 0) and (fq >= ft):
            # query_name looks more like KO IDs
            df = df.rename(columns={'query_name': 'KO id', 'target name': 'target name'})
            df['hmm_len'] = df['qlen']
        elif ft > 0:
            # target name looks more like KO IDs
            df = df.rename(columns={'target name': 'KO id', 'query_name': 'target name'})
            df['hmm_len'] = df['tlen']
        else:
            raise ValueError("Can't detect KO id column. Check input format.")

        #df = df.rename(columns={'i_score': 'score', 'i_Evalue': 'E-value'})
        df = df.rename(columns={'full_score': 'score', 'full_Evalue': 'E-value'})

        # Compute alignment/hmm segment lengths and remove zero-length hits
        df['ali_span'] = (df['ali to'] - df['ali from']).abs()
        df['hmm_span'] = (df['hmm to'] - df['hmm from']).abs()
        df["strand"] = np.where(df["ali to"] >= df["ali from"], "+", "-")

        df = df[(df['ali_span'] > 0) & (df['hmm_span'] > 0)].copy()

        df['per_hit_hmm_coverage'] = (df['hmm_span'] + 1) / df['hmm_len']
        df.drop(columns=['ali_span', 'hmm_span'], inplace=True)

        # after building df, hmm_len, etc.
        df['group_id'] = df.groupby(['strand', 'target name', 'KO id'], sort=False).ngroup()

        n_groups = df['group_id'].max() + 1

        starts = df[['hmm from', 'hmm to']].min(axis=1).to_numpy(np.int64)
        ends   = df[['hmm from', 'hmm to']].max(axis=1).to_numpy(np.int64)
        gids   = df['group_id'].to_numpy(np.int64)

        order = np.lexsort((starts, gids))  # sort by gid, then start

        starts = starts[order]
        ends   = ends[order]
        gids   = gids[order]

        # union length per group (Numba)
        #covered_per_group = _hmm_union_len_per_group(gids, starts, ends, n_groups)

        # DEBUG: verify Numba vs Python
        covered_py = _hmm_union_len_per_group_py(gids, starts, ends, n_groups)


        # hmm_len per group
        hmm_len_per_group = (df.groupby('group_id', sort=False)['hmm_len'].first().to_numpy())

        coverage_per_group = covered_py / hmm_len_per_group

        # broadcast back
        gid_full = df['group_id'].to_numpy()
        df['hmm_covered_len']       = covered_py[gid_full]
        df['hmm_coverage_fraction'] = coverage_per_group[gid_full]
        return df


@njit
def _assign_groups_numba(starts, ends, frac_thresh):
    n = len(starts)
    g_st = np.empty(n, dtype=np.float64)
    g_en = np.empty(n, dtype=np.float64)
    gcount = 0
    grp_ids = np.empty(n, dtype=np.int32)

    for i in range(n):
        s = starts[i]; e = ends[i]
        assigned = False
        for g in range(gcount - 1, -1, -1):
            gs = g_st[g]; ge = g_en[g]
            if s > ge:
                break
            overlap = min(e, ge) - max(s, gs)
            if overlap <= 0.0:
                continue
            short_len = (e - s) if (e - s) < (ge - gs) else (ge - gs)
            if (overlap / short_len) >= frac_thresh:
                grp_ids[i] = g + 1
                if s < gs: g_st[g] = s
                if e > ge: g_en[g] = e
                assigned = True
                break
        if not assigned:
            g_st[gcount] = s
            g_en[gcount] = e
            gcount += 1
            grp_ids[i] = gcount
    return grp_ids


## Overlap grouping: HMM hits at the same genomic locus are grouped together so that competing annotations (e.g., two different KOs hitting the same ORF) can be adjudicated rather than treated as independent detections.
class Overlap:
    @staticmethod
    def cluster_strand(df, from_col="ali from", to_col="ali to", frac_thresh=0.6):
        # Compute strand-agnostic interval endpoints so forward and reverse hits are clustered correctly regardless of coordinate orientation.
        starts = df[[from_col, to_col]].min(axis=1).to_numpy(dtype=np.float64)
        ends   = df[[from_col, to_col]].max(axis=1).to_numpy(dtype=np.float64)
        orig_idx = df.index.to_numpy()

        # sort by start
        order = np.argsort(starts, kind="mergesort")
        starts = starts[order]; ends = ends[order]; orig_idx = orig_idx[order]

        # fast assignment
        grp_ids = _assign_groups_numba(starts, ends, float(frac_thresh))

        # return same shape/column as before
        out = pd.Series(grp_ids, index=orig_idx).sort_index()
        return out.to_frame("grp_id")

    @staticmethod
    def assign_overlap_groups(df_hits):
        print("Grouping overlapping HMM hits by genomic position...\n")
        df = df_hits.copy()

        out = []
        for (tgt, strand), sub in df.groupby(["target name","strand"], sort=False):
            clustered = Overlap.cluster_strand(sub)      # ← now uses the Numba path
            sub = sub.join(clustered, how="left")
            out.append(sub)

        result = pd.concat(out).sort_index()
        result["overlap_group"] = (
            result["target name"].astype(str)
            + "_" + result["grp_id"].astype(str)
            + "_" + result["strand"].astype(str)
        )
        return result


class File_Helpers:
    @staticmethod
    def load_kofamdb_file(kofampath):
        if not os.path.exists(kofampath):
            raise ImportError(f"KOfam DB file not found: {kofampath}")

        df = pd.read_csv(
        kofampath,
        sep=r"\s+",
        header=None,
        skiprows=1,
        usecols=[0, 1, 2],
        names=['KO id', 'kofam_score_threshold', 'score_type'])

        if df.shape[1] < 2:
            raise ValueError(f"KOfam file {kofampath} has <2 columns; can't parse thresholds.")

        df['KO id'] = df['KO id'].astype(str).str.strip().str.upper()
        df['kofam_score_threshold'] = pd.to_numeric(df['kofam_score_threshold'], errors='coerce')
        df = df.dropna(subset=['kofam_score_threshold'])

        if 'score_type' not in df.columns:
            df['score_type'] = 'full'
        else:
            df['score_type'] = (
                df['score_type']
                .fillna('full')
                .astype(str).str.strip().str.lower()
            )
            df.loc[~df['score_type'].isin(['full', 'domain']), 'score_type'] = 'full'

        
        return dict(zip(df['KO id'], zip(df['kofam_score_threshold'], df['score_type'])))
    
    @staticmethod
    def load_module_eq(module_eq_path):
        with open(module_eq_path, "r") as fh:
            module_equations = json.load(fh) 
        return module_equations
        
    @staticmethod
    def load_module_freq(module_fre_paths):
        module_freq = {}

        with open(module_fre_paths) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                module_id, freq = line.split("\t")
                module_freq[module_id] = float(freq)
        return module_freq
        
    @staticmethod
    def load_cooccurrence_lookup_targeted(path, needed_kos: set) -> dict:
        # Stream the co-occurrence JSON and extract only entries where KO_i is in needed_kos. Much lower memory than full load. Returns {KO_i: {KO_j: Nij}} for needed KOs only.

        if path is None:
            return {}
        path = Path(path)
        if not path.exists():
            logging.warning("Co-occurrence lookup not found: %s — pass 2 disabled.", path)
            return {}

        ko_i_pattern = re.compile(r'^\s{2}"(K\d{5})":\s*\{')
        ko_j_pattern = re.compile(r'^\s{4}"(K\d{5})":\s*([0-9.e+\-]+),?$')

        result = {}
        current_ko_i = None
        current_block = {}
        in_needed = False

        logging.info("Streaming co-occurrence lookup for %d KOs: %s", len(needed_kos), path.name)

        with open(path, "r") as fh:
            for line in fh:
                m_i = ko_i_pattern.match(line)
                if m_i:
                    # save previous block if needed
                    if in_needed and current_block:
                        result[current_ko_i] = current_block
                    current_ko_i = m_i.group(1)
                    current_block = {}
                    in_needed = current_ko_i in needed_kos
                    continue

                if in_needed and current_ko_i:
                    m_j = ko_j_pattern.match(line)
                    if m_j:
                        try:
                            current_block[m_j.group(1)] = float(m_j.group(2))
                        except ValueError:
                            pass

            # save last block
            if in_needed and current_block:
                result[current_ko_i] = current_block

        logging.info("Co-occurrence lookup: extracted %d KOs (of %d needed)", len(result), len(needed_kos))
        return result
        
    @staticmethod
    def read_ko_occurrence(kooccpath):
        df = pd.read_csv(kooccpath, sep="\t", skiprows=1, header=None, names=['KO id','count','occurences'])
        print("The prior taxonomical class chosen is from file: ", kooccpath,"\n")
        df['KO_freq'] = df['occurences'].astype(float)
        return df
    
    @staticmethod
    def lineage_paths(taxonomy: str, paths: Paths):
        val = (taxonomy or "").strip().lower()
        if val in PHYLUM: level, name = "phylum", val
        elif val in KINGDOM: level, name = "kingdom", val
        else: level, name = "domain", "bacteria"
        tag = "domain_level_priors" if level == "domain" else f"{name}_{level}_level_priors"
        want_counts = paths.counts_dir / f"ko_freq_ko_matrix_sampleids_{tag}.tsv"
        want_one    = paths.onehop_dir  / f"One_Hop_Refilled_{tag}.json"
        want_two    = paths.twohop_dir  / f"Two_Hop_Refilled_{tag}.json"
        want_all = paths.module_neighbor_dir / f"Module_AllHop_Refilled_{tag}.json"
        want_cooc   = (                                                        
        paths.cooccurrence_lookup_dir / f"CoOccurrence_AllVsAll_{tag}.json"
        if paths.cooccurrence_lookup_dir is not None else None)
        if not (want_counts.exists() and want_one.exists() and want_two.exists()):
            if tag != "domain_level_priors":
                print(f"[taxonomy] Using domain-level fallbacks for '{tag}'.", file=sys.stderr)
            tag = "domain_level_priors"
            want_counts = paths.counts_dir / f"ko_counts_ko_matrix_sampleids_{tag}.tsv"
            want_one    = paths.onehop_dir  / f"One_Hop_Refilled_{tag}.json"
            want_two    = paths.twohop_dir  / f"Two_Hop_Refilled_{tag}.json"
            want_all = paths.module_neighbor_dir / f"Module_AllHop_Refilled_{tag}.json"
            want_cooc   = (                                                     
            paths.cooccurrence_lookup_dir / f"CoOccurrence_AllVsAll_{tag}.json"
            if paths.cooccurrence_lookup_dir is not None else None)
        
        return want_counts, want_one, want_two, want_all, want_cooc, tag
    
    @staticmethod
    def modules_to_kos(module_json_dir):
        #Build comma-separated 'Mxxxxx,Myyyyy' mapping

        ko_to_modules: dict[str, list[str]] = {}
        pattern = os.path.join(module_json_dir, "module_*_nodes.json")

        for filepath in glob.glob(pattern):
            filename = os.path.basename(filepath)                          
            module_name = filename.replace("_nodes.json", "").replace("module_", "")  # like 'M00001'

            with open(filepath, "r") as f:
                module_nodes = json.load(f)                                 # list of node ids

            # Node names that start with K; keep KO part before first underscore
            module_kos = {n.split("_", 1)[0] for n in module_nodes if isinstance(n, str) and n.startswith("K")}

            for ko in module_kos:
                ko_to_modules.setdefault(ko, []).append(module_name)

        return {ko: ",".join(sorted(mods)) for ko, mods in ko_to_modules.items()}
    
    @staticmethod
    def ko_to_reactions_dict(ko_reaction_path):
        #Parse a flat KO-to-reaction file into dict
        ko_to_reactions: Dict[str, Set[str]] = {}
        with open(ko_reaction_path, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                
                ko  = parts[0].replace("ko:", "").strip().upper()
                rxn = parts[1].replace("rn:", "").strip().upper()
                
                if not KO_RE.match(ko):
                    continue
                
                ko_to_reactions.setdefault(ko, set()).add(rxn)
        logging.info("Loaded reactions for %d KOs from %s", len(ko_to_reactions), ko_reaction_path)
        return ko_to_reactions
    
    @staticmethod
    def load_module_reactions(module_reaction_dir):
        with open(module_reaction_dir, "r") as fh:
            module_reaction = json.load(fh) 
        return module_reaction
    
    @staticmethod
    def reactions_for_module_bestpath(module_map, module_id: str, best_path: str):
        if not module_id or not best_path:
            return ""

        mod = module_map.get(module_id)
        if not mod:
            return ""

        kos = re.findall(r"K\d{5}", str(best_path))

        rxns = set()
        for ko in kos:
            for r in mod.get(ko, []):
                rxns.add(r)

        return ",".join(sorted(rxns))
    
    @staticmethod
    def load_module_descriptions(module_descriptions_path) -> Dict[str, str]:
        #Parse kegg_bacteria_modules.json to {module_id: description_string}
        
        desc_path = str(module_descriptions_path)
        if not os.path.exists(desc_path):
            logging.warning("Module descriptions file not found: %s skipping.", desc_path)
            return {}

        with open(desc_path, "r", encoding="utf-8") as fh:
            raw = json.load(fh)

        module_desc: Dict[str, str] = {}
        for category, modules in raw.items():
            if not isinstance(modules, dict):
                continue
            for mod_id, info in modules.items():
                if isinstance(info, dict) and "Description" in info:
                    module_desc[mod_id] = info["Description"]

        logging.info("Loaded %d module descriptions from %s", len(module_desc), desc_path)
        return module_desc


  

class PositionScores:
    @staticmethod
    def compute_perposition_overlapgroup_softmax(df):
        #Calculate the softmax for each annotation in a group
        #Computation:
        # log_w_i = score_i * ln(2)
        # group_log_sum = logsumexp(log_w_i for each overlap_group)
        # X_i = exp(log_w_i - group_log_sum)
        df = df.copy().astype({'score': float})
        # Compute log(2**score) = score * ln(2)
        df['log_per_hit_weight'] = df['score'] * np.log(2.0)
        # log-sum-exp within each overlap_group
        df['group_log_sum'] = df.groupby('overlap_group')['log_per_hit_weight'].transform(lambda x: np.logaddexp.reduce(x.values))
        with np.errstate(divide='ignore', invalid='ignore'):
            df['overlap_relative_position_confidence'] = np.exp(df['log_per_hit_weight'] - df['group_log_sum']).fillna(0.0)
        return df
    

    @staticmethod
    def calculate_best_hit_with_noise(df,e_threshold=1e-4):
        df = df.copy()
        # noise term: noise_weight = 2**(-log2(e_threshold)) noise_logw = -ln(e_threshold)
        noise_logw = -np.log(e_threshold)
        # per-hit log-weight: ln(2**score) = score * ln(2)
        df['log_per_hit_weight'] = df['score'] * np.log(2)
        # group log-sum of per-hit weights
        df['group_log_sum'] = df.groupby('overlap_group')['log_per_hit_weight'].transform(lambda x: np.logaddexp.reduce(x.values))
        # total log-weight = log(group_sum + noise_weight)
        df['total_log_weight'] = np.logaddexp(df['group_log_sum'], noise_logw)
        # hit confidence = per_hit_weight / total_weight
        # log-space: exp(log_w − total_log_w)
        df['hit_conf'] = np.exp(df['log_per_hit_weight'] - df['total_log_weight'])
        # debug: print any nans
        nan_rows = df[df['hit_conf'].isna()][['score','log_per_hit_weight','group_log_sum','total_log_weight','hit_conf']]
        if not nan_rows.empty:
            print("Rows with NaN hit_conf:\n", nan_rows)
        # Pick the max-confidence row per overlap_group
        idx = (df.groupby('overlap_group')['hit_conf'].idxmax().dropna().astype(int))
        winners = df.loc[idx, ['overlap_group', 'KO id', 'score', 'hit_conf']].rename(
        columns={
            'KO id': 'overlapgroup_winner',
            'score': 'overlapgroup_winner_score',
            'hit_conf': 'overlapgroup_winner_hit_conf'})
        return winners.loc[idx].reset_index(drop=True)
    
    @staticmethod
    def winner_info_and_flags(df, kofampath):
        print("Calculating hit confidence for each KO ids...\n")
        df_soft = PositionScores.compute_perposition_overlapgroup_softmax(df)
        winners = PositionScores.calculate_best_hit_with_noise(df)
        keep_cols = ['overlap_group', 'overlapgroup_winner', 'overlapgroup_winner_score', 'overlapgroup_winner_hit_conf']
        df_new = df_soft.merge(winners[keep_cols], on='overlap_group', how='left', validate='many_to_one')

        #Merge the kofam score threshold here
        kofam_map = None
        try:
            # load_kofamdb_file
            kofam_map = File_Helpers.load_kofamdb_file(kofampath)
            if len(kofam_map) > 0:
                # map thresholds to rows
                df_new['kofam_score_threshold'] = df_new['KO id'].map(lambda ko: kofam_map.get(ko, (np.nan, None))[0])
                df_new['kofam_score_type'] = df_new['KO id'].map(lambda ko: kofam_map.get(ko, (np.nan, None))[1])

                # conditions to pass
                compare_score = pd.Series(np.where(df_new['kofam_score_type'].eq('domain'),df_new['i_score'],df_new['score']),index=df_new.index)
                conditions = (df_new['kofam_score_threshold'].notna() & compare_score.notna() & (compare_score >= df_new['kofam_score_threshold']))

                # Outcompeted flag
                df_new['is_outcompeted'] = (df_new['KO id'] != df_new['overlapgroup_winner'])

                #Below threshold flag 
                # ? = Outcompeted and above kofam threshold
                # ! = Not outcompeted (winner) and below kofam threshold
                # below threshold, has score, not a winner
                # below threshold no score

                # hit_conf = 1.0 when above threshold, else relative confidence
                df_new['hit_conf'] = np.where(conditions,1.0,df_new['overlap_relative_position_confidence'])
                # flag is_dubious = True if passes threshold but NOT the overlapgroup winner
                df_new['flag_is_dubious'] = conditions & (df_new['KO id'] != df_new['overlapgroup_winner'])

                # flag is_below_kofam_threshold = True when threshold missing OR score < threshold
                df_new['flag_is_below_kofam_threshold'] = ~conditions
        except ImportError:
            kofam_map = None

        return df_new


class NeighborCalculations:
    @staticmethod
    def make_neighbor_dictionary(NEIGHBOR_TXT, df=None):
        with open(NEIGHBOR_TXT) as f:
            neighbor_data = json.load(f)
            
        #Nj_count = 0
        #Nij_count=0

        adj_raw = {}
        ko_counts = {}
        
        for ko, nbrs in neighbor_data.items():
            ko_counts[ko] = float(nbrs.get("_count", 0.0))
            
            
        for ko, nbrs in neighbor_data.items():
            out = {}
            for nb, val in nbrs.items():
                if nb == "_count":
                    continue
                
                Nij = float(val)
                Nj = ko_counts.get(nb)
                
                if Nj is None:
                    raise ValueError(f"Missing _count for neighbor KO '{nb}' " f"(referenced from '{ko}')")
                    #Nj_count += 1

                #if Nij > Nj:
                #    raise ValueError(f"Inconsistent counts: Nij > Nj "f"(i='{ko}', j='{nb}', Nij={Nij}, Nj={Nj})")
                #    Nij_count += 1
                
                out[nb] = Nij
            if out:
                adj_raw[ko] = out

        #print(len(neighbor_data))
        #print(Nj_count,Nij_count)
        #print(len(ko_counts))
        return adj_raw, ko_counts
    

class CalculateKOProbabilities:
    @staticmethod
    def sigma_completeness_alteration(df, sigma_val):
        sigma_val_update = 1 - ((np.exp(3 * sigma_val) - 1) / np.exp(3))
        df['sigma'] = sigma_val_update
        return df


    @staticmethod
    def calculate_dk_per_ko(df, ko_occ, verbose: bool = False, M=0.7):
        logging.info("Calculating per-KO probabilities...")

        df = pd.merge(df, ko_occ, on="KO id", how='left')
        df['KO_freq'] = df['KO_freq'].fillna(0.0)
        
        # Check how many hit_conf values are missing
        #missing_conf = df['hit_conf'].isna().sum()
        #print(f"Missing hit_conf values: {missing_conf} out of {total_rows} ({100*missing_conf/total_rows:.2f}%)")
        # Calculate Dk
        seed = df['hit_conf'].fillna(0.0)
        df['Dk'] = seed + (1 - seed) * 0.7 * df['KO_freq'] * df['KO_freq']
        
        if verbose:
            for _, row in df.iterrows():
                term = (1 - row['hit_conf']) * M * (row['KO_freq'] ** 2)
                logging.debug(
                    f"KO id: {row['KO id']}, "
                    f"O_i (hit_conf): {row['hit_conf']:.4f}, "
                    f"F_i (KO_freq): {row['KO_freq']:.4f}, "
                    f"(1-O_i)*M*F_i²: {term:.4f} = P_i: {row['Dk']:.4f}"
                )

        
        return df
    
    @staticmethod
    def calculate_reliable_conditional_prob(i, j, neighbor_map, ko_counts, lambda_param=50, verbose=False):
        #Calculate R(i,j) = N_{i and j} / (N_j + λ)
        #min permissible influencer
        min_frac = 0.25

        # Get co-occurrence count: N_{i and j}
        n_i_and_j = neighbor_map.get(i, {}).get(j, 0)
        
        # Get count of j: N_j
        n_j = ko_counts.get(j, 0)
        
        # Calculate R(i,j) with lambda regularization
        if n_j <= 0:
            if verbose:
                logging.debug(
                    "KO %s <- buddy %s excluded: Nj=0 (no support)",
                    i, j
                )
            return 0.0

        frac = n_i_and_j /n_j

        if frac < min_frac:
            if verbose:
                logging.debug(
                    "KO %s <- buddy %s excluded: Nij/Nj = %.3f < %.2f "
                    "(Nij=%d Nj=%d)",
                    i, j, frac, min_frac, n_i_and_j, n_j
                )
            return 0.0
        
        r_ij = n_i_and_j / (n_j + lambda_param)
        
        return r_ij



    @staticmethod
    def dk_neighbor_update(df, neighbor_map, ko_counts, alpha=0.6, lambda_param=50, return_used=False, verbose: bool = False):
        logging.info("Updating per-KO probabilities based on the influence neighborhood...\n")
        s = df[['KO id','Dk','hit_conf','count']].copy()
         
        s['KO id']    = s['KO id'].astype(str).str.strip().str.upper()
        s['Dk']       = pd.to_numeric(s['Dk'], errors='coerce').fillna(0.0)
        s['hit_conf'] = pd.to_numeric(s['hit_conf'], errors='coerce').fillna(0.0)
        s['count']    = pd.to_numeric(s['count'], errors='coerce').fillna(0.0)

        dk_dict                = dict(zip(s['KO id'], s['Dk']))
        hit_conf_map_current   = dict(zip(s['KO id'], s['hit_conf']))

        new_dk = {}
        used_neighbors = {}

        buddy_stats_map = {}   #for viz
        

        for i, p_i in dk_dict.items():
            module_families = neighbor_map.get(i, {})
            #Default
            buddy_stats_map[i] = {
                "alpha": float(alpha),
                "lambda": float(lambda_param),
                "S": 0.0,
                "X": None,
                "shift": 0.0,
                "buddy_count_used": 0,
                "buddies": []
            }

            if not module_families:
                new_dk[i] = p_i
                used_neighbors[i] = []
                if verbose:
                    logging.debug(
                        "KO %s: no neighbors, C_i = P_i = %.4f",
                        i, p_i
                    )
                continue


            # Calculate R(i,j) for each family j
            rij_map = {}
            for j in module_families.keys():
                if j == i:
                    continue
                r_ij = CalculateKOProbabilities.calculate_reliable_conditional_prob(
                    i, j, neighbor_map, ko_counts, lambda_param, verbose
                )
                if r_ij > 0.0 and math.isfinite(r_ij):
                    #squared for dampening
                    rij_map[j] = r_ij #*r_ij

            used_neighbors[i] = sorted(rij_map.keys())
            
            if not rij_map:
                new_dk[i] = p_i
                if verbose:
                    logging.debug("KO %s: neighbors found but no reliable R(i,j); C_i = P_i = %.4f",i, p_i)
                continue

            # Sum of reliable conditional probabilities (weights)
            S = sum(rij_map.values())

            if S <= 0.0:
                new_dk[i] = p_i
                if verbose:
                    logging.debug("KO %s: S <= 0 (S=%.4f); C_i = P_i = %.4f",i, S, p_i)
                continue
                
            # Spring calculation with R(i,j) as weights
            X = alpha ** (1.0 / S)
            a_i = 1.0 - p_i

            # Update using buddy influence
            shift = 0.0

            r_pj_values = []
            buddies_list = []
            for j, r_ij in rij_map.items():
                pj = dk_dict.get(j, 0.0)
                weight = r_ij / S
                contrib = a_i * weight * X * pj
                shift += contrib
                r_pj_values.append(r_ij * pj)
                
                
                Nij = neighbor_map[i].get(j, 0.0)
                Nj  = ko_counts.get(j, 0.0)
                
                if Nij < 0 or Nj < 0:
                    raise ValueError(f"Negative counts: i={i} j={j} Nij={Nij} Nj={Nj}")

                if Nij > Nj:
                    raise ValueError(f"Inconsistent counts (Nij > Nj): i={i} j={j} Nij={Nij} Nj={Nj}. " f"Check neighbor_map counts vs ko_counts definition.")
                    
                if r_ij > 1.0 + 1e-12:
                    raise ValueError(f"R(i,j) > 1: i={i} j={j} R={r_ij} Nij={Nij} Nj={Nj} lambda={lambda_param}")

                if verbose:
                    Nj_l = Nj + lambda_param
                    logging.debug("KO %s <- buddy %s | Nij=%.4f Nj=%.4f Nj+λ=%.4f R(i,j)=%.4f "
                        "weight=%.4f pj=%.4f contrib=%.6f",i, j, Nij, Nj, Nj_l, r_ij, weight, pj, contrib)
                
                buddies_list.append({
                    "ko": j,
                    "Nij": float(Nij),
                    "Nj": float(Nj),
                    "Rij": float(r_ij),
                    "weight": float(weight),
                    "pj": float(pj),          
                    "contrib": float(a_i * weight * X * pj)  
                })    

            new_val = min(p_i + shift, 1.0)
            new_dk[i] = new_val

            # write buddy stats for THIS i (inside loop)
            buddy_stats_map[i] = {
                "alpha": float(alpha),
                "lambda": float(lambda_param),
                "S": float(S),
                "X": float(X),
                "shift": float(shift),
                "buddy_count_used": int(len(buddies_list)),
                "buddies": sorted(buddies_list, key=lambda x: x["weight"], reverse=True),
            }

            if verbose:
                strong_buddies = sum(
                    1 for j, r_ij in rij_map.items()
                    if r_ij > 0.7 and dk_dict.get(j, 0.0) > 0.8
                )
                max_signal = max(r_pj_values) if r_pj_values else 0.0

                logging.debug(
                    "KO %s summary | P_i=%.4f C_i=%.4f shift=%.4f buddies=%d "
                    "| strong_buddies=%d | max(R*P_j)=%.4f",
                    i, p_i, new_val, shift, len(rij_map),
                    strong_buddies, max_signal
                )


        # Create output DataFrame with C_i
        df = df.copy()
        df['Dk_Neighbor'] = df['KO id'].map(new_dk).fillna(df['Dk'])
        for ko in dk_dict.keys():
            buddy_stats_map.setdefault(ko, None)
        
        df['buddy_stats'] = df['KO id'].map(buddy_stats_map)
        if return_used:
            return df, used_neighbors
        return df
        

def logm(level, mod_id, msg, *args):
    logging.log(level, "module=%s " + msg, mod_id, *args)

def logmk(level, mod_id, ko, msg, *args):
    logging.log(level, "module=%s ko=%s " + msg, mod_id, ko, *args)



class CalculateModuleProbabilities:
    _ALLOWED_FUNCS  = {"max": max, "min": min}
    _ALLOWED_BINOPS = {ast.Add: op.add, ast.Sub: op.sub, ast.Mult: op.mul, ast.Div: op.truediv}
    _ALLOWED_UNARY  = {ast.UAdd: op.pos, ast.USub: op.neg}

    KO_TOKEN_EXTRACT = re.compile(r'K(\d{5})(?:\w+)?')
    KO_TOKEN_STRICT  = re.compile(r'^K\d{5}$')

    # Matches "-K12345" only when "-" is UNARY:
    # start of string OR preceded by "(" or "," or "+" or "*" or "/" (possibly with spaces)
    _OPT_KO = re.compile(r'(?:(?<=^)|(?<=[(,+*/]))\s*-\s*(K\d{5})\b')

    @staticmethod
    def _ignore_optional_kos(eq: str) -> str:
        """
        Double fix. Replace unary optional KO tokens like '-K12345' with a neutral value so the
        expression remains valid. Does NOT touch '1 - K12345' (binary subtraction).
        """
        CMP = CalculateModuleProbabilities  
        s = eq

        def repl(m: re.Match) -> str:
            # Decide neutral element by nearby operator context.
            # Look backwards from match start to find the last non-space char.
            start = m.start()
            j = start - 1
            while j >= 0 and s[j].isspace():
                j -= 1
            prev = s[j] if j >= 0 else ''

            # If it appears right after a '*', neutral factor is 1 (ignore in product).
            # Otherwise, use 0 (ignore in sums / function args).
            return "1" if prev == "*" else "0"

        return CMP._OPT_KO.sub(repl, s)

    @staticmethod
    def _eval_ast(node, env):
            CMP = CalculateModuleProbabilities  

            if isinstance(node, ast.Expression):
                return CMP._eval_ast(node.body, env)

            if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
                return float(node.value)

            if isinstance(node, ast.Name):
                name = node.id
                if not CMP.KO_TOKEN_STRICT.fullmatch(name):
                    raise ValueError(f"Unknown variable '{name}'")
                return float(env.get(name, 0.0))

            if isinstance(node, ast.UnaryOp) and type(node.op) in CMP._ALLOWED_UNARY:
                return CMP._ALLOWED_UNARY[type(node.op)](CMP._eval_ast(node.operand, env))

            if isinstance(node, ast.BinOp) and type(node.op) in CMP._ALLOWED_BINOPS:
                left = CMP._eval_ast(node.left, env)
                right = CMP._eval_ast(node.right, env)
                return CMP._ALLOWED_BINOPS[type(node.op)](left, right)

            if isinstance(node, ast.Call):
                if isinstance(node.func, ast.Name) and node.func.id in CMP._ALLOWED_FUNCS:
                    func = CMP._ALLOWED_FUNCS[node.func.id]
                    args = [CMP._eval_ast(a, env) for a in node.args]
                    return float(func(*args))
                raise ValueError("Only max(...) and min(...) calls are allowed")

            if isinstance(node, ast.Tuple):
                return tuple(CMP._eval_ast(elt, env) for elt in node.elts)

            raise ValueError(f"Unsupported expression element: {ast.dump(node)}")


    @staticmethod
    def _normalize_equation(eq):
        # Replace any Kxxxxx_suffix with plain Kxxxxx (e.g., K00844_xyz to K00844)
        return CalculateModuleProbabilities.KO_TOKEN_EXTRACT.sub(
            lambda m: f"K{m.group(1)}", eq
        )

    @staticmethod
    def eval_equation(eq, dk_map, debug=False):
        CMP = CalculateModuleProbabilities
        eq_norm = CMP._normalize_equation(eq)
        # ignore unary optional -Kxxxxx tokens
        eq_clean = CMP._ignore_optional_kos(eq_norm)
        # Build env only for tokens present
        tokens = {f"K{m}" for m in CMP.KO_TOKEN_EXTRACT.findall(eq_clean)}
        env = {k: float(dk_map.get(k, 0.0)) for k in tokens}
        if debug:
            print("\neval_equation: Evaluating:", eq)
            print("Tokens:", tokens)
            print("Env:", env)


        val = float(CMP._eval_ast(ast.parse(eq_clean, mode="eval"), env))
        val_clamped = max(0.0, min(1.0, val))
        if debug: print(f"→ Result {val} (clamped {val_clamped})\n")
        return val_clamped

    @staticmethod
    def build_dk_maps_from_df(df: pd.DataFrame,before_col: str = "Dk",after_col: str = "Dk_Neighbor") -> tuple[dict, dict]:
        # Normalize KO ids
        s = df[["KO id", before_col]].copy()
        s["KO id"] = s["KO id"].astype(str).str.strip().str.upper()
        Dk = dict(zip(s["KO id"], pd.to_numeric(s[before_col], errors="coerce").fillna(0.0)))

        if after_col in df.columns:
            t = df[["KO id", after_col]].copy()
            t["KO id"] = t["KO id"].astype(str).str.strip().str.upper()
            Dk_Neighbor = dict(zip(t["KO id"], pd.to_numeric(t[after_col], errors="coerce").fillna(0.0)))
        else:
            Dk_Neighbor = Dk  # fallback

        return Dk, Dk_Neighbor

    @staticmethod
    def evaluate_step_probabilities(module_dict: dict,
                                   df: pd.DataFrame,
                                   before_col: str = "Dk",
                                   after_col: str = "Dk_Neighbor",  verbose: bool = False) -> pd.DataFrame:
        CMP = CalculateModuleProbabilities
        Dk, Dk_Neighbor = CMP.build_dk_maps_from_df(df, before_col, after_col)

        log_lines = []
        if verbose:
            log_lines.append("=" * 80)
            log_lines.append("DEBUG: STEP-LEVEL PROBABILITIES")
            log_lines.append("=" * 80)
            log_lines.append(f"Total modules in module_dict: {len(module_dict)}")

        rows_steps = []
        for mod_id, entry in module_dict.items():

            mod_eq = entry.get("module_equation", "")
            steps  = entry.get("steps", [])

            if verbose:
                logm(logging.DEBUG, mod_id, "module_equation=%s", mod_eq)
                logm(logging.DEBUG, mod_id, "n_steps=%d", len(steps))

            # Step-level
            for s in steps:
                idx = int(s["step"])
                eqn = s["equation"]

                p_b = CMP.eval_equation(eqn, Dk)
                p_a = CMP.eval_equation(eqn, Dk_Neighbor)

                # Extract KO IDs from equation (simple regex)
                kos_in_eq = sorted(set(re.findall(r"K\d{5}", eqn)))


                for ko in kos_in_eq:
                    v_b = Dk.get(ko, 0.0)
                    v_a = Dk_Neighbor.get(ko, 0.0)
                    log_lines.append(f"  {ko:<10} {v_b:<14.4f} {v_a:<20.4f}")    

                # Step-level probabilities
                rows_steps.append({
                    "module": mod_id,
                    "multiline": False,
                    "step": idx,
                    "equation": eqn,
                    "p_before": p_b,
                    "p_after":  p_a,
                })

                if verbose:
                    logm(logging.DEBUG, mod_id, "step=%d equation=%s", idx, eqn)
                    logm(logging.DEBUG, mod_id, "step=%d kos=%s", idx, ",".join(kos_in_eq))
                    for ko in kos_in_eq:
                        logmk(logging.DEBUG, mod_id, ko, "step=%d Dk=%.4f Dk_Neighbor=%.4f",
                            idx, Dk.get(ko, 0.0), Dk_Neighbor.get(ko, 0.0))
                    logm(logging.DEBUG, mod_id, "step=%d p_before=%.6f p_after=%.6f", idx, p_b, p_a)


        
        steps_df   = pd.DataFrame(rows_steps).sort_values(["module","step"]).reset_index(drop=True)
        if verbose:
            all_modules  = set(module_dict.keys())
            step_modules = set(steps_df["module"].unique())
            missing      = sorted(all_modules - step_modules)

            logging.debug("=" * 80)
            logging.debug("MODULE COVERAGE CHECK")
            logging.debug("Total modules in module_dict: %d", len(all_modules))
            logging.debug("Total modules in steps_df:     %d", len(step_modules))
            logging.debug("Missing modules:              %d", len(missing))
            if missing:
                logging.debug("First few missing: %s", ", ".join(missing[:20]))
            logging.debug("=" * 80)
        return steps_df
    
    @staticmethod
    def calculate_confidence(E: float, n: int, freq: float, thresh: float, prior_str: float) -> float:
        
        alpha_prior = freq * prior_str
        beta_prior = (1 - freq) * prior_str
        alpha_post = alpha_prior + E
        beta_post = beta_prior + (n - E)
        confidence = 1 - beta.cdf(thresh, alpha_post, beta_post)

        return confidence
    
    @staticmethod
    def calculate_hdi(alpha_post: float, beta_post: float, ci_level: float = 0.95, resolution: int = 1000) -> Tuple[float, float]:
        # Degenerate cases: if posterior is essentially a point mass
        pdf_sum = beta.pdf(np.array([0.5]), alpha_post, beta_post).sum()
        if not np.isfinite(pdf_sum) or pdf_sum == 0:
            mean = alpha_post / (alpha_post + beta_post) if (alpha_post + beta_post) > 0 else 0.0
            return (mean, mean)

        x = np.linspace(0, 1, resolution)
        pdf_vals = beta.pdf(x, alpha_post, beta_post)
        total = pdf_vals.sum()
        if not np.isfinite(total) or total == 0:
            mean = alpha_post / (alpha_post + beta_post) if (alpha_post + beta_post) > 0 else 0.0
            return (mean, mean)

        sorted_idx = np.argsort(-pdf_vals)
        cumulative = np.cumsum(pdf_vals[sorted_idx]) / total
        in_hdi = sorted_idx[cumulative <= ci_level]

        if len(in_hdi) == 0:
            mean = alpha_post / (alpha_post + beta_post) if (alpha_post + beta_post) > 0 else 0.0
            return (mean, mean)

        return float(x[in_hdi].min()), float(x[in_hdi].max())
    
    @staticmethod
    def calculate_credible_stats(alpha_post: float, beta_post: float, ci_level: float = 0.95) -> dict:
        tail = (1 - ci_level) / 2
        ab = alpha_post + beta_post
        return {
            "posterior_mean": alpha_post / ab,
            "posterior_variance": (alpha_post * beta_post) / (ab ** 2 * (ab + 1)),
            "ci_low": beta.ppf(tail, alpha_post, beta_post),
            "ci_high": beta.ppf(1 - tail, alpha_post, beta_post),}
    
    @staticmethod
    def annotate_module_ko_warnings(modules_df, ko_df, severe_ratio: float = 0.5):
        ko_flags = ko_df.groupby("KO id").agg({"is_outcompeted": "any","flag_is_below_kofam_threshold": "any","score": "max","kofam_score_threshold": "first",}).to_dict(orient="index")
        warning_texts = []
        flagged_ko_counts = []
        severe_ko_counts = []
        total_ko_counts = []

        for _, row in modules_df.iterrows():
            best_kos_str = row.get("module_best_path_kos", "")
            if not best_kos_str or pd.isna(best_kos_str):
                warning_texts.append("")
                flagged_ko_counts.append(0)
                severe_ko_counts.append(0)
                total_ko_counts.append(0)
                continue

            best_kos = [k.strip() for k in str(best_kos_str).split(",") if k.strip()]
            total_ko_counts.append(len(best_kos))
            warnings = []
            n_flagged = 0
            n_severe = 0


            for ko in best_kos:
                info = ko_flags.get(ko)
                if info is None:
                    continue

                ko_score = info["score"]
                ko_thresh = info["kofam_score_threshold"]
                parts = []

                if info["is_outcompeted"]:
                    parts.append("outcompeted at locus")
                    n_flagged += 1

            
                if info["flag_is_below_kofam_threshold"] and pd.notna(ko_thresh) and pd.notna(ko_score):
                    ratio = ko_score / ko_thresh if ko_thresh > 0 else 0.0
                    gap = ko_thresh - ko_score
                    parts.append(f"below KOfam threshold by {gap:.1f} "
                                f"(score={ko_score:.1f}, threshold={ko_thresh:.1f}, "
                                f"ratio={ratio:.2f})")
                    n_flagged += 1
                    if ratio < severe_ratio:
                        n_severe += 1

                
                elif info["flag_is_below_kofam_threshold"]:
                    parts.append("below KOfam threshold (no score available)")
                    n_flagged += 1
                    n_severe += 1

                if parts:
                    warnings.append(f"{ko}: {'; '.join(parts)}")

            flagged_ko_counts.append(n_flagged)
            severe_ko_counts.append(n_severe)
            warning_texts.append(" | ".join(warnings) if warnings else "")


        modules_df = modules_df.copy()
        modules_df["best_path_warnings"] = warning_texts
        modules_df["best_path_flagged_kos"] = flagged_ko_counts
        modules_df["best_path_severe_kos"] = severe_ko_counts
        modules_df["best_path_total_kos"] = total_ko_counts
        return modules_df

    
    @staticmethod
    def classify_module(confidence: float, ci_low: float, ci_high: float, hdi_low: float, hdi_high: float, threshold: float, posterior_mean: float, n_steps: int, effect_size: float, best_path_warnings: str = "", best_path_flagged_kos: int = 0, best_path_severe_kos: int = 0, best_path_total_kos=0):

        ci_width = ci_high - ci_low

        # base tier
        if confidence >= 0.95 and ci_low >= threshold:
            tier = "present"
        elif confidence >= 0.9:
            if n_steps <= 2 and ci_width > 0.50:
                tier = "prior_driven"
            else:
                tier = "likely_present"
        elif confidence >= 0.3:
            tier = "uncertain"
        elif confidence >= 0.05:
            tier = "likely_absent"
        else:
            tier = "absent"

        # Downgrading when KOs are severely below threshold
        downgraded = False
        if best_path_severe_kos > 0 and best_path_total_kos > 0:
            severe_fraction = best_path_severe_kos / best_path_total_kos
            if severe_fraction >= (1/3):
                # Significant fraction of best path is dubious then downgrade
                if tier == "present":
                    tier = "likely_present"
                    downgraded = True
                elif tier in ("likely_present", "prior_driven"):
                    tier = "dubious"
                    downgraded = True

        if tier == "present":
            status = (f"Present: confidence={confidence:.3f} with 95% credible interval "
                    f"[{ci_low:.2f}, {ci_high:.2f}] and HDI [{hdi_low:.2f}, {hdi_high:.2f}] "
                    f"entirely above threshold={threshold:.2f}. "
                    f"Posterior mean={posterior_mean:.3f}, effect size={effect_size:.2f}.")

        elif tier == "likely_present":
            status = (f"Likely present: confidence={confidence:.3f} and posterior mean={posterior_mean:.3f} "
                    f"above threshold={threshold:.2f}, but credible interval "
                    f"[{ci_low:.2f}, {ci_high:.2f}] extends below threshold. "
                    f"HDI [{hdi_low:.2f}, {hdi_high:.2f}], effect size={effect_size:.2f}.")

        elif tier == "prior_driven":
            status = (f"Likely present (prior-driven): confidence={confidence:.3f} but based on "
                    f"only {n_steps} step(s) with wide credible interval "
                    f"[{ci_low:.2f}, {ci_high:.2f}] (width={ci_width:.2f}). "
                    f"HDI [{hdi_low:.2f}, {hdi_high:.2f}]. "
                    f"High confidence may be due to prior (frequency and/or genomic evidence) instead of influence.")

        elif tier == "dubious":
            status = (f"Dubious: confidence={confidence:.3f} but {best_path_severe_kos} KO(s) on "
                    f"best path score below 50% of their KOfam threshold. "
                    f"Credible interval [{ci_low:.2f}, {ci_high:.2f}], "
                    f"HDI [{hdi_low:.2f}, {hdi_high:.2f}]. "
                    f"Module detection may be based on weak or spurious HMM hits.")

        elif tier == "uncertain":
            status = (f"Uncertain: confidence={confidence:.3f} with credible interval "
                    f"[{ci_low:.2f}, {ci_high:.2f}] straddling threshold={threshold:.2f}. "
                    f"Posterior mean={posterior_mean:.3f}, HDI [{hdi_low:.2f}, {hdi_high:.2f}]. "
                    f"Module may be present but evidence is inconclusive.")

        elif tier == "likely_absent":
            status = (f"Likely absent: confidence={confidence:.3f}, posterior mean={posterior_mean:.3f} "
                    f"below threshold={threshold:.2f}. Credible interval "
                    f"[{ci_low:.2f}, {ci_high:.2f}] mostly below threshold. "
                    f"Some genomic evidence exists but is insufficient.")

        else:
            status = (f"Absent: confidence={confidence:.3f}, posterior mean={posterior_mean:.3f}. "
                    f"Credible interval [{ci_low:.2f}, {ci_high:.2f}] and "
                    f"HDI [{hdi_low:.2f}, {hdi_high:.2f}] well below threshold={threshold:.2f}. "
                    f"No meaningful genomic evidence.")

        if downgraded:
            logging.debug("Module downgraded — %s", best_path_warnings)
            status += f" DOWNGRADED due to: {best_path_warnings}"
        elif best_path_flagged_kos > 0:
            status += f" Note: {best_path_flagged_kos} KO(s) flagged but within acceptable range."

        return status




    @staticmethod
    def calculate_module_confidence(steps_df: pd.DataFrame, module_dict: dict, genome_completeness: float, module_frequencies: dict = None, default_frequency: float = 0.5,  prior_strength: float = 1.0, default_beta_thresh: float = 0.65,  verbose: bool = False) -> pd.DataFrame:

        CMP = CalculateModuleProbabilities

        if module_frequencies is None:
            module_frequencies = {}

        # Validate genome completeness (warn but proceed)
        if genome_completeness < 0.4:
            warnings.warn(
                f"Genome completeness ({genome_completeness:.2f}) is very low (< 0.4). "
                "This may lead to spurious results. Proceed with caution.", UserWarning
            )
        
        if genome_completeness < 0 or genome_completeness > 1.0:
            raise ValueError(
                f"Genome completeness must be between 0 and 1.0, got {genome_completeness}"
            )
    
        beta_threshold = default_beta_thresh * genome_completeness

        rows_modules = []

        for mod_id in steps_df['module'].unique():
            # Get steps for this module
            module_steps = steps_df[steps_df['module'] == mod_id]
            
            # Calculate E (sum of step probabilities) and n_steps
            E_before = module_steps['p_before'].sum()
            E_after = module_steps['p_after'].sum()
            n_steps = len(module_steps)

            # Get module frequency (enforce minimum of 0.001 for zero values)
            freq = module_frequencies.get(mod_id, default_frequency)
            if freq == 0:
                freq = 0.001
            if freq == 1.0:
                freq = 0.999
        
            # Compute priors and posteriors explicitly (for debug visibility)
            prior_strength = n_steps*0.1
            alpha_prior = freq * prior_strength
            beta_prior = (1.0 - freq) * prior_strength

            alpha_post_before = alpha_prior + E_before
            beta_post_before  = beta_prior + (n_steps - E_before)

            alpha_post_after  = alpha_prior + E_after
            beta_post_after   = beta_prior + (n_steps - E_after)

            # Calculate confidence using Bayesian approach
            conf_before = CMP.calculate_confidence(E_before, n_steps, freq, beta_threshold, prior_strength)
            conf_after = CMP.calculate_confidence(E_after, n_steps, freq, beta_threshold, prior_strength)
                    
            credible_stats = CMP.calculate_credible_stats(alpha_post_after, beta_post_after)
            hdi_low, hdi_high = CMP.calculate_hdi(alpha_post_after, beta_post_after)
            effect_size = (credible_stats["posterior_mean"] - beta_threshold) / math.sqrt(credible_stats["posterior_variance"]) if credible_stats["posterior_variance"] > 0 else 0.0


            mod_eq = module_dict.get(mod_id, {}).get("module_equation", "")

            #module_status = CMP.classify_module(conf_after, credible_stats["ci_low"], credible_stats["ci_high"],hdi_low, hdi_high, beta_threshold, credible_stats["posterior_mean"], n_steps, effect_size)
            
            rows_modules.append({
                "module": mod_id,
                "module_equation": mod_eq,
                "n_steps": n_steps,
                "E_before": E_before,
                "E_after": E_after,
                "module_frequency": freq,
                "module_probability_before": conf_before,
                "module_probability_after": conf_after,
                "posterior_mean": credible_stats["posterior_mean"],
                "posterior_variance": credible_stats["posterior_variance"],
                "beta_threshold": beta_threshold,
                "ci_low": credible_stats["ci_low"],
                "ci_high": credible_stats["ci_high"],
                "hdi_low": hdi_low,
                "hdi_high": hdi_high,
                "effect_size": effect_size,
                #"Module_Status":module_status
            })

            if verbose:
                logm(logging.DEBUG, mod_id, "n_steps=%d", n_steps)
                logm(logging.DEBUG, mod_id, "E_before=%.4f E_after=%.4f", E_before, E_after)
                logm(logging.DEBUG, mod_id, "freq=%.4f prior_strength=%.4f", freq, prior_strength)
                logm(logging.DEBUG, mod_id, "alpha_prior=%.4f beta_prior=%.4f", alpha_prior, beta_prior)
                logm(logging.DEBUG, mod_id, "posterior_before alpha=%.4f beta=%.4f", alpha_post_before, beta_post_before)
                logm(logging.DEBUG, mod_id, "posterior_after  alpha=%.4f beta=%.4f", alpha_post_after, beta_post_after)
                logm(logging.DEBUG, mod_id, "beta_threshold=%.4f", beta_threshold)
                logm(logging.DEBUG, mod_id, "confidence_before=%.3f confidence_after=%.3f", conf_before, conf_after)
                logm(logging.DEBUG, mod_id, "posterior_mean=%.4f posterior_var=%.6f ci=[%.4f, %.4f]", credible_stats["posterior_mean"], credible_stats["posterior_variance"], credible_stats["ci_low"], credible_stats["ci_high"])
             
        
        modules_df = pd.DataFrame(rows_modules).sort_values(["module"]).reset_index(drop=True)
        
        
        return modules_df
        
    @staticmethod
    def evaluate_multiline_step_probabilities(
        module_dict_multiline: dict,
        df: pd.DataFrame,
        before_col: str = "Dk",
        after_col: str = "Dk_Neighbor",
        step_format: str = "path.step",   # "path.step" or "path_step"
        verbose: bool = False ) -> pd.DataFrame:
        CMP = CalculateModuleProbabilities
        Dk, Dk_Neighbor = CMP.build_dk_maps_from_df(df, before_col, after_col)

        rows_steps = []

        for mod_id, entry in module_dict_multiline.items():
            lines = entry.get("lines", [])
            if not lines:
                continue

            for line_obj in lines:
                path = int(line_obj.get("line", 0))  # your JSON uses "line"
                steps = line_obj.get("steps", [])

                for s in steps:
                    step_idx = int(s["step"])
                    eqn = s["equation"]

                    p_b = CMP.eval_equation(eqn, Dk)
                    p_a = CMP.eval_equation(eqn, Dk_Neighbor)

                    if step_format == "path_step":
                        step_label = f"{path}_{step_idx}"
                    else:
                        step_label = f"{path}.{step_idx}"

                    rows_steps.append({
                        "module": mod_id,
                        "multiline": True,
                        "step": step_label,
                        "equation": eqn,
                        "p_before": p_b,
                        "p_after":  p_a,
                    })

                    if verbose:
                        #logging.getLogger().info(f"--- Module {mod_id} ---")
                        kos = sorted(set(re.findall(r"K\d{5}", eqn)))
                        logging.debug("module=%s multiline=%s step=%s equation=%s", mod_id, True, step_label, eqn)

                        for ko in kos:
                            logging.debug(
                                "module=%s ko=%s Dk=%.4f Dk_Neighbor=%.4f", mod_id,
                                ko, Dk.get(ko, 0.0), Dk_Neighbor.get(ko, 0.0)
                            )
                        logging.debug(
                            "       p_before=%.6f  p_after=%.6f",
                            p_b, p_a
                        )


        steps_df = pd.DataFrame(rows_steps)

        if steps_df.empty:
            # return empty with expected columns for safety
            return pd.DataFrame(columns=["module","multiline","step","equation","p_before","p_after"])
        def _step_key(x):
            try:
                if isinstance(x, str) and "." in x:
                    a,b = x.split(".", 1)
                    return (int(a), int(b))
                if isinstance(x, str) and "_" in x:
                    a,b = x.split("_", 1)
                    return (int(a), int(b))
            except Exception:
                pass
            return (10**9, 10**9)

        steps_df["_sort"] = steps_df["step"].map(_step_key)
        steps_df = steps_df.sort_values(["module","_sort"]).drop(columns=["_sort"]).reset_index(drop=True)
        return steps_df




    @staticmethod
    def calculate_multiline_module_confidence_from_steps(
        multiline_steps_df: pd.DataFrame,
        module_dict_multiline: dict,
        genome_completeness: float,
        module_frequencies: dict = None,
        default_frequency: float = 0.5,
        prior_strength: float = 1.0,
        default_beta_thresh: float = 0.65,
        verbose: bool = False
    ) -> pd.DataFrame:

        CMP = CalculateModuleProbabilities
        if module_frequencies is None:
            module_frequencies = {}

        if multiline_steps_df.empty:
            return pd.DataFrame(columns=[
                "module","module_equation","n_steps","E_before","E_after",
                "module_frequency","module_probability_before","module_probability_after"
            ])

        def _parse_path(step_label):
            if isinstance(step_label, str) and "." in step_label:
                a, _ = step_label.split(".", 1)
                return int(a)
            if isinstance(step_label, str) and "_" in step_label:
                a, _ = step_label.split("_", 1)
                return int(a)
            # fallback: treat as path 1
            return 1

        df = multiline_steps_df.copy()
        df["path"] = df["step"].map(_parse_path)

        beta_threshold = default_beta_thresh * genome_completeness

        rows_modules = []

        for mod_id, gmod in df.groupby("module", sort=False):
            # module frequency guardrails (avoid alpha/beta degeneracy)
            freq = float(module_frequencies.get(mod_id, default_frequency))
            if freq <= 0.0: freq = 0.001
            if freq >= 1.0: freq = 0.999

            # compute confidence for each path independently
            path_rows = []
            for path_id, gpath in gmod.groupby("path", sort=True):
                E_before = float(gpath["p_before"].sum())
                E_after  = float(gpath["p_after"].sum())
                n_steps  = int(len(gpath))
                ps = n_steps * 0.1
                alpha_p = freq * ps
                beta_p  = (1 - freq) * ps
                a_post  = alpha_p + E_after
                b_post  = beta_p + (n_steps - E_after)
               

                conf_before = CMP.calculate_confidence(E_before, n_steps, freq, beta_threshold, ps)
                conf_after  = CMP.calculate_confidence(E_after,  n_steps, freq, beta_threshold, ps)
                cred = CMP.calculate_credible_stats(freq * ps + E_after,(1 - freq) * ps + (n_steps - E_after))
                hdi_l, hdi_h = CMP.calculate_hdi(a_post, b_post)
                eff = (cred['posterior_mean'] - beta_threshold) / math.sqrt(cred['posterior_variance']) if cred['posterior_variance'] > 0 else 0.0

                #path_status = CMP.classify_module(conf_after, cred["ci_low"], cred["ci_high"], hdi_l, hdi_h, beta_threshold, cred["posterior_mean"], n_steps, eff)



                path_rows.append({
                    "path": path_id,
                    "n_steps": n_steps,
                    "E_before": E_before,
                    "E_after": E_after,
                    "conf_before": conf_before,
                    "conf_after": conf_after,
                    "posterior_mean": cred["posterior_mean"],
                    "posterior_variance": cred["posterior_variance"],
                    "ci_low": cred["ci_low"],
                    "ci_high": cred["ci_high"],
                    "beta_threshold": beta_threshold,
                    "hdi_low": hdi_l,
                    "hdi_high": hdi_h,
                    "effect_size": eff,
                    #"Module_Status": path_status
                })
                if verbose:

                    logging.debug(
                        "[PATH] module=%s path=%s n_steps=%d",
                        mod_id, path_id, n_steps
                    )
                    logging.debug(
                        "       E_before=%.4f E_after=%.4f",
                        E_before, E_after
                    )
                    logging.debug(
                        "       conf_before=%.6f conf_after=%.6f",
                        conf_before, conf_after
                    )


            # choose winning path by max conf_after
            best = max(path_rows, key=lambda r: r["conf_after"])
            mod_eq = module_dict_multiline.get(mod_id, {}).get("module_equation", "")

            rows_modules.append({
                "module": mod_id,
                "module_equation": mod_eq,
                "n_steps": best["n_steps"],
                "E_before": best["E_before"],
                "E_after": best["E_after"],
                "module_frequency": freq,
                "module_probability_before": best["conf_before"],
                "module_probability_after":  best["conf_after"],
                "posterior_mean": best["posterior_mean"],
                "posterior_variance": best["posterior_variance"],
                "beta_threshold": beta_threshold,
                "ci_low": best["ci_low"],
                "ci_high": best["ci_high"],
                "hdi_low": best["hdi_low"],
                "hdi_high": best["hdi_high"],
                "effect_size": best["effect_size"],
            })
            if verbose:
                logging.debug(
                    "[PATH-SELECT] module=%s best_conf=%.6f",
                    mod_id, best["conf_after"]
                )


        return pd.DataFrame(rows_modules).sort_values("module").reset_index(drop=True)
    




class ModuleBestPath:
    class _Node: pass

    class _KO(_Node):
        def __init__(self, kid: str):
            self.kid = kid

    class _AND(_Node):
        def __init__(self, kids: List["_Node"]):
            self.kids = kids

    class _OR(_Node):
        def __init__(self, kids: List["_Node"]):
            self.kids = kids

    def __init__(
        self,
        module_eq: Dict[str, Any],
        ko_df: pd.DataFrame,
        *,
        ko_id_col: str = "KO id",
        ko_prob_col: str = "Dk_Neighbor",
        score_col: str = "score",
        keep_duplicate_kos: str = "max",
    ):
        self.module_eq = module_eq
        self.pKO, self.ko_score = self._build_pko_and_score(
            ko_df,
            ko_id_col=ko_id_col,
            ko_prob_col=ko_prob_col,
            score_col=score_col,
            keep=keep_duplicate_kos,
        )

    def run_all(self):
        rows = []
        failures = []

        for mid, payload in self.module_eq.items():

            # single-line modules
            for st in payload.get("steps", []):
                try:
                    node = self._parse_step_equation(st["equation"])
                    score, kos = self._eval_best_path(node)
                    rows.append({
                        "module": mid,
                        "multiline": False,
                        "step": st["step"],
                        "best_path_score": score,
                        "best_path_kos": ",".join(sorted(kos)),
                    })
                except Exception as e:
                    failures.append({
                        "module": mid,
                        "kind": "step_equation",
                        "step": st["step"],
                        "error": str(e),
                        "equation_head": st["equation"][:200],
                    })

            # multiline modules
            for line in payload.get("lines", []):
                line_no = line.get("line")
                for st in line.get("steps", []):
                    try:
                        node = self._parse_step_equation(st["equation"])
                        score, kos = self._eval_best_path(node)
                        rows.append({
                            "module": mid,
                            "multiline": True,
                            "step": f"{line_no}.{st['step']}",
                            "best_path_score": score,
                            "best_path_kos": ",".join(sorted(kos)),
                        })
                    except Exception as e:
                        failures.append({
                            "module": mid,
                            "kind": "line_step_equation",
                            "step": f"{line_no}.{st['step']}",
                            "error": str(e),
                            "equation_head": st["equation"][:200],
                        })


        return pd.DataFrame(rows), pd.DataFrame(failures)

    _KO_RE = re.compile(r"K\d{5}")

    def _parse_step_equation(self, expr: str) -> "_Node":
        s = expr.strip()
        s = self._unwrap_parens(s)

        # atomic KO
        if self._KO_RE.fullmatch(s):
            return self._KO(s)

        # max(A,B,...)  → OR
        if s.startswith("max(") and s.endswith(")"):
            inner = s[4:-1]
            parts = self._split_top_level(inner, ",")
            return self._OR([self._parse_step_equation(p) for p in parts])

        # noisy-OR: (1 - ((1 - A)*(1 - B)*...))
        if s.startswith("1 - ("):
            inner = self._unwrap_parens(s[4:])
            terms = self._split_top_level(inner, "*")
            kids = []
            for t in terms:
                t = self._unwrap_parens(t)
                if not t.startswith("1 - "):
                    raise ValueError(f"Invalid noisy-OR term: {t}")
                kids.append(self._parse_step_equation(t[4:].strip()))
            return self._OR(kids)

        # AND / sequence
        if "*" in s:
            parts = self._split_top_level(s, "*")
            return self._AND([self._parse_step_equation(p) for p in parts])

        raise ValueError(f"Unrecognized step equation format: {s}")

    def _unwrap_parens(self, s: str) -> str:
        s = s.strip()
        while s.startswith("(") and s.endswith(")"):
            depth = 0
            ok = True
            for i, ch in enumerate(s):
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
                if depth == 0 and i < len(s) - 1:
                    ok = False
                    break
            if not ok:
                break
            s = s[1:-1].strip()
        return s

    def _split_top_level(self, s: str, sep: str) -> List[str]:
        parts, buf, depth = [], [], 0
        for ch in s:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            if ch == sep and depth == 0:
                parts.append("".join(buf).strip())
                buf = []
            else:
                buf.append(ch)
        tail = "".join(buf).strip()
        if tail:
            parts.append(tail)
        return parts

    def _tie_score(self, kos: Set[str]) -> float:
        # Use max score among KOs used in that branch (usually 1 KO)
        best = float("-inf")
        for k in kos:
            best = max(best, self.ko_score.get(k, float("-inf")))
        return best


    def _eval_best_path(self, node: "_Node") -> Tuple[float, Set[str]]:
        if isinstance(node, self._KO):
            return float(self.pKO.get(node.kid, 0.0)), {node.kid}

        if isinstance(node, self._AND):
            score = 1.0
            path: Set[str] = set()
            for ch in node.kids:
                s, ks = self._eval_best_path(ch)
                score *= s
                path |= ks
            return score, path

        if isinstance(node, self._OR):
            best_s = -1.0
            best_path: Set[str] = set()
            best_tie = float("-inf")

            for ch in node.kids:
                s, ks = self._eval_best_path(ch)

                if s > best_s:
                    best_s, best_path = s, ks
                    best_tie = self._tie_score(ks) if abs(s - 1.0) < 1e-12 else float("-inf")
                    continue

                # tie on probability
                if abs(s - best_s) < 1e-12 and abs(s - 1.0) < 1e-12:
                    tie = self._tie_score(ks)
                    if tie > best_tie:
                        best_s, best_path, best_tie = s, ks, tie
                    elif abs(tie - best_tie) < 1e-12:
                        # deterministic fallback: pick lexicographically smallest KO-set
                        if sorted(ks) < sorted(best_path):
                            best_s, best_path, best_tie = s, ks, tie

            return (0.0 if best_s < 0 else best_s), best_path


        raise TypeError(type(node))


    def _build_pko_and_score(
        self,
        df: pd.DataFrame,
        *,
        ko_id_col: str,
        ko_prob_col: str,
        score_col: str,
        keep: str,
    ) -> Tuple[Dict[str, float], Dict[str, float]]:

        if ko_id_col not in df.columns:
            raise ValueError(f"KO id col not found: {ko_id_col!r}")
        if ko_prob_col not in df.columns:
            raise ValueError(f"KO prob col not found: {ko_prob_col!r}")
        if score_col not in df.columns:
            raise ValueError(f"Score col not found: {score_col!r}")

        tmp = df[[ko_id_col, ko_prob_col, score_col]].copy()
        tmp["KO_base"] = tmp[ko_id_col].astype(str).str.extract(r"(K\d{5})", expand=False)
        tmp["p"] = pd.to_numeric(tmp[ko_prob_col], errors="coerce")
        tmp["score"] = pd.to_numeric(tmp[score_col], errors="coerce")
        tmp = tmp.dropna(subset=["KO_base", "p"])

        # If score missing for some rows, treat as -inf so it never wins a tie
        tmp["score"] = tmp["score"].fillna(float("-inf"))

        if keep == "max":
            # probability: max over duplicates
            p_df = tmp.groupby("KO_base", as_index=False)["p"].max()
            # score: max over duplicates (best hit) — consistent with your HMM dedup logic
            s_df = tmp.groupby("KO_base", as_index=False)["score"].max()
        elif keep == "last":
            tmp = tmp.drop_duplicates("KO_base", keep="last")
            p_df = tmp[["KO_base", "p"]]
            s_df = tmp[["KO_base", "score"]]
        else:
            raise ValueError("keep_duplicate_kos must be 'max' or 'last'")

        pKO = dict(zip(p_df["KO_base"], p_df["p"].astype(float)))
        ko_score = dict(zip(s_df["KO_base"], s_df["score"].astype(float)))
        return pKO, ko_score
    

    @staticmethod
    def compute_module_best_paths(steps_df):
        module_best_rows = []

        for module, g in steps_df.groupby("module", sort=False):
            if not g["multiline"].any():
                all_kos = []

                for _, r in g.sort_values("step").iterrows():
                    if pd.notna(r["best_path_kos"]) and r["best_path_kos"] != "":
                        all_kos.extend(r["best_path_kos"].split(","))

                module_best_rows.append({"module": module,"module_best_path_kos": ",".join(sorted(set(all_kos))),})
            else:
                best_score = -1.0
                best_kos = ""
                best_ko_count = 10**9

                g = g.copy()
                g["line"] = g["step"].astype(str).str.split(".").str[0]

                for line, gl in g.groupby("line", sort=False):

                    score = 1.0
                    kos = []

                    for _, r in gl.iterrows():
                        score *= float(r["best_path_score"]) if pd.notna(r["best_path_score"]) else 0.0
                        if pd.notna(r["best_path_kos"]) and r["best_path_kos"] != "":
                            kos.extend(r["best_path_kos"].split(","))

                    kos = sorted(set(kos))
                    ko_count = len(kos)

                    if score > best_score:
                        best_score = score
                        best_kos = ",".join(kos)
                        best_ko_count = ko_count

                    elif score == best_score:
                        if ko_count < best_ko_count:
                            best_kos = ",".join(kos)
                            best_ko_count = ko_count

                module_best_rows.append({
                    "module": module,
                    "module_best_path_kos": best_kos,
                })

        return pd.DataFrame(module_best_rows)



class SubstitutionPipeline:

    @staticmethod
    def get_step_reaction_anchor(steps_df: pd.DataFrame, module: str, step: str):
        mask = (steps_df["module"] == module) & (steps_df["step"].astype(str) == str(step))
        rows = steps_df[mask]
        if rows.empty:
            return set()
        rxn_str = rows.iloc[0].get("best_path_reactions", "")
        if not rxn_str or pd.isna(rxn_str):
            return set()
        return {r.strip().upper() for r in str(rxn_str).split(",") if r.strip()}


    @staticmethod
    def compute_zscore(score: float, threshold: float) -> float:
        if threshold <= 0 or pd.isna(threshold) or pd.isna(score):
            return np.nan
        return (score - threshold) / threshold

    
    @staticmethod
    def find_substitutes(steps_df: pd.DataFrame, dk_full: pd.DataFrame, ko_to_reactions: Dict[str, Set[str]], z_threshold: float = -0.25) -> pd.DataFrame:
        """For each module step, find functional substitute KOs for:
          Case 1 - Outcompeted: canonical KO lost at locus to overlapgroup_winner
          Case 2 - Severely below threshold: score < z_threshold * kofam_score_threshold,substitute must be at same ORF
        """
        logging.info("Finding functional substitutes...")

        # All KOs with hit_conf > 0 (have real evidence)
        dk_with_hits = dk_full[dk_full["hit_conf"].fillna(0.0) > 0].copy()
        dk_with_hits["KO id"] = dk_with_hits["KO id"].astype(str).str.strip().str.upper()
        # overlap_group -> list of (ko, hit_conf) for all KOs with hits
        locus_ko_map: Dict[str, List[tuple]] = {}
        for _, row in dk_with_hits.iterrows():
            og = str(row.get("overlap_group", ""))
            if not og:
                continue
            locus_ko_map.setdefault(og, []).append((
                row["KO id"],
                float(row.get("hit_conf", 0.0)),
                float(row.get("score", 0.0)),))
            
        # KO -> overlap_group (for canonical KO locus lookup)
        ko_to_locus: Dict[str, str] = {}
        for _, row in dk_full.iterrows():
            ko = str(row.get("KO id", "")).strip().upper()
            og = str(row.get("overlap_group", ""))
            if ko and og:
                ko_to_locus[ko] = og

        # KO -> overlapgroup_winner
        ko_to_winner: Dict[str, str] = {}
        for _, row in dk_full.iterrows():
            ko = str(row.get("KO id", "")).strip().upper()
            winner = str(row.get("overlapgroup_winner", "")).strip().upper()
            if ko and winner:
                ko_to_winner[ko] = winner
        
        # KO -> is_outcompeted, flag_is_below_kofam_threshold, score, threshold
        ko_flags: Dict[str, dict] = {}
        for _, row in dk_full.iterrows():
            ko = str(row.get("KO id", "")).strip().upper()
            ko_flags[ko] = {
                "is_outcompeted":               bool(row.get("is_outcompeted", False)),
                "flag_is_below_kofam_threshold": bool(row.get("flag_is_below_kofam_threshold", False)),
                "score":                         float(row.get("score", 0.0)),
                "kofam_score_threshold":         float(row.get("kofam_score_threshold", np.nan))
                    if pd.notna(row.get("kofam_score_threshold")) else np.nan,
                "hit_conf":                      float(row.get("hit_conf", 0.0)),
                "winner":                        str(row.get("overlapgroup_winner", "")).strip().upper(),
            }

        # Walk every step in steps_df     
        records = []

        for _, step_row in steps_df.iterrows():
            module  = str(step_row["module"])
            step    = str(step_row["step"])
            equation = str(step_row.get("equation", ""))

            # Reaction anchor for this step
            anchor_reactions = SubstitutionPipeline.get_step_reaction_anchor(
                steps_df, module, step
            )
            if not anchor_reactions:
                continue

            # Canonical KOs in this step equation
            canonical_kos = set(re.findall(r"K\d{5}", equation))
            if not canonical_kos:
                continue

            # pre-compute all KOs in this module's equations once
            module_kos_all_steps = set(
                ko
                for _, other_step in steps_df[steps_df["module"] == module].iterrows()
                for ko in re.findall(r"K\d{5}", str(other_step.get("equation", "")))
            )

            for canonical_ko in canonical_kos:
                flags = ko_flags.get(canonical_ko)
                if flags is None:
                    # KO not in genome at all — skip
                    continue

                # Case 1: Outcompeted  
                if flags["is_outcompeted"]:
                    winner = flags["winner"]
                    if not winner or winner == canonical_ko:
                        continue

                    # Eq. 27 condition 1 (implicit): substitute must be at the SAME locus
                    # as the canonical enzyme, not just win some locus elsewhere
                    canonical_locus = ko_to_locus.get(canonical_ko)
                    winner_locus = ko_to_locus.get(winner)
                    if not canonical_locus or winner_locus != canonical_locus:
                        continue

                    # skip if substitute already evaluated anywhere in this module
                    if winner in module_kos_all_steps:
                        continue

                    # Eq. 27 condition 2: rxn(ko†) ∩ rxn(ko) ∩ rxn_s ≠ ∅
                    # — must overlap the CANONICAL enzyme's own reactions, not just the step's anchor set
                    canonical_reactions = ko_to_reactions.get(canonical_ko, set())
                    winner_reactions = ko_to_reactions.get(winner, set())
                    shared = canonical_reactions & winner_reactions & anchor_reactions

                    if shared:
                        records.append({
                            "module":              module,
                            "step":                step,
                            "canonical_ko":        canonical_ko,
                            "substitute_ko":       winner,
                            "case":                "outcompeted",
                            "shared_reactions":    ",".join(sorted(shared)),
                            "n_shared":            len(shared),
                            "z_score":             np.nan,
                            "substitute_hit_conf": float(
                                dk_full.loc[
                                    dk_full["KO id"].astype(str).str.upper() == winner,
                                    "hit_conf"
                                ].max() if (dk_full["KO id"].astype(str).str.upper() == winner).any()
                                else 0.0
                            ),
                        })

                # ---------------------------------------------------------- #
                # Case 2: Severely below threshold                           
                # elif flags["flag_is_below_kofam_threshold"]:
                #     z = SubstitutionPipeline.compute_zscore(
                #         flags["score"], flags["kofam_score_threshold"]
                #     )
                #     if np.isnan(z) or z >= z_threshold:
                #         continue

                #     # Same locus only
                #     locus = ko_to_locus.get(canonical_ko)
                #     if not locus:
                #         continue

                #     locus_candidates = locus_ko_map.get(locus, [])

                #     for (cand_ko, cand_hit_conf, cand_score) in locus_candidates:
                #         if cand_ko == canonical_ko:
                #             continue
                #         if cand_hit_conf <= 0:
                #             continue

                #         cand_reactions = ko_to_reactions.get(cand_ko, set())
                #         shared = anchor_reactions & cand_reactions

                #         if shared:
                #             records.append({
                #                 "module":              module,
                #                 "step":                step,
                #                 "canonical_ko":        canonical_ko,
                #                 "substitute_ko":       cand_ko,
                #                 "case":                "below_threshold",
                #                 "shared_reactions":    ",".join(sorted(shared)),
                #                 "n_shared":            len(shared),
                #                 "z_score":             round(z, 4),
                #                 "substitute_hit_conf": cand_hit_conf,
                #             })

        subs_df = pd.DataFrame(records) if records else pd.DataFrame(columns=[
            "module", "step", "canonical_ko", "substitute_ko",
            "case", "shared_reactions", "n_shared", "z_score", "substitute_hit_conf"
        ])

        print(f"Found {len(subs_df)} functional substitute(s) across " f"{subs_df['module'].nunique() if not subs_df.empty else 0} module(s).\n")

        return subs_df
    @staticmethod
    def _norm_step_key(x) -> str:
        """Canonicalize a step identifier for comparison across sources.

        subs_df's "step" comes from steps_df (already int-cast upstream, see
        evaluate_step_probabilities: idx = int(s["step"])), while module_eq's
        own s["step"] is whatever numeric type survived a raw json.load —
        JSON does not distinguish 8 from 8.0, so Python can hand back either
        an int or a float depending on how the source file was serialized.
        str(8) != str(8.0), so a plain string comparison between the two
        sides silently fails and a substitute never gets injected even
        though find_substitutes correctly accepted it. Route both sides
        through float() first so "8", 8, and 8.0 all normalize the same way,
        and only fall back to a raw string compare for genuinely
        non-numeric step ids (e.g. multiline "line.step" composites, which
        get their own dotted-string comparison below and never reach here).
        """
        try:
            f = float(x)
            return str(int(f)) if f.is_integer() else str(f)
        except (TypeError, ValueError):
            return str(x).strip()

    @staticmethod
    def build_augmented_equations(
        module_eq: dict,
        subs_df: pd.DataFrame,
    ) -> dict:
        if subs_df.empty:
            return module_eq

        aug_eq = copy.deepcopy(module_eq)

        # Group substitutes by (module, step)
        for (module, step), grp in subs_df.groupby(["module", "step"]):
            module = str(module)
            step   = str(step)

            if module not in aug_eq:
                continue

            entry = aug_eq[module]

            # canonical -> sorted unique substitutes, for THIS step only
            sub_map = {}
            for _, r in grp.iterrows():
                can = str(r["canonical_ko"]).strip().upper()
                sub = str(r["substitute_ko"]).strip().upper()
                if can == sub:
                    continue
                sub_map.setdefault(can, set()).add(sub)
            sub_map = {k: sorted(v) for k, v in sub_map.items()}

            if not sub_map:
                continue

            # ---- single-line steps ----
            for s in entry.get("steps", []):
                if str(s["step"]) == step:
                    s["equation"] = SubstitutionPipeline._inject_substitutes(
                        s["equation"], sub_map
                    )

            # ---- multiline steps ----
            for line in entry.get("lines", []):
                for s in line.get("steps", []):
                    if str(f"{line['line']}.{s['step']}") == step:
                        s["equation"] = SubstitutionPipeline._inject_substitutes(
                            s["equation"], sub_map
                        )

        return aug_eq

    
    @staticmethod
    def _inject_substitutes(equation: str, sub_map: dict) -> str:

        if not sub_map:
            return equation

        # normalise keys/values
        clean = {}
        for can, subs in sub_map.items():
            can = str(can).strip().upper()
            uniq = sorted({str(s).strip().upper() for s in subs} - {can})
            if uniq:
                clean[can] = uniq
        if not clean:
            return equation

        def _repl(m):
            ko = m.group(0).upper()
            subs = clean.get(ko)
            if not subs:
                return m.group(0)          # untouched
            return f"max({ko}," + ",".join(subs) + ")"
        return re.sub(r"K\d{5}", _repl, equation)

    @staticmethod
    def get_substitute_initial_details(
        subs_df: pd.DataFrame,
        dk_full: pd.DataFrame,
        neighbor_map: dict,
        ko_counts: dict,
        cooccurrence_lookup: dict,
    ) -> pd.DataFrame:

        # index dk_full by KO id for fast lookup
        ko_idx = (
            dk_full
            .drop_duplicates(subset="KO id", keep="first")
            .set_index("KO id")
        )

        rows = []
        for _, row in subs_df.iterrows():
            canonical  = str(row["canonical_ko"]).strip().upper()
            substitute = str(row["substitute_ko"]).strip().upper()
            module     = str(row["module"])
            step       = str(row["step"])
            case       = str(row["case"])

            # substitute's current values
            sub_row = ko_idx.loc[substitute] if substitute in ko_idx.index else None
            can_row = ko_idx.loc[canonical]  if canonical  in ko_idx.index else None

            sub_dk          = float(sub_row["Dk"])          if sub_row is not None else None
            sub_dk_neighbor = float(sub_row["Dk_Neighbor"]) if sub_row is not None else None
            sub_hit_conf    = float(sub_row["hit_conf"])    if sub_row is not None else None
            sub_ko_freq     = float(sub_row["KO_freq"])     if sub_row is not None else None

            can_dk          = float(can_row["Dk"])          if can_row is not None else None
            can_dk_neighbor = float(can_row["Dk_Neighbor"]) if can_row is not None else None

            # Nij between substitute and canonical directly
            nij_sub_can = cooccurrence_lookup.get(substitute, {}).get(canonical, 0.0)
            nij_can_sub = cooccurrence_lookup.get(canonical,  {}).get(substitute, 0.0)

            # canonical's current neighbors in neighbor_map
            can_neighbors = set(neighbor_map.get(canonical, {}).keys())
            n_can_neighbors = len(can_neighbors)

            # how many of canonical's neighbors have real Nij with substitute
            n_with_real_nij = sum(
                1 for nb in can_neighbors
                if cooccurrence_lookup.get(substitute, {}).get(nb, 0.0) > 0.0
            )

            # does substitute already have ANY edges to canonical's neighbors
            sub_existing_neighbors = set(neighbor_map.get(substitute, {}).keys())
            n_already_connected = len(sub_existing_neighbors & can_neighbors)

            rows.append({
                "module":               module,
                "step":                 step,
                "case":                 case,
                "canonical_ko":         canonical,
                "substitute_ko":        substitute,
                # substitute current state
                "sub_Dk":               sub_dk,
                "sub_Dk_Neighbor":      sub_dk_neighbor,
                "sub_hit_conf":         sub_hit_conf,
                "sub_KO_freq":          sub_ko_freq,
                # canonical current state
                "can_Dk":               can_dk,
                "can_Dk_Neighbor":      can_dk_neighbor,
                # co-occurrence between the pair
                "Nij_sub_to_can":       nij_sub_can,
                "Nij_can_to_sub":       nij_can_sub,
                # neighbor overlap
                "n_canonical_neighbors":        n_can_neighbors,
                "n_with_real_nij_in_lookup":    n_with_real_nij,
                "n_already_in_neighbor_map":    n_already_connected,
            })

        return pd.DataFrame(rows)

    @staticmethod
    def augment_neighbor_map(
        neighbor_map: dict,
        ko_counts: dict,
        subs_df: pd.DataFrame,
        cooccurrence_lookup: dict,
        default_nij: float = 0.0,
    ) -> tuple[dict, pd.DataFrame]:

            aug_map = copy.deepcopy(neighbor_map)
            pairs = (
                subs_df[["canonical_ko", "substitute_ko"]]
                .drop_duplicates()
                .values.tolist()
            )
            audit_rows = []

            for canonical, substitute in pairs:
                canonical  = str(canonical).strip().upper()
                substitute = str(substitute).strip().upper()
                can_neighbors = neighbor_map.get(canonical, {})

                if not can_neighbors:
                    logging.debug(
                        "SubstitutionPipeline: canonical %s has no neighbors, skipping.",
                        canonical
                    )
                    continue

                all_neighbors = set(can_neighbors.keys()) | {canonical}

                for nb in all_neighbors:
                    if nb == substitute:
                        continue

                    nij_fwd = cooccurrence_lookup.get(substitute, {}).get(nb, default_nij)

                    aug_map.setdefault(substitute, {})
                    if nb not in aug_map[substitute]:
                        aug_map[substitute][nb] = nij_fwd


                    audit_rows.append({
                        "canonical_ko":     canonical,
                        "substitute_ko":    substitute,
                        "neighbor":         nb,
                        "nij_fwd":          nij_fwd,
                        "nij_fwd_real":     nij_fwd > 0.0,
                        "nb_in_ko_counts":  nb in ko_counts,
                        "sub_in_ko_counts": substitute in ko_counts,
                    })

            audit_df = pd.DataFrame(audit_rows)

            if not audit_df.empty:
                n_pairs   = len(pairs)
                n_edges   = len(audit_df)
                n_real    = int(audit_df["nij_fwd_real"].sum())
                n_default = n_edges - n_real
                logging.info(
                    "SubstitutionPipeline: %d pairs | %d edges added | "
                    "%d real Nij | %d defaulted to 0.0",
                    n_pairs, n_edges, n_real, n_default
                )
                logging.info(
                    "Neighbor map augmented: %d substitute–canonical pairs | "
                    "%d edges added to the neighbor map | "
                    "%d edges have real co-occurrence counts | "
                    "%d edges defaulted to Nij=0.0 (no observed co-occurrence)",
                    n_pairs, n_edges, n_real, n_default
                )
            else:
                logging.info("SubstitutionPipeline: no edges added")


            return aug_map, audit_df

    @staticmethod
    def summarise_audit(audit_df: pd.DataFrame) -> pd.DataFrame:
        """
        Per (canonical, substitute) summary of the augmentation.
        Useful for the output CSV and debugging.
        """
        if audit_df.empty:
            return pd.DataFrame()
 
        return (
            audit_df
            .groupby(["canonical_ko", "substitute_ko"])
            .agg(
                n_neighbors_inherited = ("neighbor",      "count"),
                n_with_real_nij_fwd   = ("nij_fwd_real", "sum"),
                mean_nij_fwd          = ("nij_fwd",      "mean"),
                max_nij_fwd           = ("nij_fwd",      "max"),
            )
            .reset_index()
        )




    @staticmethod
    def run(
        *,
        steps_df: pd.DataFrame,
        dk_full: pd.DataFrame,
        dk_calculations: pd.DataFrame,
        dk_update_calculations: pd.DataFrame,
        neighbor_map: dict,
        ko_counts: dict,
        ko_to_reactions: Dict[str, Set[str]],
        cooc_json,
        single_only: dict,
        multi_only: dict,
        module_eq: dict,                    # needed to build augmented equations
        genome_completeness: float,         # sigma, for module confidence
        module_frequencies: dict,           # module priors, for module confidence
        module_reaction_map: dict = None,   # for steps best_path_reactions
        module_descriptions: dict = None,   # for steps module_description
        modules_df_pre: pd.DataFrame = None,  # pass-1 modules, for with/without compare
        output_prefix: str,
        verbose: bool = False,
        z_threshold: float = -0.25, ) -> dict:


        # Find functional substitutes using pass-1 best_path_reactions
        subs_df = SubstitutionPipeline.find_substitutes(
            steps_df=steps_df,
            dk_full=dk_full,           # full unfiltered KO set
            ko_to_reactions=ko_to_reactions,
            z_threshold=z_threshold,
        )
        #bare "subs_df.csv" dumped 
        subs_df.to_csv(f"{output_prefix}_substitutes.csv", index=False)
        logging.info("[SUBST] Wrote %s (%d substitute pairs)",
                     f"{output_prefix}_substitutes.csv", len(subs_df))

        sub_details = None
        neighbor_map_aug = neighbor_map
        steps_v2 = None
        modules_v2 = None          #init
        module_best_v2 = None      #init
        module_eq_aug = module_eq  #init

        # GATE 1 diagnostic (why no output files appear)
        if subs_df.empty:
            logging.warning(
                "[SUBST] GATE 1 BLOCKED: no substitutes found. "
                "Nothing after subs_df.csv will be written."
            )

        # Build the set of KOs we actually need lookups for:
        # substitute KOs + canonical KOs (to find their neighbors' Nij with substitute)
        needed_kos = set()
        if not subs_df.empty:
            needed_kos.update(subs_df["substitute_ko"].str.upper())
            needed_kos.update(subs_df["canonical_ko"].str.upper())
            # also need neighbors of canonicals
            for canonical in subs_df["canonical_ko"].str.upper():
                needed_kos.update(neighbor_map.get(canonical, {}).keys())

        cooccurrence_lookup = File_Helpers.load_cooccurrence_lookup_targeted(cooc_json, needed_kos)

        # GATE 2 diagnostic (why no output files appear)
        logging.info("[SUBST] cooc_json path = %s", cooc_json)
        logging.info("[SUBST] cooc file exists = %s",
                     Path(cooc_json).exists() if cooc_json else False)
        logging.info("[SUBST] cooccurrence_lookup entries = %d",
                     len(cooccurrence_lookup) if cooccurrence_lookup else 0)
        if not subs_df.empty and not cooccurrence_lookup:
            logging.warning(
                "[SUBST] GATE 2 BLOCKED: co-occurrence lookup is empty. "
                "Pass 2 is skipped entirely -- no _substitute_v2_check.csv, "
                "no _steps_v2.csv, no _modules_v2.csv will be written."
            )

        if not subs_df.empty:
            sub_details = SubstitutionPipeline.get_substitute_initial_details(
                subs_df=subs_df,
                dk_full=dk_full,
                neighbor_map=neighbor_map,
                ko_counts=ko_counts,
                cooccurrence_lookup=cooccurrence_lookup,
            )
            logging.info(
                "Substitute details: %d substitutes across %d modules",
                len(sub_details), sub_details["module"].nunique()
            )
            if verbose:
                logging.debug("Substitute details head:\n%s", sub_details.head())

            # only if write intermediates is true
            if verbose:
                sub_details_path = f"{output_prefix}_substitute_details.csv"
                sub_details.to_csv(sub_details_path, index=False)
                logging.info("Substitute details written to %s", sub_details_path)

            if cooccurrence_lookup:
                neighbor_map_aug, aug_audit_df = SubstitutionPipeline.augment_neighbor_map(
                    neighbor_map=neighbor_map,
                    ko_counts=ko_counts,
                    subs_df=subs_df,
                    cooccurrence_lookup=cooccurrence_lookup,
                )

                # intermediate, only with write_intermediates
                aug_summary = SubstitutionPipeline.summarise_audit(aug_audit_df)
                if verbose:
                    aug_summary.to_csv(f"{output_prefix}_substitute_augmentation_audit.csv", index=False)
            else:
                neighbor_map_aug = neighbor_map
                logging.info("No substitutes or co-occurrence lookup, skipping augmentation")

            # Pass 2: update Dk_Neighbor for substitute KOs
            if cooccurrence_lookup:

                dk_v2, used_neighbors_v2 = CalculateKOProbabilities.dk_neighbor_update(
                    dk_calculations,
                    neighbor_map_aug,
                    ko_counts,
                    alpha=0.6,
                    return_used=True,
                    verbose=verbose,
                )

                dk_v2 = dk_v2.rename(columns={"Dk_Neighbor": "Dk_Neighbor_v2"})

                # Record the post-augmentation associate set per KO.
                # For a substitute this is its own associates UNION the ones
                # inherited from the canonical it substitutes (whatever survived
                # the R(i,j) reliability filter inside dk_neighbor_update).
                dk_v2["KO_Neighbors_v2"] = dk_v2["KO id"].map(lambda k: ",".join(used_neighbors_v2.get(k, [])))
                dk_v2["KO_Neighbor_Count_v2"] = dk_v2["KO id"].map(lambda k: len(used_neighbors_v2.get(k, []))).fillna(0).astype(int)

                dk_full = dk_full.merge(
                    dk_v2[["KO id", "Dk_Neighbor_v2", "KO_Neighbors_v2", "KO_Neighbor_Count_v2"]], on="KO id", how="left",)
                dk_full["Dk_Neighbor_v2"] = dk_full["Dk_Neighbor_v2"].fillna(dk_full["Dk_Neighbor"])

                logging.info(
                    "Pass 2 max abs diff (Dk_Neighbor_v2 - Dk_Neighbor): %.6f",
                    (dk_full["Dk_Neighbor_v2"] - dk_full["Dk_Neighbor"]).abs().max()
                )

                sub_kos = set(subs_df["substitute_ko"].str.upper())
                # Superseded by _BLIMMP_substituted_dk.csv, which
                # carries all of these columns plus the full original schema.
                if verbose:
                    dk_full[dk_full["KO id"].str.upper().isin(sub_kos)][
                        ["KO id", "Dk", "Dk_Neighbor", "Dk_Neighbor_v2",
                         "KO_Neighbors_v2", "KO_Neighbor_Count_v2"]
                    ].to_csv(f"{output_prefix}_substitute_v2_check.csv", index=False)

                # Pass 2: re-run steps for affected modules
                # Inject substitutes into step equations as
                # max(canonical, substitute) so the substitute can actually
                # compete. Without this the substitute's updated Dk has
                # nowhere to go and every downstream number is unchanged.
                module_eq_aug = SubstitutionPipeline.build_augmented_equations(
                    module_eq, subs_df
                )
                single_aug = {m: e for m, e in module_eq_aug.items() if e.get("steps")}
                multi_aug  = {m: e for m, e in module_eq_aug.items() if e.get("lines")}

                affected_modules = set(subs_df["module"].unique())
                # Filter the AUGMENTED equations, not single_only/multi_only
                single_affected  = {m: e for m, e in single_aug.items() if m in affected_modules}
                multi_affected   = {m: e for m, e in multi_aug.items()  if m in affected_modules}

                logging.info("Pass 2: re-evaluating %d affected modules (%d single, %d multiline)",len(affected_modules), len(single_affected), len(multi_affected))

                substitute_rows = dk_full[dk_full["KO id"].str.upper().isin(sub_kos)].copy()
                dk_update_calculations = dk_update_calculations.merge(dk_v2[["KO id", "Dk_Neighbor_v2"]], on="KO id", how="left",)
                dk_update_calculations["Dk_Neighbor_v2"] = (dk_update_calculations["Dk_Neighbor_v2"].fillna(dk_update_calculations["Dk_Neighbor"]))
                dk_for_v2_eval = pd.concat([dk_update_calculations, substitute_rows], ignore_index=True,).drop_duplicates(subset="KO id", keep="first")

                steps_v2_single = CalculateModuleProbabilities.evaluate_step_probabilities( single_affected, dk_for_v2_eval, before_col="Dk", after_col="Dk_Neighbor_v2", verbose=verbose,)
                steps_v2_multi = CalculateModuleProbabilities.evaluate_multiline_step_probabilities(
                    multi_affected, dk_for_v2_eval,
                    before_col="Dk", after_col="Dk_Neighbor_v2",
                    step_format="path.step",
                    verbose=verbose,
                )
                steps_v2 = pd.concat([steps_v2_single, steps_v2_multi], ignore_index=True)
                steps_v2["step"] = steps_v2["step"].astype(str)

                # This is what lets a substitute actually win its step and
                # appear in the reported best path. Note ko_prob_col is
                # Dk_Neighbor_v2 here (pass 1 used Dk_Neighbor).
                best_paths_v2 = ModuleBestPath(
                    module_eq=module_eq_aug,
                    ko_df=dk_for_v2_eval,
                    ko_id_col="KO id",
                    ko_prob_col="Dk_Neighbor_v2",
                    score_col="score",
                    keep_duplicate_kos="max",
                )
                best_steps_v2, best_fails_v2 = best_paths_v2.run_all()
                best_steps_v2["step"] = best_steps_v2["step"].astype(str)

                steps_v2 = steps_v2.merge(
                    best_steps_v2[["module", "multiline", "step",
                                   "best_path_score", "best_path_kos"]],
                    on=["module", "multiline", "step"], how="left",
                )
                module_best_v2 = ModuleBestPath.compute_module_best_paths(steps_v2)

                # Module probabilities AFTER substitution
                # BUGFIX (Aug 2026): cast to bool first -- "multiline" arrives as object
                # dtype, and ~ on object does bitwise negation (True -> -1),
                # which blows up as a column indexer.
                _is_multi = steps_v2["multiline"].astype(bool)
                modules_v2_single = CalculateModuleProbabilities.calculate_module_confidence(
                    steps_df=steps_v2[~_is_multi],
                    module_dict=single_affected,
                    genome_completeness=genome_completeness,
                    module_frequencies=module_frequencies,
                    verbose=verbose,
                )
                modules_v2_multi = CalculateModuleProbabilities.calculate_multiline_module_confidence_from_steps(
                    steps_v2[_is_multi],
                    module_dict_multiline=multi_affected,
                    genome_completeness=genome_completeness,
                    module_frequencies=module_frequencies,
                    verbose=verbose,
                )
                modules_v2 = pd.concat(
                    [modules_v2_single, modules_v2_multi], ignore_index=True
                )
                modules_v2 = modules_v2.merge(module_best_v2, on="module", how="left")

                # Module probabilities WITH vs WITHOUT substitution
                # modules_v2 alone only shows the "with" side. Merge the pass-1
                # numbers so each module row carries both, plus the deltas.
                modules_cmp = modules_v2.copy()
                modules_cmp["module_description"] = (
                    modules_cmp["module"].map(module_descriptions).fillna("")
                    if module_descriptions is not None else ""
                )

                if modules_df_pre is not None and not modules_df_pre.empty:
                    _pre_m = modules_df_pre.set_index("module")

                    def _pre_col(col):
                        if col not in _pre_m.columns:
                            logging.warning(
                                "[SUBST] pass-1 modules_df has no '%s'; "
                                "with/without comparison for it will be blank.", col
                            )
                            return {}
                        return _pre_m[col].to_dict()

                    _m = modules_cmp["module"]
                    modules_cmp["module_confidence_no_sub"]  = _m.map(_pre_col("module_probability_after"))
                    modules_cmp["posterior_mean_no_sub"]     = _m.map(_pre_col("posterior_mean"))
                    modules_cmp["ci_low_no_sub"]             = _m.map(_pre_col("ci_low"))
                    modules_cmp["ci_high_no_sub"]            = _m.map(_pre_col("ci_high"))
                    modules_cmp["E_after_no_sub"]            = _m.map(_pre_col("E_after"))
                    modules_cmp["best_path_no_sub"]          = _m.map(_pre_col("module_best_path_kos"))

                    modules_cmp["module_confidence_with_sub"] = modules_cmp["module_probability_after"]
                    modules_cmp["module_confidence_delta"] = (
                        modules_cmp["module_confidence_with_sub"]
                        - modules_cmp["module_confidence_no_sub"]
                    )
                    modules_cmp["module_confidence_improved"] = (
                        modules_cmp["module_confidence_delta"] > 1e-12
                    )
                    modules_cmp["E_after_delta"] = (
                        modules_cmp["E_after"] - modules_cmp["E_after_no_sub"]
                    )
                    modules_cmp["best_path_changed"] = [
                        set(str(a or "").upper().split(",")) != set(str(b or "").upper().split(","))
                        for a, b in zip(
                            modules_cmp["best_path_no_sub"],
                            modules_cmp["module_best_path_kos"],
                        )
                    ]
                else:
                    logging.warning(
                        "[SUBST] modules_df_pre not supplied; "
                        "_module_probabilities_sub_v2.csv will have no 'no_sub' columns."
                    )

                _cmp_cols = [c for c in [
                    "module", "module_description", "n_steps",
                    "module_confidence_no_sub", "module_confidence_with_sub",
                    "module_confidence_delta", "module_confidence_improved",
                    "E_after_no_sub", "E_after", "E_after_delta",
                    "best_path_no_sub", "module_best_path_kos", "best_path_changed",
                    "posterior_mean_no_sub", "posterior_mean",
                    "ci_low_no_sub", "ci_low", "ci_high_no_sub", "ci_high",
                    "beta_threshold", "effect_size", "module_frequency",
                ] if c in modules_cmp.columns]
                _cmp_path = f"{output_prefix}_module_probabilities_sub_v2.csv"
                modules_cmp[_cmp_cols].to_csv(_cmp_path, index=False)
                logging.info("[SUBST] module with/without comparison -> %s", _cmp_path)

                # Superseded by _BLIMMP_substituted_module_steps.csv
                # and _module_probabilities_sub_v2.csv respectively.
                if verbose:
                    steps_v2.to_csv(f"{output_prefix}_steps_v2.csv", index=False)
                    modules_v2.to_csv(f"{output_prefix}_modules_v2.csv", index=False)

                # Dk CSV in the SAME schema as the original ──
                # Base is dk_full (unfiltered) so substitute KOs that belong
                # to no KEGG module are still present.
                dk_v2_out = dk_full.copy()

                # Modules column lives on dk_update_calculations, not dk_full
                _mod_map = dict(zip(
                    dk_update_calculations["KO id"],
                    dk_update_calculations.get("Modules", pd.Series(dtype=object)),
                ))
                dk_v2_out["Modules"] = dk_v2_out["KO id"].map(_mod_map)

                # the substituted run's final values become the "final" columns
                dk_v2_out["ko_probability_pass1"] = dk_v2_out["Dk_Neighbor"]
                dk_v2_out["Dk_Neighbor"] = dk_v2_out["Dk_Neighbor_v2"]
                
                dk_v2_out["KO_Neighbors"] = dk_v2_out["KO_Neighbors_v2"]
                dk_v2_out["KO_Neighbor_Count"] = dk_v2_out["KO_Neighbor_Count_v2"]

                # KO_Neighbors above is only what survived the R(i,j) filter in
                # dk_neighbor_update. An associate is dropped when Nj <= 0, or
                # when Nij/Nj < 0.25 (the min_frac floor in
                # calculate_reliable_conditional_prob). These columns expose the
                # full candidate set from the augmented map and what fell out.
                _ko_upper = dk_v2_out["KO id"].astype(str).str.upper()

                _cand_map = {
                    k: sorted(set(v.keys()) - {k})
                    for k, v in neighbor_map_aug.items()
                }
                dk_v2_out["KO_Neighbors_all"] = _ko_upper.map(
                    lambda k: ",".join(_cand_map.get(k, []))
                )
                dk_v2_out["KO_Neighbor_All_Count"] = _ko_upper.map(
                    lambda k: len(_cand_map.get(k, []))
                ).fillna(0).astype(int)

                dk_v2_out["KO_Neighbors_dropped"] = _ko_upper.map(
                    lambda k: ",".join(sorted(
                        set(_cand_map.get(k, [])) - set(used_neighbors_v2.get(k, []))
                    ))
                )
                dk_v2_out["KO_Neighbor_Dropped_Count"] = _ko_upper.map(
                    lambda k: len(
                        set(_cand_map.get(k, [])) - set(used_neighbors_v2.get(k, []))
                    )
                ).fillna(0).astype(int)

                def _drop_reasons(i_ko):
                    """Per-dropped-associate reason string: j(Nij/Nj=frac,reason)."""
                    out = []
                    cands = _cand_map.get(i_ko, [])
                    kept = set(used_neighbors_v2.get(i_ko, []))
                    for j in cands:
                        if j in kept:
                            continue
                        nij = neighbor_map_aug.get(i_ko, {}).get(j, 0.0)
                        nj = ko_counts.get(j, 0.0)
                        if nj <= 0:
                            out.append(f"{j}(Nj=0)")
                        else:
                            frac = nij / nj
                            out.append(f"{j}(Nij/Nj={frac:.3f}<0.25)")
                    return ";".join(out)

                dk_v2_out["KO_Neighbor_Drop_Reasons"] = _ko_upper.map(_drop_reasons)

                # Associates added by augmentation (in neighbor_map_aug but not original neighbor_map)
                _orig_assoc = {k: set(v.keys()) for k, v in neighbor_map.items()}
                _aug_assoc  = {k: set(v.keys()) for k, v in neighbor_map_aug.items()}

                dk_v2_out["substituted_associates"] = _ko_upper.map(lambda k: ",".join(sorted( _aug_assoc.get(k, set()) - _orig_assoc.get(k, set()) - {k})))
                dk_v2_out["substituted_associate_count"] = _ko_upper.map(lambda k: len(_aug_assoc.get(k, set()) - _orig_assoc.get(k, set()) - {k})).fillna(0).astype(int)

                # the two new columns
                _sub_map = (
                    subs_df.groupby(subs_df["substitute_ko"].str.upper())["canonical_ko"]
                    .apply(lambda s: ",".join(sorted(set(s.str.upper()))))
                    .to_dict()
                )
                dk_v2_out["is_substitute"] = _ko_upper.isin(sub_kos)
                dk_v2_out["substitutes_for"] = _ko_upper.map(_sub_map).fillna("")

                # Enrich steps_v2 so it matches the ORIGINAL
                # steps schema. Adds the two columns the original has that
                # the pass-2 rebuild does not produce on its own, plus two
                # substitution markers mirroring the Dk file.
                _base_cols = ["module", "multiline", "step", "equation",
                              "p_before", "p_after", "best_path_score", "best_path_kos"]
                _base_steps = steps_df[[c for c in _base_cols if c in steps_df.columns]].copy()
                _base_steps["step"] = _base_steps["step"].astype(str)

                _v2 = steps_v2.copy()
                _v2["step"] = _v2["step"].astype(str)

                _base_key = _base_steps["module"].astype(str) + "|" + _base_steps["step"]
                _v2_key = _v2["module"].astype(str) + "|" + _v2["step"]

                steps_v2_out = pd.concat(
                    [_base_steps[~_base_key.isin(set(_v2_key))], _v2],
                    ignore_index=True,
                ).sort_values(["module", "step"]).reset_index(drop=True)

                if module_reaction_map is not None:
                    steps_v2_out["best_path_reactions"] = steps_v2_out.apply(
                        lambda r: File_Helpers.reactions_for_module_bestpath(
                            module_reaction_map, r["module"], r["best_path_kos"]
                        ),
                        axis=1,
                    )
                else:
                    steps_v2_out["best_path_reactions"] = ""
                if module_descriptions is not None:
                    steps_v2_out["module_description"] = (
                        steps_v2_out["module"].map(module_descriptions).fillna("")
                    )
                else:
                    steps_v2_out["module_description"] = ""

                # which steps had a substitute injected, and which KO(s)
                # which steps had a substitute injected, and which KO(s)
                _step_sub_map = (
                    subs_df.assign(_k=subs_df["module"].astype(str) + "|" + subs_df["step"].astype(str))
                    .groupby("_k")["substitute_ko"]
                    .apply(lambda s: ",".join(sorted(set(s.str.upper()))))
                    .to_dict()
                )
                # ...and which canonical KO each substitute stood in for.
                # Without this, once a substitute wins its step, best_path_kos
                # and substitute_kos both show only the winner, and the
                # displaced canonical enzyme disappears from the row entirely.
                _step_canonical_map = (
                    subs_df.assign(_k=subs_df["module"].astype(str) + "|" + subs_df["step"].astype(str))
                    .groupby("_k")["canonical_ko"]
                    .apply(lambda s: ",".join(sorted(set(s.str.upper()))))
                    .to_dict()
                )
                _step_key = (
                    steps_v2_out["module"].astype(str) + "|" + steps_v2_out["step"].astype(str)
                )
                steps_v2_out["substitute_kos"] = _step_key.map(_step_sub_map).fillna("")
                steps_v2_out["canonical_ko"] = _step_key.map(_step_canonical_map).fillna("")
                steps_v2_out["has_substitute"] = steps_v2_out["substitute_kos"] != ""
                # Reactions for substitute KOs
                # A substitute is by definition NOT a member KO of the module,
                # so module_reaction_map has no entry for it and
                # reactions_for_module_bestpath returns "" whenever a
                # substitute wins the path. The correct reactions are already
                # in subs_df["shared_reactions"] that shared reaction is the
                # whole reason the KO qualified as a substitute. Fill from
                # there for any winning substitute.
                _sub_rxn_map = {}
                if "shared_reactions" in subs_df.columns:
                    for _, r in subs_df.iterrows():
                        k = (f"{r['module']}|{r['step']}", str(r["substitute_ko"]).upper())
                        rxns = str(r.get("shared_reactions", "") or "")
                        if rxns:
                            _sub_rxn_map.setdefault(k, set()).update(rxns.split(","))

                def _fill_sub_reactions(idx):
                    #Union existing reactions with the shared reactions of any substitute that appears in this step's best path.
                    row = steps_v2_out.loc[idx]
                    existing = str(row.get("best_path_reactions", "") or "")
                    have = {r for r in existing.split(",") if r}
                    subs = [s for s in str(row["substitute_kos"]).split(",") if s]
                    if not subs:
                        return existing
                    bp = {k.strip().upper()
                          for k in str(row.get("best_path_kos", "") or "").split(",") if k.strip()}
                    key_prefix = f"{row['module']}|{row['step']}"
                    for s in subs:
                        if s.upper() in bp:
                            have |= _sub_rxn_map.get((key_prefix, s.upper()), set())
                    return ",".join(sorted(r for r in have if r))

                steps_v2_out["best_path_reactions"] = [
                    _fill_sub_reactions(i) for i in steps_v2_out.index
                ]

                # did a substitute actually win the reported best path?
                steps_v2_out["substitute_won_best_path"] = [
                    bool(subs) and any(
                        k in str(bp).upper().split(",") for k in subs.split(",")
                    )
                    for subs, bp in zip(
                        steps_v2_out["substitute_kos"],
                        steps_v2_out["best_path_kos"].fillna(""),
                    )
                ]

                # Pass-1 vs pass-2 comparison
                # Merge the real pass-1 numbers from steps_df (which used Dk_Neighbor and the unaugmented equations) and diff them.
                # Key on module+step: step ids are "1","2" for single-line and "1.1","1.2" for multiline, so the pair is unique without relying on the multiline dtype.
                _pre = steps_df.copy()
                _pre["_k"] = (
                    _pre["module"].astype(str) + "|" + _pre["step"].astype(str)
                )

                def _pre_lookup(col, default=np.nan):
                    if col not in _pre.columns:
                        logging.warning(
                            "[SUBST] pass-1 steps_df has no '%s'; "
                            "before/after comparison for it will be blank.", col
                        )
                        return {}
                    return dict(zip(_pre["_k"], _pre[col]))

                _pre_p   = _pre_lookup("p_after")
                _pre_bpk = _pre_lookup("best_path_kos")
                _pre_bps = _pre_lookup("best_path_score")
                _pre_eq  = _pre_lookup("equation")

                steps_v2_out["step_probability_pre_sub"] = _step_key.map(_pre_p)
                steps_v2_out["step_probability_delta"] = (
                    steps_v2_out["p_after"] - steps_v2_out["step_probability_pre_sub"]
                )
                steps_v2_out["step_probability_improved"] = (
                    steps_v2_out["step_probability_delta"] > 1e-12
                )

                steps_v2_out["equation_pre_sub"] = _step_key.map(_pre_eq)
                steps_v2_out["best_path_kos_pre_sub"] = _step_key.map(_pre_bpk)
                steps_v2_out["best_path_score_pre_sub"] = _step_key.map(_pre_bps)
                steps_v2_out["best_path_changed"] = [
                    set(str(a or "").upper().split(",")) != set(str(b or "").upper().split(","))
                    for a, b in zip(
                        steps_v2_out["best_path_kos_pre_sub"],
                        steps_v2_out["best_path_kos"],
                    )
                ]

                # did substitution actually help? summary
                _sub_steps = steps_v2_out[steps_v2_out["has_substitute"]]
                _n_sub     = len(_sub_steps)
                _n_improved = int(_sub_steps["step_probability_improved"].sum())
                _n_pathchg  = int(_sub_steps["best_path_changed"].sum())
                _n_subwon   = int(_sub_steps["substitute_won_best_path"].sum())
                _mean_delta = (
                    float(_sub_steps["step_probability_delta"].mean())
                    if _n_sub else 0.0
                )
                _max_delta = (
                    float(_sub_steps["step_probability_delta"].max())
                    if _n_sub else 0.0
                )
                logging.info("--- Substitution summary ---")
                logging.info("  Steps with at least one substitute candidate : %d", _n_sub)
                logging.info("  Steps where substitution improved probability : %d", _n_improved)
                logging.info("  Steps where best-path KOs changed             : %d", _n_pathchg)
                logging.info("  Steps where substitute won the best path      : %d", _n_subwon)
                logging.info("  Mean probability gain (substituted steps)     : %+.6f", _mean_delta)
                logging.info("  Max probability gain  (substituted steps)     : %+.6f", _max_delta)
                logging.info("----------------------------")

                dk_v2_out = dk_v2_out.rename(columns={
                    "KO_Neighbors_all":          "all_associates",
                    "KO_Neighbor_All_Count":     "all_associate_count",
                    "KO_Neighbors_dropped":      "dropped_associates",
                    "KO_Neighbor_Dropped_Count": "dropped_associate_count",
                    "KO_Neighbor_Drop_Reasons":  "drop_reasons",
                })


        # ── Summary of everything this stage wrote ──
        _expected = [
            f"{output_prefix}_substitutes.csv",
            f"{output_prefix}_BLIMMP_substituted_dk.csv",
            f"{output_prefix}_BLIMMP_substituted_module_steps.csv",
            f"{output_prefix}_module_probabilities_sub_v2.csv",
        ]
        if verbose:
            _expected += [
                f"{output_prefix}_substitute_details.csv",
                f"{output_prefix}_substitute_augmentation_audit.csv",
                f"{output_prefix}_substitute_v2_check.csv",
                f"{output_prefix}_steps_v2.csv",
                f"{output_prefix}_modules_v2.csv",
            ]

        return {
            "subs_df": subs_df,
            "sub_details": sub_details,
            "neighbor_map_aug": neighbor_map_aug,
            "dk_full": dk_full,
            "dk_update_calculations": dk_update_calculations,
            "steps_v2": steps_v2,
            "modules_v2": modules_v2,            
            "module_best_v2": module_best_v2,     
            "module_eq_aug": module_eq_aug,
            "dk_v2_out": dk_v2_out,
            "steps_v2_out": steps_v2_out,   
        }


class FileWriters:
    @staticmethod
    def ensure_dir(path):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)


    @staticmethod
    def write_csv_outputs(df_ko: pd.DataFrame,steps_df: pd.DataFrame, modules_df: pd.DataFrame,output_prefix: str,*,basename = "blimmp", extra_ko_cols=None, extra_step_cols=None):
        """
        Write 3 CSVs:
        1) KO CSV
        2) Steps CSV
        3) Modules CSV

        extra_ko_cols / extra_step_cols: optional lists of column names to
        append to the KO / steps CSVs after the standard schema (used by the
        substitution stage to tack on is_substitute, has_substitute etc.
        without forking the schema).
        Pass steps_df=None or modules_df=None to skip writing those tables.
        """
        ko_path = f"{output_prefix}_{basename}_dk.csv"
        FileWriters.ensure_dir(ko_path)
        RENAME_KO_TABLE = {
            "KO id":                          "ko_id",
            "Dk_Neighbor":                    "ko_probability",
            "Dk":                             "ko_probability_prior",
            "hmm_len":                        "hmm_len",
            "target name":                    "orf_name",
            "qlen":                           "orf_len",
            "E-value":                        "evalue",
            "score":                          "hmm_score",
            "i_Evalue":                       "domain_evalue",
            "i_score":                        "domain_score",
            "hmm from":                       "hmm_from",
            "hmm to":                         "hmm_to",
            "ali from":                       "ali_from",
            "ali to":                         "ali_to",
            "overlapgroup_winner":            "orf_best_match",
            "overlapgroup_winner_score":      "orf_best_match_score",
            "overlap_relative_position_confidence": "locus_competition_score",
            "hit_conf":                       "detection_confidence",
            "kofam_score_threshold":          "kofam_threshold",
            "kofam_score_type":               "kofam_score_type",
            "is_outcompeted":                 "outcompeted_at_locus",
            "flag_is_below_kofam_threshold":  "below_kofam_threshold",
            "flag_is_dubious":                "is_dubious",
            "KO_Neighbors":                   "influencing_associates",
            "KO_Neighbor_Count":              "influencing_associate_count",
            "Modules":                        "participating_modules",
            "count":                          "ko_count",
            "KO_freq":                        "prior_frequency",
        }

        df_ko = df_ko.rename(columns=RENAME_KO_TABLE)
        # Build the ordered column list from the rename table, but only include
        # columns that actually landed in df_ko (some HMM hits may be missing
        # optional fields like strand or i_score if the file is sparse).
        keep_ko_cols = [c for c in RENAME_KO_TABLE.values() if c in df_ko.columns]
        if extra_ko_cols:
            keep_ko_cols = keep_ko_cols + [
                c for c in extra_ko_cols if c in df_ko.columns
            ]
        df_ko = df_ko[keep_ko_cols]

        df_ko["ko_count"] = pd.to_numeric(df_ko["ko_count"], errors="coerce")
        df_ko["ko_count"] = (df_ko["ko_count"].replace([np.inf, -np.inf], np.nan).fillna(0).astype(int))
        df_ko["influencing_associate_count"] = df_ko["influencing_associate_count"].astype(int)
        df_ko["ko_probability"] = df_ko["ko_probability"].round(3)
        df_ko["ko_probability_prior"] = df_ko["ko_probability_prior"].round(3)
        df_ko["prior_frequency"] = df_ko["prior_frequency"].round(3)
        df_ko.to_csv(ko_path, index=False)

        # allow steps table to be skipped (pass steps_df=None)
        steps_csv_path = None
        if steps_df is not None:
            steps_csv_path = f"{output_prefix}_{basename}_module_steps.csv"
            FileWriters.ensure_dir(steps_csv_path)
            rename_step_prob = {
                "module": "module",
                "module_description": "module_description",
                "step": "step_id",
                "p_after": "step_probability", #
                "equation": "equation",
                "best_path_kos" : "best_path_kos",
                "best_path_score": "best_path_ko_confidence",
                "best_path_reactions": "best_path_reactions"
            }

            steps_df = steps_df.rename(columns=rename_step_prob)
            keep_step_cols = [c for c in rename_step_prob.values() if c in steps_df.columns]
            if extra_step_cols:
                keep_step_cols = keep_step_cols + [
                    c for c in extra_step_cols if c in steps_df.columns
                ]
            steps_df = steps_df[keep_step_cols]
            steps_df["step_probability"] = steps_df["step_probability"].round(3)
            steps_df["best_path_ko_confidence"] = steps_df["best_path_ko_confidence"].round(3)
            steps_df.to_csv(steps_csv_path, index=False)

        # allow modules table to be skipped (pass modules_df=None)
        modules_csv_path = None
        if modules_df is not None:
            modules_csv_path = f"{output_prefix}_{basename}_module_probabilities.csv"
            FileWriters.ensure_dir(modules_csv_path)

            #Rename column names

            rename_module_prob = {
                "module": "module",
                "module_description": "module_description", 
                "module_probability_after": "module_confidence",
                "module_probability_before": "module_confidence_prior",
                "module_best_path_kos": "best_path",
                "module_best_path_reactions": "best_path_reactions",
                "module_frequency": "module_frequency",
                "n_steps": "num_steps",
                "E_before": "num_steps_present_raw",
                "E_after": "num_steps_present",
                "posterior_mean": "posterior_mean",
                "posterior_variance": "posterior_variance",
                "beta_threshold": "completeness_threshold",
                "ci_low": "credible_interval_low",
                "ci_high": "credible_interval_high",
                "hdi_low": "hdi_low",
                "hdi_high": "hdi_high",
                "effect_size": "effect_size",
                "best_path_warnings": "best_path_warnings",
                "best_path_flagged_kos": "best_path_flagged_kos",
                "best_path_severe_kos": "best_path_severe_kos",
                "module_status": "module_status",
            }
            modules_df = modules_df.rename(columns=rename_module_prob)
            keep_cols = list(rename_module_prob.values())
            modules_df = modules_df[keep_cols]
            modules_df["num_steps"] = modules_df["num_steps"].astype(int)
            modules_df["num_steps_present_raw"] = modules_df["num_steps_present_raw"].round(2)
            modules_df["num_steps_present"] = modules_df["num_steps_present"].round(2)

            modules_df["module_frequency"] = modules_df["module_frequency"].round(3)
            modules_df["module_confidence_prior"] = modules_df["module_confidence_prior"].round(3)
            modules_df["module_confidence"] = modules_df["module_confidence"].round(3)
            modules_df["num_steps_present"] = modules_df["num_steps_present"].round(1)
            modules_df["posterior_mean"] = modules_df["posterior_mean"].round(4)
            modules_df["posterior_variance"] = modules_df["posterior_variance"].round(6)
            modules_df["credible_interval_low"] = modules_df["credible_interval_low"].round(4)
            modules_df["credible_interval_high"] = modules_df["credible_interval_high"].round(4)
            modules_df["hdi_low"] = modules_df["hdi_low"].round(4)
            modules_df["hdi_high"] = modules_df["hdi_high"].round(4)
            modules_df["effect_size"] = modules_df["effect_size"].round(3)
            modules_df.to_csv(modules_csv_path, index=False)

        print(f"[Done.] KO-level file written to {ko_path}")
        if steps_csv_path:
            print(f"[Done.] Step-level file written to {steps_csv_path}")
        if modules_csv_path:
            print(f"[Done.] Module-level file written to {modules_csv_path}")

        return {
        "ko_csv": ko_path,
        "steps_csv": steps_csv_path,
        "modules_csv": modules_csv_path}
    
    
    @staticmethod
    def write_module_json(
        df: pd.DataFrame,
        modules_df: pd.DataFrame,
        steps_df: pd.DataFrame,   
        output_prefix: str,
        *,
        basename: str = "blimmp",
        module_json_dir: str| None = None,
        module_reaction_map
    ):

        def _sanitize_for_json(x):
            # None / pandas NA
            if x is None or x is pd.NA:
                return None

            # numpy scalars
            if isinstance(x, (np.integer,)):
                return int(x)
            if isinstance(x, (np.bool_,)):
                return bool(x)
            if isinstance(x, (np.floating,)):
                x = float(x)
                return None if not math.isfinite(x) else x

            # python floats
            if isinstance(x, float):
                return None if not math.isfinite(x) else x

            # containers
            if isinstance(x, dict):
                return {str(k): _sanitize_for_json(v) for k, v in x.items()}
            if isinstance(x, (list, tuple)):
                return [_sanitize_for_json(v) for v in x]

            if isinstance(x, (str, int, bool)):
                return x

            return str(x)

        

        def _json_default(o):
            if o is None:
                return None
            if o is pd.NA:
                return None
            if isinstance(o, (np.integer,)):
                return int(o)
            if isinstance(o, (np.floating,)):
                x = float(o)
                return None if not np.isfinite(x) else x
            if isinstance(o, (np.bool_,)):
                return bool(o)
            return str(o)

        json_path = f"{output_prefix}_{basename}_modules.json"
        os.makedirs(os.path.dirname(json_path) or ".", exist_ok=True)

        want_cols = ["module_equation", "module_probability_before", "module_probability_after", "posterior_mean", "posterior_variance", "beta_threshold", "ci_low", "ci_high"]
        if "module_best_path_kos" in modules_df.columns:
            want_cols.append("module_best_path_kos")
        if "was_substituted" in modules_df.columns:
            want_cols.append("was_substituted")
        if "substitution_delta" in modules_df.columns:
            want_cols.append("substitution_delta")

        mod_prob = modules_df.set_index("module")[want_cols].to_dict(orient="index")



        df_nodes = df.copy()
        df_nodes["KO id"] = df_nodes["KO id"].astype(str).str.strip().str.upper()

        df_nodes = df_nodes[df_nodes["Modules"].notna()].copy()  # guard NaNs

        df_nodes["modules_present"] = (
            df_nodes["Modules"]
            .astype(str)
            .str.split(",")
            .apply(lambda xs: [s.strip() for s in xs if s and s.strip()])
        )

        df_nodes = df_nodes.assign(module_list=df_nodes["Modules"].astype(str).str.split(","))
        df_nodes = df_nodes.explode("module_list")
        df_nodes["module_list"] = df_nodes["module_list"].astype(str).str.strip()
        node_fields = [
            "KO id","target name","E-value","score",
            "overlapgroup_winner","overlapgroup_winner_score",
            "overlapgroup_winner_hit_conf","overlap_relative_position_confidence",
            "kofam_score_threshold","hit_conf","flag_is_dubious","is_outcompeted",
            "flag_is_below_kofam_threshold","KO_freq","Dk","Dk_Neighbor",
            "KO_Neighbors","KO_Neighbor_Count", "buddy_stats","modules_present",
        ]
        present = [c for c in node_fields if c in df_nodes.columns]

        # Unfiltered, KO-indexed lookup — built AFTER modules_present so its
        # schema matches df_nodes/`present` exactly (same columns, same
        # dtypes), and keeping "KO id" as a real column (drop=False) rather
        # than consuming it into the index. A KO can be relevant to a
        # module's graph without being a native member of that module's own
        # equation — e.g. it outcompetes a canonical enzyme at a shared
        # locus (overlapgroup_winner), or it was accepted as a functional
        # substitute (Eq. 27) for one. Neither case requires the KO to have
        # any Modules membership of its own, so such KOs must be reachable
        # here even though df_nodes above has already dropped them.
        ko_lookup_base = df.copy()
        ko_lookup_base["KO id"] = ko_lookup_base["KO id"].astype(str).str.strip().str.upper()
        if "modules_present" not in ko_lookup_base.columns:
            ko_lookup_base["modules_present"] = [[] for _ in range(len(ko_lookup_base))]
        ko_lookup = (
            ko_lookup_base.drop_duplicates(subset="KO id", keep="first")
                          .set_index("KO id", drop=False)
        )

        aggregated = {}


        # 
        steps_pack_by_module: dict[str, dict] = {}

        if steps_df is not None and not steps_df.empty:
            rename_map = {
                "module_id": "module",
                "step_no": "step",
                "step_equation": "equation",
                "step_prob_before": "p_before",
                "step_prob_after": "p_after",
                "best_path_kos": "best_path_kos",
                "best_path_score": "best_path_score",
                "best_path_reactions": "best_path_reactions",

            }
            steps_norm = steps_df.rename(
                columns={k: v for k, v in rename_map.items() if k in steps_df.columns}
            ).copy()

            required = {"module", "step", "equation", "p_before", "p_after"}
            missing = sorted(required - set(steps_norm.columns))
            if missing:
                raise ValueError(f"steps_df missing required columns: {missing}. Have: {list(steps_norm.columns)}")

            steps_norm["module"]   = steps_norm["module"].astype(str).str.strip()
            steps_norm["step"] = steps_norm["step"].astype(str).str.strip()
            steps_norm["step"] = pd.to_numeric(steps_norm["step"], errors="coerce")

            steps_norm["equation"] = steps_norm["equation"].astype(str)
            steps_norm["p_before"] = pd.to_numeric(steps_norm["p_before"], errors="coerce")
            steps_norm["p_after"]  = pd.to_numeric(steps_norm["p_after"], errors="coerce")

            # optional cols
            if "best_path_kos" in steps_norm.columns:
                steps_norm["best_path_kos"] = steps_norm["best_path_kos"].astype(str).str.strip()

            if "best_path_score" in steps_norm.columns:
                steps_norm["best_path_score"] = pd.to_numeric(steps_norm["best_path_score"], errors="coerce")
            if "best_path_reactions" in steps_norm.columns:
                steps_norm["best_path_reactions"] = steps_norm["best_path_reactions"].astype(str).fillna("").str.strip()


            def _split_step_id(x: float) -> tuple[int, float]:
                if not np.isfinite(x):
                    return (1, float("nan"))
                line = int(np.floor(x))
                inner = round((x - line) * 10, 6)  # .1 -> 1.0, .2 -> 2.0
                return (line, inner)

            for mod_id, sub in steps_norm.groupby("module", sort=False):
                sub = sub.dropna(subset=["step"]).sort_values("step")

                # detect multiline: any non-integer step id (like 1.1) OR multiple integer "lines"
                has_decimal = bool(((sub["step"] % 1) != 0).any())

                if has_decimal:
                   
                    sub = sub.copy()
                    sub[["line", "inner_step"]] = sub["step"].apply(lambda v: pd.Series(_split_step_id(float(v))))
                    lines_out = []
                    for line_no, gline in sub.groupby("line", sort=True):
                        gline = gline.sort_values("inner_step")
                        steps_list = []
                        for _, row in gline.iterrows():
                            step_dict = {
                                "step": float(row["inner_step"]),     # inner step within the line (e.g., 1.0)
                                "equation": row["equation"],
                                "p_before": float(row["p_before"]) if np.isfinite(row["p_before"]) else None,
                                "p_after":  float(row["p_after"]) if np.isfinite(row["p_after"]) else None,
                                
                            }
                            if "best_path_kos" in gline.columns:
                                step_dict["best_path_kos"] = row.get("best_path_kos", None)

                            if "has_substitute" in gline.columns:
                                step_dict["has_substitute"] = bool(row.get("has_substitute", False))
                            if "canonical_ko" in gline.columns:
                                step_dict["canonical_ko"] = row.get("canonical_ko", "") or ""
                            if "substitute_kos" in gline.columns:
                                step_dict["substitute_kos"] = row.get("substitute_kos", "") or ""

                            if "best_path_reactions" in gline.columns:
                                step_dict["best_path_reactions"] = row.get("best_path_reactions", "") or ""

                            
                            if "best_path_score" in gline.columns:
                                v = row.get("best_path_score", None)
                                step_dict["best_path_score"] = float(v) if (v is not None and np.isfinite(v)) else None
                            steps_list.append(step_dict)

                        steps_inline = "; ".join([f"step {int(s['step']) if s['step']%1==0 else s['step']}: {s['equation']}" for s in steps_list])

                        lines_out.append({
                            "line": int(line_no),
                            "steps": steps_list,
                            "steps_inline": steps_inline
                        })

                    steps_pack_by_module[mod_id] = {"lines": lines_out}

                else:
                    steps_list = []
                    for _, row in sub.iterrows():
                        step_dict = {
                            "step": float(row["step"]),
                            "equation": row["equation"],
                            "p_before": float(row["p_before"]) if np.isfinite(row["p_before"]) else None,
                            "p_after":  float(row["p_after"]) if np.isfinite(row["p_after"]) else None,
                        }
                        if "best_path_kos" in sub.columns:
                            step_dict["best_path_kos"] = row.get("best_path_kos", None)

                        if "has_substitute" in sub.columns:
                            step_dict["has_substitute"] = bool(row.get("has_substitute", False))
                        if "canonical_ko" in sub.columns:
                            step_dict["canonical_ko"] = row.get("canonical_ko", "") or ""
                        if "substitute_kos" in sub.columns:
                            step_dict["substitute_kos"] = row.get("substitute_kos", "") or ""

                        if "best_path_reactions" in sub.columns:
                            step_dict["best_path_reactions"] = row.get("best_path_reactions", "") or ""

                        if "best_path_score" in sub.columns:
                            v = row.get("best_path_score", None)
                            step_dict["best_path_score"] = float(v) if (v is not None and np.isfinite(v)) else None
                        steps_list.append(step_dict)

                    steps_pack_by_module[mod_id] = {"steps": steps_list}
        

        nodes_by_module: dict[str, pd.DataFrame] = {}
        if "module_list" in df_nodes.columns:
            for mod_id, grp in df_nodes.groupby("module_list", sort=False):
                mod_id = str(mod_id).strip()
                if mod_id:
                    nodes_by_module[mod_id] = grp

        # Write all modules that have module-level metadata OR step info
        all_mod_ids = set(mod_prob.keys())
        all_mod_ids |= set(steps_pack_by_module.keys())  

        for mod_id in sorted(all_mod_ids):
            if not mod_id or mod_id not in mod_prob:
                continue

            meta = mod_prob.get(mod_id, {
                "module_equation": "",
                "module_probability_before": None,
                "module_probability_after": None,
                "module_best_path_kos": None,
            })

            grp = nodes_by_module.get(mod_id, None)

            # base nodes from df (limited to columns actually present)
            if grp is not None and not grp.empty:
                nodes = grp[present].to_dict(orient="records")
            else:
                nodes = []  # will become START/SINK only

            _base_defaults = {
                "KO id": None,
                "target name": "NA",
                "E-value": None,
                "score": None,
                "overlapgroup_winner": "NA",
                "overlapgroup_winner_score": None,
                "overlapgroup_winner_hit_conf": None,
                "overlap_relative_position_confidence": None,
                "kofam_score_threshold": None,
                "hit_conf": None,
                "flag_is_dubious": False,
                "is_outcompeted": False,
                "flag_is_below_kofam_threshold": False,
                "KO_freq": None,
                "Dk": None,
                "Dk_Neighbor": None,
                "KO_Neighbors": "NA",
                "KO_Neighbor_Count": 0,
                "buddy_stats": None,
            }

            _start = _base_defaults.copy()
            _start.update({"KO id": "START", "Dk": 1.0, "Dk_Neighbor": 1.0})
            start_row = {k: v for k, v in _start.items() if k in present}

            _sink = _base_defaults.copy()
            _sink.update({"KO id": "SINK", "Dk": 1.0, "Dk_Neighbor": 1.0})
            sink_row = {k: v for k, v in _sink.items() if k in present}

            # Substitute KOs that won a best path in this module are not members
            # of the module's own equation (Eq. 27 condition 3), so nodes_by_module
            # never includes them — they need to be pulled in by id from
            # ko_lookup (unfiltered by Modules membership) and given the same
            # field set as every other node.
            pack_for_subs = steps_pack_by_module.get(mod_id, {})

            def _iter_step_dicts(pack):
                for step in pack.get("steps", []):
                    yield step
                for line in pack.get("lines", []):
                    for step in line.get("steps", []):
                        yield step

            canonical_for_sub = {}   # substitute_ko -> canonical_ko it replaced
            step_for_sub = {}        # substitute_ko -> step id it was placed at
            for step in _iter_step_dicts(pack_for_subs):
                sub_field = step.get("substitute_kos", "") or ""
                can_field = step.get("canonical_ko", "") or ""
                for sub_ko in [k for k in sub_field.split(",") if k]:
                    canonical_for_sub[sub_ko] = can_field
                    step_for_sub[sub_ko] = step.get("step")

            existing_ko_ids = {n.get("KO id") for n in nodes}
            for sub_ko, canonical_ko in canonical_for_sub.items():
                if sub_ko in existing_ko_ids:
                    continue
                if sub_ko not in ko_lookup.index:
                    continue
                # Same field set as every canonical node — no stripped-down stub
                sub_node = ko_lookup.loc[sub_ko, present].to_dict()
                sub_node["KO id"] = sub_ko
                sub_node["is_substitute_node"] = True
                sub_node["substitutes_for"] = canonical_ko
                sub_node["substitute_step"] = step_for_sub.get(sub_ko)
                nodes.append(sub_node)
                existing_ko_ids.add(sub_ko)

            # Outcompeting KOs (overlapgroup_winner) are, by the same logic,
            # not necessarily members of this module's own equation — a
            # winner can belong to a different module entirely, or to none
            # (e.g. K00812 outcompeting K00832 in M00034 while K00812 itself
            # has no KEGG module membership). Every node already carries its
            # own overlapgroup_winner id/score/hit_conf inline, but the
            # winner itself should also get a full node so it can be
            # inspected (hovered, expanded) the same way a substitute is,
            # rather than existing only as a name string on the loser's row.
            for n in list(nodes):
                winner_ko = str(n.get("overlapgroup_winner", "") or "").strip().upper()
                if not winner_ko or winner_ko in ("NA", "NAN"):
                    continue
                if winner_ko == str(n.get("KO id", "")).strip().upper():
                    continue  # KO wasn't outcompeted; it's its own winner
                if winner_ko in existing_ko_ids:
                    continue
                if winner_ko not in ko_lookup.index:
                    continue
                winner_node = ko_lookup.loc[winner_ko, present].to_dict()
                winner_node["KO id"] = winner_ko
                winner_node["is_outcompeting_node"] = True
                winner_node["outcompetes"] = n.get("KO id")
                nodes.append(winner_node)
                existing_ko_ids.add(winner_ko)

            nodes.insert(0, start_row)
            nodes.append(sink_row)

            entry = {
                "module_equation": meta.get("module_equation", ""),
                "module_probability_before": meta.get("module_probability_before", None),
                "module_probability_after": meta.get("module_probability_after", None),
                "best_path": meta.get("module_best_path_kos", None),  
                "nodes": nodes,
            }

            pack = steps_pack_by_module.get(mod_id) 
            if pack:
                entry.update(pack)

            aggregated[mod_id] = entry



        missing_modules = sorted(set(modules_df["module"]) - set(aggregated.keys()))
        if missing_modules:
            print(f"{len(missing_modules)} modules omitted (not present in df or no nodes):")
            print(", ".join(missing_modules[:50]))

        aggregated = _sanitize_for_json(aggregated)

        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(aggregated, f, indent=2, ensure_ascii=False, allow_nan=False)


        print(f"[Done.] Wrote {json_path} with {len(aggregated)} modules.")
        return json_path

    
    def merge_substitution_results(modules_df_pass1, sub_result, module_descriptions=None):
        #Return a single module confidence table: For modules that had a substitute AND the substituted confidence is higher, use the pass-2 confidence. All other modules keep their pass-1 values.
        merged = modules_df_pass1.copy()
        merged["was_substituted"] = False
        merged["substitution_delta"] = 0.0
        merged["confidence_source"] = "pass1"

        modules_v2 = sub_result.get("modules_v2")
        if modules_v2 is None or modules_v2.empty:
            return merged

        v2_idx = modules_v2.set_index("module")

        for i, row in merged.iterrows():
            mid = row["module"]
            if mid not in v2_idx.index:
                continue
            v2_row = v2_idx.loc[mid]
            conf_v2 = float(v2_row["module_probability_after"])
            conf_v1 = float(row["module_probability_after"])
            delta = conf_v2 - conf_v1

            # Always reflect the accepted substitution's best path, regardless
            # of whether it raised the module's overall confidence enough to
            # matter — a substitution can be a real, correct enzyme-level
            # change even when other, unrelated failing steps keep the
            # module's overall call unchanged.
            if "module_best_path_kos" in v2_row.index:
                merged.at[i, "module_best_path_kos"] = v2_row["module_best_path_kos"]

            if delta > 1e-100:   # if substitution helped at all
                merged.at[i, "module_probability_after"] = conf_v2
                merged.at[i, "E_after"]           = float(v2_row["E_after"])
                merged.at[i, "posterior_mean"]    = float(v2_row["posterior_mean"])
                merged.at[i, "posterior_variance"]= float(v2_row["posterior_variance"])
                merged.at[i, "ci_low"]            = float(v2_row["ci_low"])
                merged.at[i, "ci_high"]           = float(v2_row["ci_high"])
                merged.at[i, "hdi_low"]           = float(v2_row["hdi_low"])
                merged.at[i, "hdi_high"]          = float(v2_row["hdi_high"])
                merged.at[i, "effect_size"]       = float(v2_row["effect_size"])
                merged.at[i, "was_substituted"]   = True
                merged.at[i, "substitution_delta"]= round(delta, 4)
                merged.at[i, "confidence_source"] = "pass2"

        return merged


## Main Callers

class BlimmpPipeline:
    def __init__(self, cfg: RunConfig, paths: Paths):
        self.cfg = cfg
        self.paths = paths

    def run(self):

        #File loaders
        sample_name = os.path.basename(self.cfg.input_file).split('.')[0]
        logging.info(f"Processing sample {sample_name} with sigma={self.cfg.sigma}")
        if not (0.0 <= self.cfg.sigma <= 1.0):
            raise ValueError("--sigma must be between 0 and 1")
        
        # taxonomy-driven paths
        counts_tsv, onehop_json, twohop_json, all_neighbor_json, cooc_json, tag = File_Helpers.lineage_paths(self.cfg.taxonomy, self.paths)
        logging.info(f"Taxonomic level chosen: {self.cfg.taxonomy}")

        ko_occ = File_Helpers.read_ko_occurrence(str(counts_tsv))
        if self.cfg.verbose:
            logging.debug("KO occurrence head:\n%s", ko_occ.head())

        #Ignore, warming up JIT
        _ = _assign_groups_numba(np.array([0.,1.]), np.array([1.,2.]), 0.6)

        hmm_hits = HMMParsers.process_domtblout(self.cfg.input_file)

        hmm_groups = Overlap.assign_overlap_groups(hmm_hits)
        # KO-level de-dup (best score per KO)
        hmm_groups_dedup = (hmm_groups.sort_values('score', ascending=False).drop_duplicates(subset='KO id', keep='first').reset_index(drop=True))

        # Reporting values
        counts = hmm_groups_dedup[["overlap_group", "KO id"]].agg(pd.Series.nunique)
        logging.info("Resolving overlapping HMM hits — grouping annotations at the same genomic locus...")
        logging.info("Found %d overlapping locus groups covering %d unique KO annotations.",counts['overlap_group'], counts['KO id'])



        #Loading the kofamdb file
        positionscored = PositionScores.winner_info_and_flags(hmm_groups_dedup, self.paths.kofam_ko_list_path)
        if self.cfg.verbose:
            logging.debug("positionscored head:\n%s", positionscored.head())


        # Build KO universe straight from my module JSONs
        ko_to_modules_str = File_Helpers.modules_to_kos(self.paths.module_json_dir)
        annotated_kos = set(positionscored['KO id'].astype(str).str.strip().str.upper())
        module_kos    = set(ko_to_modules_str.keys())
        ko_universe   = pd.DataFrame({'KO id': sorted(module_kos | annotated_kos)})

        # Left-join observed annotations; missing KOs get defaults
        positionscored_full = ko_universe.merge(positionscored, on='KO id', how='left')

        # Defaults for KOs not present in HMM hits
        positionscored_full['hit_conf'] = positionscored_full['hit_conf'].fillna(0.0)
        positionscored_full['E-value']  = positionscored_full['E-value'].fillna(1000.0)

        for col in ('score', 'overlap_group'):
            if col in positionscored_full:
                positionscored_full[col] = positionscored_full[col].fillna(0)

        logging.debug("positionscored_full KO sample: %s", positionscored_full['KO id'].head(10).tolist())
        logging.debug("ko_occ KO sample: %s", ko_occ['KO id'].head(10).tolist())
        logging.debug("positionscored_full shape: %s", positionscored_full.shape)

        #Calculating Dk from positionscored
        dk_calculations = CalculateKOProbabilities.calculate_dk_per_ko(positionscored_full, ko_occ, verbose=self.cfg.verbose)
        logging.info("KO probability table built: %d KOs (%d in modules, %d detected by HMM).", len(dk_calculations), len(module_kos), len(annotated_kos))
        if self.cfg.verbose:
            logging.debug("Dk calculations head:\n%s", dk_calculations.head())

        #In the future you can choose two/one hop
        neighbor_map, ko_counts = NeighborCalculations.make_neighbor_dictionary(all_neighbor_json)

        #Update the Dk calculations based on the neighbors
        dk_update_calculations, used_neighbors = CalculateKOProbabilities.dk_neighbor_update(dk_calculations,neighbor_map, ko_counts, alpha=0.6, return_used=True,verbose=self.cfg.verbose)
        print("Pass 1 complete — largest neighborhood-driven probability shift: %.6f",(dk_update_calculations["Dk_Neighbor"] - dk_update_calculations["Dk"]).abs().max())
        # attach neighbor lists/counts per KO
        dk_update_calculations["KO_Neighbors"] = dk_update_calculations["KO id"].map(lambda k: ",".join(used_neighbors.get(k, [])))
        dk_update_calculations["KO_Neighbor_Count"] = dk_update_calculations["KO id"].map(lambda k: len(used_neighbors.get(k, []))).fillna(0).astype(int)
        dk_full = dk_update_calculations.copy()

        # Add module details
        ko_to_modules_str = File_Helpers.modules_to_kos(self.paths.module_json_dir)
        dk_update_calculations["Modules"] = dk_update_calculations["KO id"].map(ko_to_modules_str)

        # Add this print after dk_full is set, to confirm new KOs are there
        new_kos_in_full = set(dk_full["KO id"]) - module_kos
        logging.info("%d detected KOs are not members of any KEGG module — included in dk_full for substitute lookup but excluded from module scoring.",len(new_kos_in_full))

        dk_update_calculations = dk_update_calculations[dk_update_calculations["Modules"].notna()].copy()
        logging.info("Module-member KO table: %d KOs retained for step probability evaluation.", dk_update_calculations["KO id"].nunique())

        module_eq = File_Helpers.load_module_eq(self.paths.module_eq_json)
        module_freq_dict = File_Helpers.load_module_freq(self.paths.module_frequencies)

        single_only = {m:e for m,e in module_eq.items() if e.get("steps")}
        multi_only  = {m:e for m,e in module_eq.items() if e.get("lines")}


        steps_df_single = CalculateModuleProbabilities.evaluate_step_probabilities(single_only,dk_update_calculations,before_col="Dk",after_col="Dk_Neighbor",verbose=self.cfg.verbose)
        #steps_df["multiline"] = False


        modules_df_single = CalculateModuleProbabilities.calculate_module_confidence(
        steps_df=steps_df_single,
        module_dict=single_only,
        genome_completeness=self.cfg.sigma,  
        module_frequencies=module_freq_dict,
        verbose=self.cfg.verbose)


        # 2) multiline (new)
        steps_df_multi = CalculateModuleProbabilities.evaluate_multiline_step_probabilities(
            multi_only, dk_update_calculations,
            before_col="Dk", after_col="Dk_Neighbor",
            step_format="path.step",
            verbose=self.cfg.verbose
        )

        modules_df_multi = CalculateModuleProbabilities.calculate_multiline_module_confidence_from_steps(
            steps_df_multi,
            module_dict_multiline=multi_only,
            genome_completeness=self.cfg.sigma,
            module_frequencies=module_freq_dict,
            verbose=self.cfg.verbose
        )


        best_paths = ModuleBestPath(
            module_eq=module_eq,
            ko_df=dk_update_calculations,
            ko_id_col="KO id",
            ko_prob_col="Dk_Neighbor",
            score_col="score",         
            keep_duplicate_kos="max",
        )


        best_steps_df, best_failures_df = best_paths.run_all()

        if self.cfg.verbose:
            best_failures_df.to_csv(
                f"{self.cfg.output_prefix}_best_path_fails.csv", index=False)

        steps_df   = pd.concat([steps_df_single, steps_df_multi], ignore_index=True).sort_values(["module","multiline","step"]).reset_index(drop=True)

        steps_df["step"] = steps_df["step"].astype(str)
        best_steps_df["step"] = best_steps_df["step"].astype(str)

        steps_df = steps_df.merge(
            best_steps_df[["module","multiline","step","best_path_score","best_path_kos"]],
            on=["module","multiline","step"],
            how="left"
        )

        ko_to_reactions = File_Helpers.ko_to_reactions_dict(str(self.paths.ko_reaction_path))

        module_reaction_map = File_Helpers.load_module_reactions(self.paths.module_reaction_dir)
        steps_df["best_path_reactions"] = steps_df.apply(lambda row: File_Helpers.reactions_for_module_bestpath(module_reaction_map, row["module"], row["best_path_kos"]),axis=1)


        module_best_df = best_paths.compute_module_best_paths(steps_df)
        modules_df = pd.concat([modules_df_single, modules_df_multi], ignore_index=True).sort_values(["module"]).reset_index(drop=True)
        modules_df = modules_df.merge(module_best_df, on="module", how="left")
        #print(modules_df.columns)
        modules_df["module_best_path_reactions"] = modules_df.apply(lambda row: File_Helpers.reactions_for_module_bestpath(module_reaction_map, row["module"], row["module_best_path_kos"]),axis=1)
        # Load module descriptions and add to steps_df and modules_df
        module_descriptions = File_Helpers.load_module_descriptions(self.paths.module_descriptions_path)       
        steps_df["module_description"] = steps_df["module"].map(module_descriptions).fillna("")
        modules_df["module_description"] = modules_df["module"].map(module_descriptions).fillna("")

        modules_df = CalculateModuleProbabilities.annotate_module_ko_warnings(modules_df, dk_update_calculations)

        modules_df["module_status"] = modules_df.apply(
            lambda r: CalculateModuleProbabilities.classify_module(
                r["module_probability_after"], r["ci_low"], r["ci_high"],
                r["hdi_low"], r["hdi_high"], r["beta_threshold"],
                r["posterior_mean"], r["n_steps"], r["effect_size"], r.get("best_path_warnings", ""), r.get("best_path_flagged_kos", 0),r.get("best_path_severe_kos", 0), r.get("best_path_total_kos", 0)),axis=1)

        if self.cfg.no_substitutes or self.cfg.verbose:
            FileWriters.write_csv_outputs(dk_update_calculations, steps_df, modules_df, self.cfg.output_prefix, basename="BLIMMP",)
            FileWriters.write_module_json(dk_update_calculations, modules_df, steps_df, self.cfg.output_prefix, basename="BLIMMP",module_json_dir=str(self.paths.module_json_dir),module_reaction_map=module_reaction_map,)

        ######### PASS 2 ########
        # Skipped only when --no-substitutes. Default and --verbose both run it.
        if not self.cfg.no_substitutes:
            substitution_results = SubstitutionPipeline.run(
                steps_df=steps_df,
                dk_full=dk_full,
                dk_calculations=dk_calculations,
                dk_update_calculations=dk_update_calculations,
                neighbor_map=neighbor_map,
                ko_counts=ko_counts,
                ko_to_reactions=ko_to_reactions,
                cooc_json=cooc_json,
                single_only=single_only,
                multi_only=multi_only,
                module_eq=module_eq,
                genome_completeness=self.cfg.sigma,
                module_frequencies=module_freq_dict,
                module_reaction_map=module_reaction_map,
                module_descriptions=module_descriptions,
                modules_df_pre=modules_df,
                output_prefix=self.cfg.output_prefix,
                verbose=self.cfg.verbose,
            )
            subs_df                = substitution_results["subs_df"]
            sub_details             = substitution_results["sub_details"]
            neighbor_map_aug        = substitution_results["neighbor_map_aug"]
            dk_full                 = substitution_results["dk_full"]
            dk_update_calculations  = substitution_results["dk_update_calculations"]
            steps_v2                = substitution_results["steps_v2"]
            modules_v2              = substitution_results["modules_v2"]
            module_best_v2          = substitution_results["module_best_v2"]
            module_eq_aug           = substitution_results["module_eq_aug"]
            dk_v2_out               = substitution_results["dk_v2_out"]
            steps_v2_out            = substitution_results["steps_v2_out"]

            

            final_modules_df = FileWriters.merge_substitution_results(modules_df, substitution_results, module_descriptions)

            # Re-run warning annotation against the POST-SUBSTITUTION best path and KO table.
            final_modules_df = CalculateModuleProbabilities.annotate_module_ko_warnings(final_modules_df, dk_v2_out)

            # Only now classify, so the status string's confidence number
            final_modules_df["module_status"] = final_modules_df.apply(
                lambda r: CalculateModuleProbabilities.classify_module(
                    r["module_probability_after"], r["ci_low"], r["ci_high"],
                    r["hdi_low"], r["hdi_high"], r["beta_threshold"],
                    r["posterior_mean"], r["n_steps"], r["effect_size"],
                    r.get("best_path_warnings", ""), r.get("best_path_flagged_kos", 0),
                    r.get("best_path_severe_kos", 0), r.get("best_path_total_kos", 0)
                ), axis=1)

            FileWriters.write_csv_outputs(dk_v2_out, steps_v2_out, final_modules_df, self.cfg.output_prefix, basename="BLIMMP_substituted", extra_ko_cols=["is_substitute", "substitutes_for", "substituted_associates"], extra_step_cols=["has_substitute", "canonical_ko", "substitute_kos"])
            #FileWriters.write_module_json(dk_update_calculations, final_modules_df, steps_v2_out, self.cfg.output_prefix, basename="BLIMMP_substituted",module_json_dir=str(self.paths.module_json_dir),module_reaction_map=module_reaction_map,)
            FileWriters.write_module_json(dk_v2_out, final_modules_df, steps_v2_out, self.cfg.output_prefix, basename="BLIMMP_substituted",module_json_dir=str(self.paths.module_json_dir),module_reaction_map=module_reaction_map,)





def main():
    p = argparse.ArgumentParser(
    description='BLIMMP: Bayesian Likelihood Inference of Metabolic Module Presence. '
                'Evaluates KEGG module completeness from HMM search results.',
    formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
        'examples:\n'
        '  # default: run substitution pass, write pass-2 outputs only\n'
        '  python BLIMMP/BLIMMP_Scripts/module_detection.py \\\n'
        '    BLIMMP/Examples/example.domtblout \\\n'
        '    -f domtblout --sigma 1.0 --output example_name\n'
        '\n'
        '  # --no-substitutes: skip substitution, write pass-1 outputs only\n'
        '  python BLIMMP/BLIMMP_Scripts/module_detection.py \\\n'
        '    BLIMMP/Examples/example.domtblout \\\n'
        '    -f domtblout --sigma 1.0 --output example_name --no-substitutes\n'
        '\n'
        '  # --verbose: run substitution, write all outputs + intermediates + debug log\n'
        '  python BLIMMP/BLIMMP_Scripts/module_detection.py \\\n'
        '    BLIMMP/Examples/example.domtblout \\\n'
        '    -f domtblout --sigma 1.0 --output example_name --verbose\n'
    ))
    
    p.add_argument('file',help='Path to the HMMER .tblout or .domtblout file')
    p.add_argument('-f', '--format',choices=['domtblout'], required=True,help='Input file format: "domtblout" for --domtblout')
    p.add_argument('-s', '--sigma',type=float, required=True,help='Genome completeness estimate (0.0-1.0). Use 1.0 if unknown or for complete genomes')
    p.add_argument('-t', '--taxonomy',default="bacteria", metavar="NAME",help='Taxonomic group for priors (default: bacteria). Options include phylum names like "cyanobacteriota" or kingdom names like "pseudomonadati"')
    p.add_argument('-o', '--output', required=True, metavar="PREFIX", help='Output prefix for result files (e.g., "results/Genomename_Result")')
    #p.add_argument('-l', '--logfile', action='store_true', help='Write a verbose debug log to <PREFIX>_debug.log')
    #p.add_argument('--write-intermediates', action='store_true', help='Write substitution debug CSVs (substitute_details, ' 'augmentation_audit, substitute_v2_check, steps_v2, modules_v2). ' 'Off by default; the main outputs already contain this data.')
    p.add_argument('-v', '--verbose', action='store_true',help=('Verbose mode: write a debug log to <PREFIX>_debug.log and save intermediate diagnostic CSVs (substitute_details, augmentation_audit, substitute_v2_check).'))
    p.add_argument('--no-substitutes', action='store_true', help=('Skip the substitution pipeline and write only pass-1 output files (KO table, module steps, module probabilities). By default BLIMMP runs the substitution pipeline and writes only the substituted versions. Use -v/--verbose to write both.'))
    
    args = p.parse_args()

    logfile_path = f"{args.output}_debug.log" if args.verbose else None
    logging.getLogger("numba").setLevel(logging.ERROR)

    logger = logging.getLogger()
    logger.handlers.clear()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

        file_handler = logging.FileHandler(logfile_path, mode="w")
        file_handler.setLevel(logging.DEBUG)

        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)

        formatter = logging.Formatter("[%(levelname)s] %(message)s")
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)

        logger.addHandler(file_handler)
        logger.addHandler(console_handler)

    else:
        logger.setLevel(logging.INFO)

        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)

        formatter = logging.Formatter("[%(levelname)s] %(message)s")


        console_handler.setFormatter(formatter)

        logger.addHandler(console_handler)


    cfg = RunConfig(input_file=args.file, fmt=args.format, sigma=args.sigma,taxonomy=args.taxonomy, output_prefix=args.output,verbose=args.verbose,logfile_path=logfile_path,no_substitutes=args.no_substitutes)
    
    paths = Paths(
        counts_dir=Path("/content/drive/MyDrive/Lab Work/Bayesian_Graph_Metabolic_Paper/Bayesian Priors: ATB Analysis/ATB_Taxonomy/ATB Frequency"),
        onehop_dir=Path("/content/drive/MyDrive/Lab Work/Bayesian_Graph_Metabolic_Paper/KEGG_Neighbors_Generation/ONE_HOP_NEIGHBOR_DATA-2"),
        twohop_dir=Path("/content/drive/MyDrive/Lab Work/Bayesian_Graph_Metabolic_Paper/KEGG_Neighbors_Generation/TWO_HOP_NEIGHBOR_DATA-2"),
        module_neighbor_dir=Path("/content/drive/MyDrive/Lab Work/Bayesian_Graph_Metabolic_Paper/KEGG_Neighbors_Generation/MODULE_ALL_NEIGHBOR_DATA"),
        module_eq_json=Path("/content/drive/MyDrive/Lab Work/Bayesian_Graph_Metabolic_Paper/BLIMMP PATH Probabilities/KEGG_Module_Equations_18AUG26.json"),
        module_json_dir=Path("/content/drive/MyDrive/Lab Work/Bayesian_Graph_Metabolic_Paper/KEGG_Graph_Generation/KEGG_Graphs_Generated_13AUG2026"),
        kofam_ko_list_path=Path("/content/drive/MyDrive/Lab Work/Bayesian_Graph_Metabolic_Paper/ko_list.txt"),
        module_frequencies = Path("/content/drive/MyDrive/Lab Work/Bayesian_Graph_Metabolic_Paper/module_freq.txt"),
        module_reaction_dir=Path("/content/drive/MyDrive/Lab Work/Bayesian_Graph_Metabolic_Paper/Reactions/module_ko_reaction.json"),
        module_descriptions_path=Path("/content/drive/MyDrive/Lab Work/Bayesian_Graph_Metabolic_Paper/kegg_bacteria_modules.json"),
        ko_reaction_path=Path("/content/drive/MyDrive/Lab Work/Bayesian_Graph_Metabolic_Paper/Reactions/ko_reaction.list"), 
        cooccurrence_lookup_dir=Path("/content/drive/MyDrive/Lab Work/Bayesian_Graph_Metabolic_Paper/"),
    )


    BlimmpPipeline(cfg, paths).run()


if __name__ == "__main__":
    main() 