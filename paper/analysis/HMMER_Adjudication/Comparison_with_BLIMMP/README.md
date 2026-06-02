# Comparison with BLIMMP Pipeline

Validation that `preprocess_domtblout.py` produces consistent adjudication results with the BLIMMP inference pipeline on the same input.

## Test genome

**SAMD00000344** — one genome from the AllTheBacteria reference set.

## Files

```
Comparison_with_BLIMMP/
├── test_SAMD00000344_preprocess.csv   # preprocess_domtblout.py output
├── new_SAMD00000344__BLIMMP_dk.csv    # BLIMMP _dk.csv output
├── plot_adjudication_comparison.py    # script to regenerate the scatter plot
└── adjudication_comparison.png        # score and confidence scatter plots
```

## Why only 259 KOs are compared

Both scripts searched against the **same full KOfam HMM database** — the difference is not about which HMMs were used. The difference is about what each script keeps in its output.

- **preprocess_domtblout.py** keeps ALL 1,782 KOs that had any hmmsearch hit in this genome, regardless of whether they belong to a KEGG module. This includes housekeeping genes, transport proteins, regulatory KOs, and other KOs that KEGG has not assigned to any metabolic module.
- **BLIMMP** only tracks KOs that appear in at least one KEGG module (1,980 KOs total in its universe). It discards non-module KOs entirely because its goal is module completeness prediction, not genome-wide annotation.

The 259 matched KOs are simply the ones that (1) had a hmmsearch hit in this genome AND (2) belong to at least one KEGG module — the only set where both scripts had something to compare.

| Group | Count |
|---|---|
| Total KOs detected by hmmsearch (preprocess) | 1,782 |
| KOs in BLIMMP module universe | 1,980 |
| **KOs in both — used for comparison** | **259** |
| Only in preprocess (real hits, not in any KEGG module) | 1,523 |
| Only in BLIMMP (in a module but no hit in this genome) | 1,721 |

The 1,721 BLIMMP-only KOs are module KOs with no hmmsearch evidence in this genome — BLIMMP assigns them `KO_annotation_conf = 0` so the module probability calculation has a complete KO vector.

## Results

| | |
|---|---|
| KOs compared | 259 |
| HMM score correlation | 1.000 |
| Confidence agrees | 236 / 259 |
| Confidence differs | 23 / 259 |
| Reason for difference | below KOfam threshold |

**Why 23 differ:** BLIMMP sets `KO_annotation_conf = 1.0` for any KO that passes the KOfam score threshold, regardless of the raw noise-floor confidence. `preprocess_domtblout.py` does not apply KOfam thresholds and always retains the raw noise-floor confidence. This is intentional — the preprocessing pipeline produces conservative confidence estimates for reference genome annotation, while BLIMMP applies KOfam thresholds as an additional evidence layer for query genome inference.

## Regenerate the plot

```bash
python plot_adjudication_comparison.py \
    --preprocess test_SAMD00000344_preprocess.csv \
    --blimmp new_SAMD00000344__BLIMMP_dk.csv \
    --out adjudication_comparison.png \
    --sample SAMD00000344
```
