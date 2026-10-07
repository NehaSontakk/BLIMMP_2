# ATB Intersections

Several large, regenerated data files are hosted as GitHub Release assets
rather than tracked in git, to keep clones small:

- `CoOccurrence_AllVsAll_domain_level_priors.json` -- all-vs-all co-occurrence
  lookup used by the substitution pipeline
- `ko_matrix.tsv` -- the master KO x genome presence matrix (25,829 KOs x
  21,380 genomes)
- `Lineage_Specific_Data/KO_Intersections_Lineage_Specific/` -- all 20
  lineages' pairwise co-occurrence matrices
- `Lineage_Specific_Data/KO_Matrices/` -- all 20 lineages' per-lineage
  presence matrices

All four are bundled into one compressed archive, split into two parts
(GitHub's per-asset limit is 2GB):
`BLIMMP_2_large_data_for_release.tar.gz.part-aa` and `...part-ab`, on the
[v0.1-data release](https://github.com/NehaSontakke/BLIMMP_2/releases/tag/v0.1-data).

Download and reassemble:
```bash
gh release download v0.1-data --pattern 'BLIMMP_2_large_data_for_release.tar.gz.part-*'
cat BLIMMP_2_large_data_for_release.tar.gz.part-aa BLIMMP_2_large_data_for_release.tar.gz.part-ab \
    > BLIMMP_2_large_data_for_release.tar.gz
tar -xzf BLIMMP_2_large_data_for_release.tar.gz
```
This extracts `CoOccurrence_AllVsAll_domain_level_priors.json` and
`ko_matrix.tsv` into the current directory, and the two `Lineage_Specific_Data`
subdirectories into their matching paths -- run this from
`paper/analysis/Prior_Generation/ATB_KO_Frequencies/` so the paths land
correctly.

**Note:** this supersedes the previous single-file asset
(`ko_intersection_ko_matrix_sampleids_domain_level_priors.tsv.gz`), which
contained domain-level intersection data affected by a `uint8` integer
overflow bug (co-occurrence counts above 255 silently wrapped modulo 256).
That asset has been removed from this release. The archive above contains
the corrected data for all 20 lineages, not just domain-level -- see the
`fix-intersection-overflow` commit for full details of the bug and fix.
