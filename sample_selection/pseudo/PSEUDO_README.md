# PSEUDO demo data (not real subjects)

**These files are synthetic.** IDs are `PSEUDO_0001` … — not UK Biobank or any real cohort.
They exist so a fresh clone can resolve a GRM / keep list and smoke-test paths without sensitive IDs.

| File | Role |
|---|---|
| `king_cutoff_over4p5.txt` | Pseudo related keep (n = 40) |
| `king_cutoff_over4p5_discovery.txt` | Pseudo discovery keep (n = 27 ≈ ⅔) |
| `sample_ids_king_over4p5_discovery.txt` | IIDs only |
| `king_over4p5_gcta_discovery.grm.*` | Tiny dense GCTA GRM (n = 27) |
| `pseudo_phenotype.csv` | Optional toy phenotypes for a path check |

Do **not** treat results from this GRM as scientific. Replace with your own KING + GCTA outputs (see parent `README.md`).

`run_sum_h2.py` picks real data from `repo_local/` or `grm_subset_king/` first; this `pseudo/` folder is the last fallback.
