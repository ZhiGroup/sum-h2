# Validation: h² vs GWAS Loci Count

Validates whether the heritability (h²) estimated by the pipeline correlates with
the number of significant GWAS loci across phenotype models (ViT, CNN, MoCo,
surface mesh, and standard IDPs from UK Biobank).

Loci are from **JAGWAS lmm_only** (linear mixed model) with FUMA positional LD clumping
(21 models with available lmm_only results).

---

## Output figure

| File | x-axis | Models |
|---|---|---|
| `scatter_3panel_hereg_jagwas.png` | JAGWAS lmm_only loci | 21 |

**3 panels:**

1. GCTA HEreg – PCA 128PCs (sum h²)
2. GCTA REML – PCA 128PCs (reml_sum_h2)
3. blockwise mvGREML – tr(G)/tr(P)

---

## Files

```
validation_h2_vs_loci_num/
├── results.json                      ← pre-computed h2 / loci / fit metrics
├── run_pipeline.py                   ← runs GCTA HEreg+REML (PCA) on all models
├── plot_scatter.py                   ← generates the 3-panel JAGWAS scatter
└── scatter_3panel_hereg_jagwas.png
```

`hereg_runs/` is created locally by `run_pipeline.py` and is not in the repo.

---

## How to reproduce

```bash
cd validation_h2_vs_loci_num
python3 run_pipeline.py      # PCA HEreg+REML for each model
python3 plot_scatter.py      # writes scatter_3panel_hereg_jagwas.png
```

`run_pipeline.py` calls `../run_sum_h2.py` with `--dim-reduction pca --n_pca 128 --h2-method both`.

---

## Correlation summary (from results.json)

Pearson r vs JAGWAS lmm_only loci (n = 21): HE PCA 0.925 · REML 0.931 · blockwise 0.670.
