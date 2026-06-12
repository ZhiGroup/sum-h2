# Validation: h² vs GWAS Loci Count

Validates whether the heritability (h²) estimated by the pipeline correlates with
the number of significant GWAS loci across 22 phenotype models (ViT, CNN, MoCo,
surface mesh, and standard IDPs from UK Biobank).

---

## Output figures

Both figures use **FUMA** positional LD clumping to define independent loci. They
differ in the upstream GWAS method:

| File | x-axis | GWAS method | Models |
|---|---|---|---|
| `scatter_5panel_hereg_fuma.png` | minP loci | minP — minimum p-value across traits per SNP | 22 |
| `scatter_5panel_hereg_jagwas.png` | JAGWAS lmm_only loci | JAGWAS lmm_only — linear mixed model | 21 |

**5 panels (same in both figures):**
1. GCTA HEreg – PCA 128PCs (sum h²)
2. GCTA HEreg – PCA EVR0.8 (mean h² up to EVR≥0.8 cutoff, stored in results.json)
3. GCTA HEreg – QR (mean h² = sum_h2 / n_dims)
4. GCTA REML – PCA 128PCs (reml_sum_h2, from updated pipeline with --covar/--qcovar)
5. blockwise mvGREML – tr(G)/tr(P)

---

## Files

```
validation_h2_vs_loci_num/
├── results.json                  ← all pre-computed results (see structure below)
├── run_pipeline.py               ← runs GCTA HEreg+REML on all 22 models
├── plot_scatter.py               ← generates scatter figures from results.json
├── scatter_5panel_hereg_fuma.png
└── scatter_5panel_hereg_jagwas.png
```

`hereg_runs/` is created locally by `run_pipeline.py` and is not in the repo.

---

## results.json structure

```
results.json
├── h2                            per-model heritability for all 22 models
│   └── <model_key>
│       ├── pca
│       │   ├── sum_h2            HEreg sum h² across 128 PCs
│       │   ├── reml_sum_h2       REML sum h² across 128 PCs
│       │   ├── mean_h2_evr08     HEreg mean h² up to EVR≥0.8 cutoff
│       │   ├── n_pcs_evr08       number of PCs retained for EVR≥0.8
│       │   └── n_pcs             number of PCs used (128)
│       ├── qr
│       │   ├── sum_h2            HEreg sum h² across all k=min(n,p) QR components
│       │   └── n_dims            number of QR dimensions
│       └── blockwise
│           └── trG_over_trP      blockwise mvGREML tr(G)/tr(P) estimate
│
├── loci                          GWAS loci counts for all 22 models
│   └── <model_key>
│       ├── minp                  minP GWAS loci (FUMA clumping), all 22 models
│       └── jagwas                JAGWAS lmm_only loci (FUMA clumping), null for 1 model
│
├── fit_metrics
│   ├── minp_loci                 Pearson r vs minP loci (22 models)
│   │   ├── he_pca_128pcs         r = 0.798, n = 22
│   │   ├── reml_pca_128pcs       r = 0.791, n = 22
│   │   ├── he_pca_evr08          r = 0.808, n = 22
│   │   ├── he_qr_mean            r = 0.783, n = 22
│   │   └── blockwise_trG_trP     r = 0.765, n = 22
│   └── jagwas_lmm_only_loci      Pearson r vs JAGWAS lmm_only loci (21 models)
│       ├── he_pca_128pcs         r = 0.925, n = 21
│       ├── reml_pca_128pcs       r = 0.931, n = 21
│       ├── he_pca_evr08          r = 0.722, n = 21
│       ├── he_qr_mean            r = 0.931, n = 21
│       └── blockwise_trG_trP     r = 0.670, n = 21
│
└── kinship_sensitivity           KING kinship threshold sensitivity (18 models × 7 variants)
    ├── rows
    ├── correlation_summary
    └── scatter_panels
```

---

## How to reproduce

### Prerequisites

- `run_idp_pipeline_king_hereg.py` must be in the repo root (one directory up from this folder) — it is already there.
- All JAGWAS lmm_only loci counts are pre-stored in `results.json["jagwas_loci"]` — no external file needed.

### Step 1 — Run the pipeline (skip if using pre-computed results)

```bash
cd validation_h2_vs_loci_num
python3 run_pipeline.py
```

Writes `hereg_runs/<model>/pca/` and `hereg_runs/<model>/qr/` (skips models already done).

### Step 2 — Regenerate plots

```bash
python3 plot_scatter.py
```

Reads from `results.json` (EVR0.8 values pre-stored; blockwise from `results.json["blockwise"]`).
Writes `scatter_5panel_hereg_fuma.png` and `scatter_5panel_hereg_jagwas.png`.

---

## Models and families

| Family | Colour | Models |
|---|---|---|
| ViT | orange | ViT-T1-PCA, ViT-T1-MAE, ViT-T1-ADNI, ViT-T1-Random, ViT-T1-replicate, ViT-T1-replicate2, ViT-T1-cv-fold1/2/3, ViT-T2-PCA |
| CNN | green | CNN_replicate1, CNN_T1_ADNI |
| MoCo | purple | MoCo_T2, MoCo_T1 |
| Surf | grey | surf_mesh_mask75, surf_mesh_mask25, z_GraphUNet, z_FusionDistill |
| IDP | black | T1_IDP_pheno_discovery, IDP_T1, IDP_T2 |

JAGWAS lmm_only figure has 21 models (`z_graph_fsaverage4_t1_gwas` has no lmm_only result).
