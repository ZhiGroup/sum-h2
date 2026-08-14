# Performance Benchmarks

Benchmarks comparing heritability-estimation methods on synthetic IDP-like phenotypes.
Results are reported as end-to-end wall time and peak process RSS, sweeping two axes:

- **Trait-count scaling** (`p` = 100 → 100 k, fixed n ≈ 2 158 subjects)
- **Sample-size scaling** (`n` = 719 / 2 158 / 6 474, fixed p = 1 000 traits)

Methods: **HE PCA 128PCs**, **REML**, **tr(P⁻¹G)**, **SVD(P⁻¹G)**, **blockwise mvGREML**.
EVR0.8, QR, and whitened kernel are not used in this repo.

---

## Output figures

| File | Description |
|---|---|
| `IDP_synthetic_trait_scaling_time_resource_v8_hereg.png` | Trait scaling — time + RSS vs p |
| `IDP_synthetic_sample_size_scaling_time_resource_v2_hereg.png` | Sample scaling — time + RSS vs n |

---

## How to reproduce

### Step 1 — Run trait-scaling benchmarks

```bash
cd performance
python3 run_benchmark.py \
    --mode trait_scaling \
    --p-list 100,1000,10000,100000 \
    --use-hereg \
    --out-json results.json
```

### Step 2 — Run sample-size scaling

```bash
python3 run_benchmark.py --mode sample_scaling
```

### Step 3 — Regenerate plots

```bash
python3 plot_scaling.py
```

Reads all data from `results.json`.

---

## Methods (summary)

| Method | Dimensionality | Heritability | Covariates |
|---|---|---|---|
| HE pipeline – PCA 128PCs | Decorrelate with PCA (128→128 for 128-d input) | GCTA `--HEreg` | `--covar` / `--qcovar` in GCTA |
| REML pipeline | same PCA decorrelation | GCTA `--reml` | see fairness note in prior docs |
| tr(P⁻¹G) / SVD(P⁻¹G) | none (full traits) | dense trait-space ops | Python residualization |
| blockwise mvGREML | none | blockwise tr(G)/tr(P) | Python residualization |

Shared front-end: load KING over4p5 GRM → load traits → align samples to covariates.

---

## Summary tables

### Trait scaling (n ≈ 2 158)

| Method | p = 100 | p = 1 k | p = 10 k | p = 100 k time | p = 100 k RSS |
|---|---|---|---|---|---|
| HE pipeline – PCA 128PCs | 6.3 s | 18.7 s | 65.5 s | 525.5 s / 0.15 h | 9.52 GiB |
| REML pipeline | 37.6 s | 55.6 s | 100.1 s | 519.9 s / 0.14 h | 7.59 GiB |
| tr(P⁻¹G) | 1.9 s | 27.4 s | 103.3 s | — | — |
| SVD(P⁻¹G) | 1.9 s | 26.4 s | 274.3 s | — | — |
| blockwise mvGREML | 1.5 s | 6.0 s | 50.5 s | 701.6 s / 0.19 h | 7.33 GiB |

### Sample scaling (p = 1 000 traits)

| Method | n = 719 | n = 2 158 | n = 6 474 |
|---|---|---|---|
| HE pipeline – PCA 128PCs | 14.3 s / 259 MiB | 19.1 s / 415 MiB | 44.1 s / 1 233 MiB |
| REML pipeline | 14.2 s / 255 MiB | 58.1 s / 407 MiB | 982.8 s / 1 049 MiB |
