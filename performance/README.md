# Performance Benchmarks

Benchmarks comparing seven heritability-estimation methods on synthetic IDP-like phenotypes.
Results are reported as end-to-end wall time and peak process RSS, sweeping two axes:

- **Trait-count scaling** (`p` = 100 → 100 k, fixed n ≈ 2 158 subjects)
- **Sample-size scaling** (`n` = 719 / 2 158 / 6 474, fixed p = 1 000 traits)

---

## Output figures

| File | Description |
|---|---|
| `IDP_synthetic_trait_scaling_time_resource_v8_hereg.png` | Trait scaling — time + RSS vs p |
| `IDP_synthetic_sample_size_scaling_time_resource_v2_hereg.png` | Sample scaling — time + RSS vs n |

---

## How to reproduce

### Step 1 — Run all trait-scaling benchmarks

```bash
cd performance
python3 run_benchmark.py \
    --mode trait_scaling \
    --p-list 100,1000,10000,100000 \
    --isolate-benchmark-processes \
    --out-json results.json
```

This runs **HE PCA 128PCs**, REML, tr(P⁻¹G), SVD(P⁻¹G), and blockwise mvGREML.

To also run EVR0.8 and QR with GCTA HEreg backend:

```bash
# EVR0.8 variant
python3 run_benchmark.py --mode trait_scaling --only-evr08 --use-hereg \
    --p-list 100,1000,10000,100000

# QR variant
python3 run_benchmark.py --mode trait_scaling --only-qr-he --use-hereg \
    --p-list 100,1000,10000,100000
```

### Step 2 — Run sample-size scaling

```bash
python3 run_benchmark.py --mode sample_scaling
```

### Step 3 — Regenerate plots

```bash
python3 plot_scaling.py
```

Reads all data from `results.json`. No external file dependencies.

---

## results.json structure

```
results.json
├── trait_scaling
│   ├── p_values            [100, 1000, 10000, 100000]
│   └── rows                one dict per p; keys: p, feature_dir, <method_key>…
│       ├── HE_pipeline_PCA_capped_e2e
│       ├── HE_EVR08_e2e
│       ├── QR_HE_e2e
│       ├── REML_PCA_capped_e2e
│       ├── tr_PinvG_e2e
│       ├── SVD_PinvG_e2e
│       └── blockwise_mvGREML_trG_trP_e2e
│
├── sample_scaling
│   ├── p_fixed             1000
│   └── rows                one dict per n (n_small_1x3, n_mid_current, n_large_3x)
│       └── <same method keys>
│
└── dense_trace_extrapolation
    └── estimate_100k       memory + time estimates for p=100k dense ops
```

---

## Methods: what is timed and how covariates are handled

All timing uses `time.perf_counter()` (wall clock from function entry to return).
Peak RSS uses `resource.getrusage(RUSAGE_SELF).ru_maxrss` — the **Python driver process
only**; GCTA subprocess memory is not captured (consistently so across all methods).
Each method runs in a **fresh isolated subprocess** and benchmarks run **sequentially**.

### Shared front-end (all pipeline methods)

Every method starts from the same on-disk inputs and runs through three setup steps:

1. **Load GRM** from GCTA binary (KING kinship, `over4p5` threshold)
2. **Load `p` trait CSVs** into a float64 numpy array (`load_Y_from_feature_dir`)
3. **Align samples** to the covariate file — filters to subjects present in both GRM
   and the covariate file; **no Python-side residualization at this step**
   (`align_Y_to_covariates` is sample filtering only)

---

### 1. HE pipeline – PCA 128PCs (`HE_pipeline_PCA_capped_e2e`)

**Function:** `bench_he_pca_capped_e2e`

**Timing breakdown** (approximate at p = 1 000):

| Step | ~% of wall | Detail |
|---|---|---|
| Load GRM | 15% | n × n binary |
| Load p trait CSVs | 35% | one file per feature, float64 |
| Align samples to covariate file | 50% | index matching only |
| PCA on raw Y | 78% | randomized SVD when n_comp < 0.8·min(n,p) |
| GCTA HEreg × min(128, p) PCs | 100% | 8 parallel jobs |

**Dimensionality:** `n_comp = min(128, n_sub − 1, p)`

**Covariate correction:** `--covar` (categorical) and `--qcovar` (quantitative) passed to
GCTA `--HEreg`. No Python-side residualization before PCA — GCTA handles it internally
on the raw PC scores.

---

### 2. HE pipeline – PCA EVR0.8 (`HE_EVR08_e2e`)

**Function:** `bench_he_evr08_e2e` with `--only-evr08 --use-hereg`

Same front-end as PCA 128PCs. Difference: a **full SVD** up to `min(n − 1, p)` components
is computed to find the exact cutoff where cumulative explained variance ratio ≥ 0.8.
Retained k varies with p. Slower than 128PCs at large p because full-rank SVD is required.

**Covariate correction:** same — `--covar/--qcovar` inside GCTA HEreg.

---

### 3. HE pipeline – QR (`QR_HE_e2e`)

**Function:** `bench_qr_he_e2e` with `--only-qr-he --use-hereg`

Same front-end. Replaces PCA with `np.linalg.qr(Y, mode="reduced")`. Retains **all**
k = min(n, p) orthonormal Q-columns with no cap. At n < p (p ≥ 10 k), k ≈ n − 5 ≈ 2 153.

**Covariate correction:** same — `--covar/--qcovar` inside GCTA HEreg on Q-columns.

---

### 4. REML pipeline (`REML_PCA_capped_e2e`)

**Function:** `bench_reml_pca_capped_e2e`

Same front-end (load GRM + traits, align samples, PCA on raw Y — identical to HE PCA
128PCs). After PCA, writes 128 PC-score CSV files to a temp directory, then calls
`_run_reml_for_feature` (from `run_idp_pipeline_king_reml.py`) in parallel (8 jobs)
which runs `gcta --reml` on each PC.

**Covariate correction: none in this benchmark.** `_run_reml_for_feature` does not
accept or pass `--covar/--qcovar` — GCTA REML runs on raw PC scores with no covariate
adjustment. This differs from the HEreg pipeline where GCTA applies `--covar/--qcovar`
correction. The timing comparison is therefore fair (same I/O and PCA cost) but the
statistical treatment is asymmetric.

**Why REML is slower at small p:** GCTA `--reml` uses iterative AI-REML (typically
5–20 Newton steps per trait, O(n²) each). `--HEreg` is a single-pass quadratic
regression. At p = 100, REML takes ~38 s vs ~7 s for HEreg because the iterative
overhead dominates; at p = 100 k, both are bottlenecked by loading 100 k trait files.

---

### 5. tr(P⁻¹G) (`tr_PinvG_e2e`)

**Function:** `bench_dense_ops_e2e(kind='trace')`

Hutchinson stochastic trace estimator of `tr(P⁻¹G)` on the full n × n GRM. **No
p = 100 k result** — the n × p working matrix exceeds available memory at that scale.

**Covariate handling:** Python-side residualization (`residualize_raw_traits`) is applied
to Y before building G and P. This is required because these methods construct explicit
covariance matrices and cannot delegate correction to an external tool like GCTA. The
residualization step is included in the timed wall clock. No PCA is performed; the full
p-trait covariance structure is used directly.

---

### 6. SVD(P⁻¹G) (`SVD_PinvG_e2e`)

**Function:** `bench_dense_ops_e2e(kind='svd')`

Truncated SVD of `P⁻¹G`. Same front-end as tr(P⁻¹G): Python residualization included
in timing, no PCA. **No p = 100 k result.**

---

### 7. blockwise mvGREML (`blockwise_mvGREML_trG_trP_e2e`)

**Function:** `bench_blockwise_mvgreml_e2e`

Blockwise estimation of `tr(G)/tr(P)`. Same Python residualization front-end as
tr/SVD. Avoids forming the full dense P⁻¹G so it scales to p = 100 k. Produces one
aggregate h² estimate per model, not per-component values. No PCA.

---

## Fairness summary

| Criterion | HE PCA/EVR/QR | REML | tr / SVD / blockwise |
|---|---|---|---|
| Same trait data source | ✓ | ✓ | ✓ |
| Same GRM | ✓ | ✓ | ✓ |
| Same sample alignment step | ✓ | ✓ | ✓ |
| Python-side covariate residualization | ✗ (GCTA handles it) | ✗ (not applied in benchmark) | ✓ (required before building covariance matrices) |
| Isolated fresh subprocess per method | ✓ | ✓ | ✓ |
| Sequential execution | ✓ | ✓ | ✓ |
| Produces per-component h² | ✓ | ✓ | ✗ (aggregate only) |

**Key asymmetry — REML covariate correction:** In the validation scatter plot
(`validation_h2_vs_loci_num/`), REML is run via `run_idp_pipeline_king_hereg.py
--h2-method both` which *does* pass `--covar/--qcovar`. In this performance benchmark,
REML is run without covariate correction. The timing comparison is still fair
(same I/O + PCA cost on both sides) but the statistical scope differs.

**Peak RSS note:** `max_rss_mb` is measured via `resource.getrusage(RUSAGE_SELF)` on
the Python driver process only. Memory used by GCTA child processes (which independently
load the GRM and phenotype files) is **not included** for any method — so all RSS
figures consistently undercount the true peak. The undercount is similar in magnitude
across all pipeline methods.

---

## Summary tables

### Trait scaling (n ≈ 2 158)

| Method | p = 100 | p = 1 k | p = 10 k | p = 100 k time | p = 100 k RSS |
|---|---|---|---|---|---|
| HE pipeline – PCA 128PCs | 6.7 s | 22.0 s | 75.0 s | 525.5 s / 0.15 h | 9.52 GiB |
| HE pipeline – PCA EVR0.8 | 3.3 s | 50.7 s | 333.5 s | 1235.9 s / 0.34 h | 12.13 GiB |
| HE pipeline – QR | 3.7 s | 27.6 s | 98.7 s | 527.1 s / 0.15 h | 7.26 GiB |
| REML pipeline | 37.6 s | 59.1 s | 107.1 s | 519.9 s / 0.14 h | 7.59 GiB |
| tr(P⁻¹G) | 2.2 s | 26.7 s | 104.8 s | — | — |
| SVD(P⁻¹G) | 2.1 s | 27.0 s | 245.8 s | — | — |
| blockwise mvGREML | 2.0 s | 7.1 s | 60.1 s | 701.6 s / 0.19 h | 7.33 GiB |

p = 100 k values for tr/SVD are extrapolated (O(p³) kernel; full dense matrices require
~320 GiB temporaries — see `dense_trace_extrapolation` in `results.json`).

### Sample scaling (p = 1 000 traits)

| Method | n = 719 | n = 2 158 | n = 6 474 |
|---|---|---|---|
| HE pipeline – PCA 128PCs | 14.3 s / 259 MiB | 19.1 s / 415 MiB | 44.1 s / 1 233 MiB |
| REML pipeline | 14.2 s / 255 MiB | 58.1 s / 407 MiB | 982.8 s / 1 049 MiB |

HE regression scales linearly in n; REML scales super-linearly (O(n²) per AI-REML
iteration).
