# Sample selection (KING over4p5)

## Historical 2,158-person analysis cohort

The historical 2,158-person cohort remains important for reproducing earlier discovery analyses. It was obtained by taking the participants in at least one KING pair with **Kinship ≥ 0.022** and intersecting them with the historical **22,985-person discovery split**. The intersection produced the 2,158-person cohort. This is historical provenance, not the rule for new temporary test cohorts.

## Temporary deterministic cohort while RAP is unavailable

For now, since the RAP platform is currently unavailable, we cannot enforce the use of a shared test cohort (which we will do after RAP opens). As a temporary solution, users should apply a kinship threshold of **≥ 0.022** and **include** all UK Biobank participants who occur in at least one qualifying KING pair, among participants with both brain MRI and genetic data (August 2020). This produces a finite and deterministic test cohort of approximately **8,000 participants** from **84,361 participants**.

| Stage | *n* | Traced from |
|---|---:|---|
| Full eligible cohort | **84,361** | Participants with both brain MRI and August-2020 genetic data |
| Related (Kinship ≥ 0.022) | **~8,000** | `king_cutoff_over4p5.txt` |

```text
84,361  full eligible cohort
        │
        ▼  include every participant in a KING pair with Kinship ≥ 0.022
 ~8,000  deterministic test cohort
```

The new final keep list is derived directly from the full eligible cohort and the KING cutoff, so the same inputs always produce the same cohort. It does not randomly sample 2,158 participants or enforce the historical discovery split.

## Filter rule

1. Identify every participant in at least one KING pair with **Kinship ≥ 0.022**.
2. Keep every identified participant who is also in the full eligible cohort.
3. Use those remaining eligible related participants as the deterministic test cohort.

| Label | KING kinship |
|---|---:|
| over5 | ≥ 0.015625 |
| **over4p5 (default)** | **≥ 0.022** |
| over4 | ≥ 0.03125 |

## Inputs

| Need | Notes |
|---|---|
| Genotypes | BGEN plus `.sample`, or PLINK bed, for KING and/or GCTA |
| Dense GCTA GRM | `prefix.grm.{id,bin,N.bin}`; the `.grm.id` defines the full eligible cohort |
| Relatedness | KING `.kin0` (default), or threshold the dense GRM directly |

## Procedure

### 1. Build the dense GCTA GRM

```bash
gcta64 --bfile your_genotypes --make-grm --out your_full_grm --thread-num 8
```

### 2. Obtain relatedness

```bash
./king_bin/king -b your_genotypes.bed --kinship --prefix king_output --degree 5
```

Or set `USE_GRM_REL=1` to use the dense-GRM equivalent, **A ≥ 0.044** for over4p5.

### 3. Build the deterministic test cohort and its GRM

```bash
cd sample_selection
KIN0=/path/to/king_output.kin0 GRM_PREFIX=/path/to/your_full_grm FORCE=1 bash build_king_over4p5_discovery.sh

# Or without KING:
USE_GRM_REL=1 GRM_PREFIX=/path/to/your_full_grm FORCE=1 bash build_king_over4p5_discovery.sh
```

| Env | Meaning |
|---|---|
| `GRM_PREFIX` | Full dense GRM prefix (**required**) |
| `KIN0` | KING `.kin0` path |
| `USE_GRM_REL=1` | Use dense-GRM relatedness; over4p5 is **A ≥ 0.044** |
| `FORCE=1` | Recreate the deterministic keep list and subset GRM |
| `SKIP_GRM=1` | Write keep lists only |

The builder writes `king_cutoff_over4p5.txt` (all eligible participants in qualifying KING pairs), `test_cohort_kinship_ge_over4p5.txt` (the final deterministic keep list), and `test_cohort_over4p5_gcta.grm.*` (the subset dense GRM).

## Mistakes to avoid

1. Dropping related participants; this workflow intentionally keeps every eligible participant occurring in at least one qualifying pair.
2. Using `DISC_TARGET_N=2158` or `DISC_ID` for the new cohort; these historical-selection controls are intentionally unsupported.
3. Treating the PSEUDO demo as real data.
