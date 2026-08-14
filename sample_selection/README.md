# KING over4p5 sample selection

**Step ①** of the pipeline: choose which subjects go into the GRM.

Haseman–Elston needs relatedness, so this step **keeps related people**. It does **not** drop one person from each related pair.

### UKB cohort sizes (this run — count the files)

| Stage | *n* | Source file |
|---|---:|---|
| **Full cohort** | **35,810** | `gwas_practice/grm/gcta/ukb_all.grm.id` (same *n* as `ukb_all.fam`) |
| **Discovery (~⅔)** | **22,985** | `gwas_practice/grm/gcta/ukb_grm_discovery.grm.id` |
| Related (Kinship ≥ 0.022) | 3,318 | `king_cutoff_over4p5.txt` |
| **Final analysis set** | **2,158** | `king_cutoff_over4p5_discovery.txt` / `king_over4p5_gcta_discovery.grm.id` |

```text
35,810  full GRM cohort          ← ukb_all.grm.id
   │
   ├── 22,985  discovery (~⅔)    ← ukb_grm_discovery.grm.id
   └──  3,318  related (Kinship ≥ 0.022 on full cohort)
              │
              └──  2,158  related ∩ discovery   ← target keep / HE GRM
```

**Target:** whatever filter you use, the **final selected *n* should stay ≈ 2,158** (same relatedness / discovery rule on a similar-sized cohort). Do not shrink or expand this set arbitrarily if you want results comparable to this UKB run.

You can start from **only PLINK bed/bim/fam or BGEN + sample**. This folder documents how to produce the pieces the builder needs, then run the filter.

---

## End-to-end: what you need and what you get

### Minimum starting point

| You have | What to do first |
|---|---|
| PLINK `.bed/.bim/.fam` | Go to [A] for `.kin0` (or skip KING — see below), [B] for dense GRM |
| BGEN + `.sample` | Convert to bed (or feed GCTA `--bgen`), then [A] + [B] |
| Already have `.kin0` + dense GRM | Skip to [C] |
| Dense GCTA GRM only | You may **skip KING** (see below) |

### Inputs the builder uses

| Input | Required? | Notes |
|---|---|---|
| **KING `.kin0`** | Default yes | Pairwise kinship (`FID1 ID1 FID2 ID2 … Kinship`). Optional if you threshold the GCTA GRM instead |
| **Full dense GRM** | **Yes** | GCTA `prefix.grm.{id,bin,N.bin}` on the same samples |
| **Discovery ID list** | **No** | If omitted, take a **random ~⅔** of the cohort (GRM `.grm.id`; see [C]) |

### Skipping KING (optional)

KING kinship φ ≈ **half** the GCTA GRM relatedness A (empirically A ≈ 2φ). You can skip KING and select related samples from the **dense GRM** alone:

- Keep every sample in ≥1 pair with **A ≥ ~0.044** (≈ 2 × 0.022), no upper bound; keep **both** members  
- Then ∩ discovery (~⅔)  
- Subset the GRM to that keep list  

You still need the **dense GCTA GRM** for HE/REML — `.kin0` cannot replace it. The default builder uses KING so the published rule stays `Kinship ≥ 0.022`.

**Whichever path you take, aim for final *n* ≈ 2,158** on this UKB-sized discovery cohort (22,985 → related ∩ discovery).

### Discovery vs KING

Discovery ≈ random **~⅔ of the cohort** (here 22,985 / 35,810). Intersect with the related keep list. What matters is ending at **related ∩ discovery ≈ 2,158**.

### Outputs

Real keep lists / GRMs have subject IDs — keep them in **`AGENT/repo_local/sample_selection/`** (not public GitHub). This public tree ships scripts plus a **PSEUDO** demo under [`pseudo/`](pseudo/) (fake `PSEUDO_*` IDs only).

| File | Content |
|---|---|
| `king_cutoff_over4p5.txt` | Related samples (Kinship ≥ 0.022) |
| `king_cutoff_over4p5_discovery.txt` | Related ∩ discovery (~⅔) |
| `sample_ids_king_over4p5_discovery.txt` | Discovery IIDs only |
| `king_over4p5_gcta_discovery.grm.{id,bin,N.bin}` | Dense GRM subset |

`run_sum_h2.py` resolves GRM in order: `repo_local/sample_selection/` → `grm_subset_king/` → `sample_selection/` → **`sample_selection/pseudo/`** (demo).

---

## A. Get `.kin0` (from bed or BGEN)

KING reads **PLINK bed**, not BGEN. Reference: [`gwas_practice/KING/run_king_pipeline.sh`](/data484_4/txia2/gwas_practice/KING/run_king_pipeline.sh) and [KING Download](https://www.kingrelatedness.com/Download.shtml).

### 1. Download KING (Linux 64-bit)

Same pattern as `gwas_practice/KING`:

```bash
mkdir -p king_bin
wget https://www.kingrelatedness.com/Linux-king.tar.gz -O Linux-king.tar.gz
tar -xzf Linux-king.tar.gz -C king_bin
chmod +x king_bin/king
```

### 2. If you only have BGEN → make a bed for KING

```bash
plink2 --bgen your.bgen ref-first \
  --sample your.sample \
  --make-bed \
  --out your_genotypes
# → your_genotypes.bed/.bim/.fam
```

(If you already have bed/bim/fam, skip this.)

### 3. Run KING kinship

```bash
./king_bin/king -b your_genotypes.bed \
  --kinship --prefix king_output --degree 5
# → king_output.kin0
```

Or reuse the helper (defaults to our UKB bed; override `BFILE`):

```bash
cd /data484_4/txia2/gwas_practice/KING
BFILE=/path/to/your_genotypes bash run_king_pipeline.sh
```

For **this** sample-selection protocol you only need Step 1 of that script (the `.kin0`). Do **not** use its “unrelated keep” step — over4p5 **keeps** relateds.

---

## B. Get the full dense GRM (**required**)

The discovery GRM is a **subset** of a full GCTA dense GRM. You must build that full GRM once from genotypes.

### From PLINK bed

```bash
gcta64 --bfile your_genotypes \
  --make-grm \
  --out your_full_grm \
  --thread-num 8
# → your_full_grm.grm.id / .grm.bin / .grm.N.bin
```

### From BGEN (+ sample)

```bash
gcta64 --bgen your.bgen \
  --sample your.sample \
  --make-grm \
  --out your_full_grm \
  --thread-num 8
```

(Optional: `--extract pruned.snplist` if you want LD-pruned SNPs.)

Same sample IDs should appear in the GRM `.grm.id` and in the KING / bed / BGEN cohort.

---

## C. Discovery ≈ random ⅔ (ID list **not** required)

**Discovery** is a roughly random ~⅔ split of the analysis cohort (rest can be replication). Not phenotype-defined.

1. `RANDOM_DISC=1` — builder draws ~⅔ of GRM `.grm.id`, or  
2. Pass an existing FID/IID list.

Intersect with the related keep (Kinship ≥ 0.022). UKB: 3,318 → **2,158**. On your data, *n* differs; the filter rule stays the same.

---

## The filter rule

Keep every sample in **at least one KING pair** with

```text
Kinship >= 0.022
```

- No upper bound (1st-degree stays).
- Keep **both** members of each such pair.
- Intersect with discovery (~⅔ of cohort, or your ID list).
- Sort by FID (numeric).

| Label | KING kinship |
|---|---|
| over5 | ≥ 0.015625 |
| **over4p5 (default)** | **≥ 0.022** |
| over4 | ≥ 0.03125 |

---

## Run the builder

**Your own genotypes (no discovery file — random ~⅔ of GRM IDs):**

```bash
cd sample_selection

KIN0=/path/to/king_output.kin0 \
GRM_PREFIX=/path/to/your_full_grm \
RANDOM_DISC=1 \
FORCE=1 bash build_king_over4p5_discovery.sh
```

Writes `discovery_ids_random_2thirds.txt`, then the over4p5 keep lists + discovery GRM.

**With an existing discovery ID list:**

```bash
KIN0=/path/to/king_output.kin0 \
DISC_ID=/path/to/discovery_ids.txt \
GRM_PREFIX=/path/to/your_full_grm \
FORCE=1 bash build_king_over4p5_discovery.sh
```

| Env var | Meaning |
|---|---|
| `KIN0` | Path to `.kin0` (default path; omit only if you select relateds from the GRM yourself) |
| `GRM_PREFIX` | Full dense GRM prefix (required for GRM subset) |
| `DISC_ID` | Optional FID/IID discovery list |
| `RANDOM_DISC=1` | Ignore `DISC_ID`; sample ~⅔ of GRM `.grm.id` |
| `DISC_FRAC` | Fraction for random split (default `0.6666667`) |
| `DISC_SEED` | RNG seed (default `42`) |
| `BFILE` | bed prefix if `.kin0` is missing (script can download KING and run kinship) |
| `SKIP_GRM=1` | Keep lists only (still need `.kin0`) |

On this machine, defaults point at `gwas_practice/KING/king_output.kin0` and `ukb_all` GRM / discovery IDs when present. If you run the builder from public `repo/sample_selection/`, **`OUT_DIR` defaults to `AGENT/repo_local/sample_selection/`** so real IDs are not written into the GitHub tree.

---

## Quick check

| Source | *n* |
|---|---:|
| UKB real (`repo_local/`) | 2,158 |
| **PSEUDO demo** (`pseudo/`) | **27** (fake `PSEUDO_*` IDs) |

Keep-list IDs and `.grm.id` must match **in the same order**.

---

## Common mistakes

| Mistake | What happens |
|---|---|
| Final *n* far from **≈ 2,158** | Not comparable to this UKB HE run — check discovery ≈ ⅔ and relatedness threshold |
| Skip dense GRM / invent GRM from `.kin0` only | Subset GCTA `--make-grm` output instead |
| Run KING’s unrelated filter | Drops relateds; HE has almost no relatedness |
| Kinship in `[0.015625, 0.125]` | Drops 1st-degree |
| Treat **PSEUDO** demo results as real | Demo GRM is synthetic — replace with your data |
| Commit real keep lists / `.grm.id` to GitHub | Use `repo_local/` for real IDs |

Software: [KING](https://www.kingrelatedness.com/Download.shtml), PLINK 2 (BGEN→bed), GCTA ≥ 1.94 (`--make-grm`).
