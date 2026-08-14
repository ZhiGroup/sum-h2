# Sample selection (KING over4p5)

**Step ①:** choose who goes into the analysis GRM. Haseman–Elston needs relatedness, so we **keep related people** (both members of each pair) — we do **not** drop one person per family.

---

## UKB sizes (this run)

| Stage | *n* | Traced from |
|---|---:|---|
| Full cohort | **35,810** | `ukb_all.grm.id` / `ukb_all.fam` |
| Discovery split | **22,985** | `ukb_grm_discovery.grm.id` |
| Related (Kinship ≥ 0.022) | 3,318 | `king_cutoff_over4p5.txt` |
| **Final keep** | **2,158** | `king_cutoff_over4p5_discovery.txt` |

```text
35,810  full cohort
   ├── 22,985  discovery split
   └──  3,318  related (Kinship ≥ 0.022)
              └──  2,158  related ∩ discovery   ← analysis GRM
```

**Aim for final *n* ≈ 2,158.** After the relatedness filter, apply an extra **random subsample** if needed so the keep list lands near that size. A ~⅔ discovery split is what this UKB run used; it is **not** required — any random fraction (or a second downsample of the related set) is fine if the end count is ≈ 2,158.

---

## Filter rule

Keep every sample in ≥1 KING pair with **Kinship ≥ 0.022** (no upper bound). Then intersect with / randomly subsample toward the discovery analysis set.

| Label | KING kinship |
|---|---|
| over5 | ≥ 0.015625 |
| **over4p5 (default)** | **≥ 0.022** |
| over4 | ≥ 0.03125 |

---

## Inputs

| Need | Notes |
|---|---|
| Genotypes | BGEN+`.sample` or PLINK bed (for KING and/or GCTA) |
| Dense GCTA GRM | `prefix.grm.{id,bin,N.bin}` — required for HE/REML; subset to the final keep list |
| Relatedness table | KING `.kin0` **or** threshold the GRM (see below) |
| Discovery / size control | Optional ID list, or random subsample so final *n* ≈ 2,158 |

PSEUDO demo (fake IDs, *n* = 27): [`pseudo/`](pseudo/).

### Skip KING (optional)

KING kinship φ ≈ half GCTA relatedness A. You can select from the GRM alone: keep pairs with **A ≥ ~0.044** (≈ 2 × 0.022), both members, then random-subsample to ≈ 2,158. You still need the dense GRM for HE/REML.

---

## Build genotypes → `.kin0` + GRM

**KING** (PLINK bed only; convert BGEN with plink2 if needed):

```bash
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xzf Linux-king.tar.gz -C king_bin && chmod +x king_bin/king
./king_bin/king -b your_genotypes.bed --kinship --prefix king_output --degree 5
# → king_output.kin0
```

Do **not** use KING’s “unrelated keep” step — we keep relateds. See also [KING Download](https://www.kingrelatedness.com/Download.shtml).

**GCTA dense GRM:**

```bash
gcta64 --bfile your_genotypes --make-grm --out your_full_grm --thread-num 8
# or: gcta64 --bgen your.bgen --sample your.sample --make-grm --out your_full_grm
```

---

## Run the builder

```bash
cd sample_selection

# Random discovery split (default fraction ≈ ⅔; change DISC_FRAC if needed)
KIN0=/path/to/king_output.kin0 \
GRM_PREFIX=/path/to/your_full_grm \
RANDOM_DISC=1 \
FORCE=1 bash build_king_over4p5_discovery.sh

# Or pass your own discovery / size-control ID list
KIN0=/path/to/king_output.kin0 \
DISC_ID=/path/to/ids.txt \
GRM_PREFIX=/path/to/your_full_grm \
FORCE=1 bash build_king_over4p5_discovery.sh
```

| Env | Meaning |
|---|---|
| `KIN0` | `.kin0` path |
| `GRM_PREFIX` | Full dense GRM prefix |
| `DISC_ID` | Optional FID/IID list |
| `RANDOM_DISC=1` | Random subset of GRM `.grm.id` |
| `DISC_FRAC` | Random fraction (default `0.6666667`) — tune so final *n* ≈ 2,158 |
| `DISC_SEED` | RNG seed (default `42`) |
| `SKIP_GRM=1` | Keep lists only |

Writes `king_cutoff_over4p5*.txt`, `sample_ids_*.txt`, and `king_over4p5_gcta_discovery.grm.*`. Keep-list order must match `.grm.id`.

---

## Mistakes to avoid

1. **Final *n* far from ≈ 2,158** — after relatedness filtering, apply an additional random subsample (any fraction, not only ⅔) so the keep list is near this UKB size.
2. **Using KING’s unrelated filter** — that drops relateds; HE then has almost no relatedness signal.
3. **Treating the PSEUDO demo as real data** — `pseudo/` is fake `PSEUDO_*` IDs for a path check only.
