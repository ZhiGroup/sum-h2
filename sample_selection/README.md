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
        │
        ▼  (1) relatedness filter — Kinship ≥ 0.022
 3,318  related keep
        │
        ▼  (2) further random filter down to ≈ 2,158
 2,158  final analysis GRM
```

In this UKB run, step (2) was done by intersecting with a **random ~⅔ discovery cohort** (22,985 / 35,810 → related ∩ discovery = **2,158**). Prefer **`DISC_TARGET_N=2158`** so the final count is exact; `DISC_FRAC` is only a rough alternative and is hard to control.

---

## Filter rule (two steps)

1. **Kinship:** keep every sample in ≥1 KING pair with **Kinship ≥ 0.022** (no upper bound; keep both members).
2. **Size:** randomly filter that related set down to **≈ 2,158**. Use **`DISC_TARGET_N=2158`** (exact). Here we originally used ∩ random ~⅔ discovery — same goal, less controllable.

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
| Discovery / size control | After kinship: **`DISC_TARGET_N=2158`** (preferred) or a random fraction / ID list |

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

# Preferred: exact final n after kinship (random sample from related keep)
KIN0=/path/to/king_output.kin0 \
GRM_PREFIX=/path/to/your_full_grm \
DISC_TARGET_N=2158 \
FORCE=1 bash build_king_over4p5_discovery.sh

# Alternative: random fraction of the full cohort (less controllable)
KIN0=/path/to/king_output.kin0 \
GRM_PREFIX=/path/to/your_full_grm \
RANDOM_DISC=1 \
DISC_FRAC=0.6666667 \
FORCE=1 bash build_king_over4p5_discovery.sh

# Or pass your own discovery ID list (optional DISC_TARGET_N to downsample)
KIN0=/path/to/king_output.kin0 \
DISC_ID=/path/to/ids.txt \
DISC_TARGET_N=2158 \
GRM_PREFIX=/path/to/your_full_grm \
FORCE=1 bash build_king_over4p5_discovery.sh
```

| Env | Meaning |
|---|---|
| `KIN0` | `.kin0` path |
| `GRM_PREFIX` | Full dense GRM prefix |
| **`DISC_TARGET_N`** | **Exact final related keep size** (e.g. `2158`) — random sample from over4p5 related |
| `DISC_ID` | Optional FID/IID list (∩ related; then downsample if `DISC_TARGET_N` set) |
| `RANDOM_DISC=1` | Random subset of full GRM `.grm.id` (use with `DISC_FRAC` if no `DISC_TARGET_N`) |
| `DISC_FRAC` | Fraction of full cohort (default `0.6666667`) — only when not using `DISC_TARGET_N` |
| `DISC_SEED` | RNG seed (default `42`) |
| `SKIP_GRM=1` | Keep lists only |

Writes `king_cutoff_over4p5*.txt`, `sample_ids_*.txt`, and `king_over4p5_gcta_discovery.grm.*`. Keep-list order must match `.grm.id`.

---

## Mistakes to avoid

1. **Stopping after kinship only** — related *n* (here 3,318) is too large; set **`DISC_TARGET_N=2158`** (or an equivalent random size filter) so the final keep is ≈ **2,158**.
2. **Using KING’s unrelated filter** — that drops relateds; HE then has almost no relatedness signal.
3. **Treating the PSEUDO demo as real data** — `pseudo/` is fake `PSEUDO_*` IDs for a path check only.
