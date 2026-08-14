# Sample selection (KING over4p5)

**Step ①:** choose who goes into the analysis GRM. Haseman–Elston needs relatedness, so we **keep related people** (both members of each pair) — we do **not** drop one person per family.

---

## UKB sizes (this run)

| Stage | *n* | Traced from |
|---|---:|---|
| Full cohort | **35,810** | `ukb_all.grm.id` / `ukb_all.fam` |
| Discovery split (historical) | **22,985** | `ukb_grm_discovery.grm.id` |
| Related (Kinship ≥ 0.022) | 3,318 | `king_cutoff_over4p5.txt` |
| **Final keep** | **2,158** | `king_cutoff_over4p5_discovery.txt` |

```text
35,810  full cohort
        │
        ▼  (1) relatedness filter — Kinship ≥ 0.022
 3,318  related keep
        │
        ▼  (2) size filter → ≈ 2,158
 2,158  final analysis GRM
```

Step (2) historically used ∩ a ~⅔ discovery ID list. For new runs prefer **`DISC_TARGET_N=2158`** or pass your own **`DISC_ID`** list.

---

## Filter rule (two steps)

1. **Kinship:** keep every sample in ≥1 KING pair with **Kinship ≥ 0.022** (no upper bound; keep both members).
2. **Size:** either **`DISC_TARGET_N=2158`** or **`DISC_ID=`** your discovery list (final keep = related ∩ list).

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
| Dense GCTA GRM | `prefix.grm.{id,bin,N.bin}` — required for HE/REML |
| Relatedness | KING `.kin0` **(optional)** — or threshold the GRM (Procedure step 2) |
| Size control | **`DISC_TARGET_N`** and/or **`DISC_ID`** |

PSEUDO demo (fake IDs, *n* = 27): [`pseudo/`](pseudo/).

---

## Procedure

### 1. Build the dense GCTA GRM (required)

```bash
gcta64 --bfile your_genotypes --make-grm --out your_full_grm --thread-num 8
# or: gcta64 --bgen your.bgen --sample your.sample --make-grm --out your_full_grm
# → your_full_grm.grm.id / .grm.bin / .grm.N.bin
```

### 2. Relatedness table — KING `.kin0` (optional)

Default path: run KING on PLINK bed (convert BGEN with plink2 if needed):

```bash
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xzf Linux-king.tar.gz -C king_bin && chmod +x king_bin/king
./king_bin/king -b your_genotypes.bed --kinship --prefix king_output --degree 5
# → king_output.kin0
```

Do **not** use KING’s “unrelated keep” step — we keep relateds. See [KING Download](https://www.kingrelatedness.com/Download.shtml).

**Skip KING:** set **`USE_GRM_REL=1`** in step 3 to threshold the dense GRM directly (**A ≥ 0.044** for over4p5 ≈ 2 × 0.022). You still need the GRM from step 1.

### 3. Run the builder (size filter + GRM subset)

Relatedness can come from KING (default) **or** directly from the dense GRM.

```bash
cd sample_selection

# A) KING .kin0 + exact final n
KIN0=/path/to/king_output.kin0 \
GRM_PREFIX=/path/to/your_full_grm \
DISC_TARGET_N=2158 \
FORCE=1 bash build_king_over4p5_discovery.sh

# B) Skip KING — filter relatedness on GCTA GRM (over4p5: A ≥ 0.044 ≈ 2×0.022)
USE_GRM_REL=1 \
GRM_PREFIX=/path/to/your_full_grm \
DISC_TARGET_N=2158 \
FORCE=1 bash build_king_over4p5_discovery.sh

# C) Discovery ID list (related ∩ DISC_ID); works with KIN0 or USE_GRM_REL=1
KIN0=/path/to/king_output.kin0 \
DISC_ID=/path/to/discovery_ids.txt \
GRM_PREFIX=/path/to/your_full_grm \
FORCE=1 bash build_king_over4p5_discovery.sh
```

| Env | Meaning |
|---|---|
| `GRM_PREFIX` | Full dense GRM prefix (**required**) |
| `KIN0` | `.kin0` path (default relatedness source) |
| **`USE_GRM_REL=1`** | Skip KING; keep pairs with GRM **A ≥ 0.044** (over4p5), **0.03125** (over5), **0.0625** (over4) |
| **`DISC_TARGET_N`** | Exact final related keep size (e.g. `2158`; default if unset) |
| **`DISC_ID`** | FID/IID discovery list; final = related ∩ list |
| `DISC_SEED` | RNG seed when sampling `DISC_TARGET_N` (default `42`) |
| `SKIP_GRM=1` | Keep lists only |

Writes `king_cutoff_over4p5*.txt`, `sample_ids_*.txt`, and `king_over4p5_gcta_discovery.grm.*`. Keep-list order must match `.grm.id`.

---

## Mistakes to avoid

1. **Stopping after kinship only** — related *n* (here 3,318) is too large; set **`DISC_TARGET_N=2158`** or provide **`DISC_ID`**.
2. **Using KING’s unrelated filter** — that drops relateds; HE then has almost no relatedness signal.
3. **Treating the PSEUDO demo as real data** — `pseudo/` is fake `PSEUDO_*` IDs for a path check only.
