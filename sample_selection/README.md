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

Step (2) historically used ∩ a ~⅔ discovery ID list. For new runs prefer **`DISC_TARGET_N=2158`** (exact random sample from the related keep) or pass your own **`DISC_ID`** list.

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
| Dense GCTA GRM | `prefix.grm.{id,bin,N.bin}` — required for HE/REML; subset to the final keep list |
| Relatedness table | KING `.kin0` **or** threshold the GRM (see below) |
| Size control | **`DISC_TARGET_N`** and/or **`DISC_ID`** (no fraction parameter) |

PSEUDO demo (fake IDs, *n* = 27): [`pseudo/`](pseudo/).

### Skip KING (optional)

KING kinship φ ≈ half GCTA relatedness A. You can select from the GRM alone: keep pairs with **A ≥ ~0.044** (≈ 2 × 0.022), both members, then apply `DISC_TARGET_N` / `DISC_ID`. You still need the dense GRM for HE/REML.

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

# Exact final n after kinship (default if you set nothing: DISC_TARGET_N=2158)
KIN0=/path/to/king_output.kin0 \
GRM_PREFIX=/path/to/your_full_grm \
DISC_TARGET_N=2158 \
FORCE=1 bash build_king_over4p5_discovery.sh

# Or your discovery ID list (related ∩ DISC_ID)
KIN0=/path/to/king_output.kin0 \
DISC_ID=/path/to/discovery_ids.txt \
GRM_PREFIX=/path/to/your_full_grm \
FORCE=1 bash build_king_over4p5_discovery.sh

# Both: sample DISC_TARGET_N from (related ∩ DISC_ID)
KIN0=/path/to/king_output.kin0 \
DISC_ID=/path/to/discovery_ids.txt \
DISC_TARGET_N=2158 \
GRM_PREFIX=/path/to/your_full_grm \
FORCE=1 bash build_king_over4p5_discovery.sh
```

| Env | Meaning |
|---|---|
| `KIN0` | `.kin0` path |
| `GRM_PREFIX` | Full dense GRM prefix |
| **`DISC_TARGET_N`** | Exact final related keep size (e.g. `2158`) |
| **`DISC_ID`** | FID/IID discovery list; final = related ∩ list |
| `DISC_SEED` | RNG seed when sampling `DISC_TARGET_N` (default `42`) |
| `SKIP_GRM=1` | Keep lists only |

Writes `king_cutoff_over4p5*.txt`, `sample_ids_*.txt`, and `king_over4p5_gcta_discovery.grm.*`. Keep-list order must match `.grm.id`.

---

## Mistakes to avoid

1. **Stopping after kinship only** — related *n* (here 3,318) is too large; set **`DISC_TARGET_N=2158`** or provide **`DISC_ID`**.
2. **Using KING’s unrelated filter** — that drops relateds; HE then has almost no relatedness signal.
3. **Treating the PSEUDO demo as real data** — `pseudo/` is fake `PSEUDO_*` IDs for a path check only.
