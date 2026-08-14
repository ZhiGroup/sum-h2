#!/usr/bin/env bash
# =============================================================================
# Build KING over4p5 (and over4 / over5) sample lists + discovery GCTA GRMs
# =============================================================================
#
# Lives in AGENT/repo/sample_selection/ (self-contained with filter_dense_grm_by_ids.py).
# See README.md in this folder for the filter rule and expected n = 2,158.
#
# What "over4p5" means
#   Keep every sample that appears in at least one KING pair with
#       Kinship >= 0.022
#   There is NO upper bound (1st-degree relatives are kept).
#
# Thresholds (KING kinship coefficient, no upper bound)
#   over5   >= 0.015625   (~5th degree)
#   over4p5 >= 0.022      (between 4th and 5th)
#   over4   >= 0.03125    (~4th degree)
#
# Expected counts (ukb_all, n=35,810)
#   over4p5        3,318   → discovery  2,158
#   over4          2,262   → discovery  1,470
#   over5         12,643   → discovery  8,126
#
# Usage
#   cd AGENT/repo/sample_selection
#   bash build_king_over4p5_discovery.sh
#   FORCE=1 bash build_king_over4p5_discovery.sh
#   RANDOM_DISC=1 KIN0=... GRM_PREFIX=... FORCE=1 bash build_king_over4p5_discovery.sh
#   SKIP_GRM=1 bash build_king_over4p5_discovery.sh
#
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Public repo/ must not hold real subject IDs. Prefer writing to sibling repo_local/.
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
AGENT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
if [[ -z "${OUT_DIR:-}" ]]; then
  if [[ "$(basename "$REPO_DIR")" == "repo" && -d "$AGENT_ROOT/repo_local/sample_selection" ]]; then
    OUT_DIR="$AGENT_ROOT/repo_local/sample_selection"
    echo "OUT_DIR default → $OUT_DIR  (public repo/ stays ID-free)"
  else
    OUT_DIR="$SCRIPT_DIR"
  fi
fi
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

KING_DIR="${KING_DIR:-/data484_4/txia2/gwas_practice/KING}"
GCTA_DIR="${GCTA_DIR:-/data484_4/txia2/gwas_practice/grm/gcta}"
BFILE="${BFILE:-$GCTA_DIR/ukb_all}"
KIN0="${KIN0:-$KING_DIR/king_output.kin0}"
# Discovery ID list is optional. If unset / missing, or RANDOM_DISC=1, draw a random ~⅔
# of IDs from the full GRM .grm.id (same idea as this repo's discovery split).
DISC_ID="${DISC_ID:-}"
RANDOM_DISC="${RANDOM_DISC:-0}"
DISC_FRAC="${DISC_FRAC:-0.6666667}"
DISC_SEED="${DISC_SEED:-42}"
GRM_PREFIX="${GRM_PREFIX:-$GCTA_DIR/ukb_all}"
PYTHON_FILTER="${PYTHON_FILTER:-$SCRIPT_DIR/filter_dense_grm_by_ids.py}"
KING_EXE="${KING_EXE:-$KING_DIR/king_bin/king}"
FORCE="${FORCE:-0}"
SKIP_GRM="${SKIP_GRM:-0}"

need() { [[ -f "$1" ]] || { echo "Missing: $1" >&2; exit 1; }; }

# Default DISC_ID on this machine only when the historical file exists and RANDOM_DISC is off.
if [[ -z "$DISC_ID" && "$RANDOM_DISC" != "1" && -f "$GCTA_DIR/ukb_grm_discovery.grm.id" ]]; then
  DISC_ID="$GCTA_DIR/ukb_grm_discovery.grm.id"
fi

# -----------------------------------------------------------------------------
# Step 0 (only if kin0 is absent): KING --kinship on bed
# Download KING as in gwas_practice/KING/run_king_pipeline.sh if needed.
# -----------------------------------------------------------------------------
if [[ ! -f "$KIN0" ]]; then
  echo "=== Step 0: KING --kinship (kin0 missing) ==="
  if [[ ! -x "$KING_EXE" ]]; then
    echo "Downloading KING (Linux 64-bit) → $KING_DIR/king_bin ..."
    mkdir -p "$KING_DIR/king_bin"
    wget -q https://www.kingrelatedness.com/Linux-king.tar.gz -O "$KING_DIR/Linux-king.tar.gz"
    tar -xzf "$KING_DIR/Linux-king.tar.gz" -C "$KING_DIR/king_bin"
    KING_EXE="$KING_DIR/king_bin/king"
    chmod +x "$KING_EXE"
  fi
  need "${BFILE}.bed"
  "$KING_EXE" -b "${BFILE}.bed" --kinship --prefix "${KIN0%.kin0}" --degree 5
fi
need "$KIN0"
need "$PYTHON_FILTER"
need "${GRM_PREFIX}.grm.id"

# -----------------------------------------------------------------------------
# Steps 1–2: related keep lists from .kin0, then ∩ discovery (file or random ⅔)
# -----------------------------------------------------------------------------
echo "=== Steps 1–2: KING cutoff lists + discovery filter ==="
python3 - "$KIN0" "${DISC_ID:-}" "$OUT_DIR" "$FORCE" "${GRM_PREFIX}.grm.id" "$RANDOM_DISC" "$DISC_FRAC" "$DISC_SEED" <<'PY'
import random
import sys
from pathlib import Path

import pandas as pd

kin0_path = Path(sys.argv[1])
disc_arg = sys.argv[2].strip()
out_dir = Path(sys.argv[3])
force = sys.argv[4] == "1"
grm_id_path = Path(sys.argv[5])
random_disc = sys.argv[6] == "1"
disc_frac = float(sys.argv[7])
disc_seed = int(sys.argv[8])

THRESHOLDS = {
    "over5": 0.015625,
    "over4p5": 0.022,
    "over4": 0.03125,
}

def load_fid_iid(path: Path):
    df = pd.read_csv(path, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
    if len(df) and str(df.iloc[0]["FID"]).upper() == "FID":
        df = df.iloc[1:].reset_index(drop=True)
    df["FID"] = df["FID"].str.strip()
    df["IID"] = df["IID"].str.strip()
    return df

def sort_key(fid: str):
    return (0, int(fid)) if fid.isdigit() else (1, fid)

def write_keep(path: Path, pairs, body_sep=" "):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("FID IID\n")
        for fid, iid in pairs:
            f.write(f"{fid}{body_sep}{iid}\n")

def write_iids(path: Path, pairs):
    with open(path, "w") as f:
        for _, iid in pairs:
            f.write(f"{iid}\n")

kin = pd.read_csv(
    kin0_path,
    sep=r"\s+",
    dtype={"FID1": str, "ID1": str, "FID2": str, "ID2": str},
)
for c in ("FID1", "ID1", "FID2", "ID2"):
    kin[c] = kin[c].str.strip()
print(f"  kin0 pairs: {len(kin):,}  Kinship [{kin['Kinship'].min():.4f}, {kin['Kinship'].max():.4f}]")

disc_path = Path(disc_arg) if disc_arg else None
use_random = random_disc or disc_path is None or not disc_path.is_file()
if use_random:
    pool_df = load_fid_iid(grm_id_path)
    pool = list(zip(pool_df["FID"], pool_df["IID"]))
    rng = random.Random(disc_seed)
    k = max(1, int(round(len(pool) * disc_frac)))
    disc = set(rng.sample(pool, k=min(k, len(pool))))
    disc_out = out_dir / "discovery_ids_random_2thirds.txt"
    write_keep(disc_out, sorted(disc, key=lambda t: sort_key(t[0])), body_sep="\t")
    print(f"  discovery: random {disc_frac:.3f} of GRM IDs → {len(disc):,}  (seed={disc_seed})")
    print(f"           wrote {disc_out.name}")
else:
    disc_df = load_fid_iid(disc_path)
    disc = set(zip(disc_df["FID"], disc_df["IID"]))
    print(f"  discovery IDs ({disc_path}): {len(disc):,}")

for label, thr in THRESHOLDS.items():
    sub = kin.loc[kin["Kinship"] >= thr]
    ids = set(zip(sub["FID1"], sub["ID1"])) | set(zip(sub["FID2"], sub["ID2"]))
    full = sorted(ids, key=lambda t: sort_key(t[0]))
    disc_keep = [p for p in full if p in disc]

    cutoff = out_dir / f"king_cutoff_{label}.txt"
    cutoff_d = out_dir / f"king_cutoff_{label}_discovery.txt"
    sids = out_dir / f"sample_ids_king_{label}_discovery.txt"

    if force or not cutoff.is_file():
        write_keep(cutoff, full, body_sep=" ")
    if force or not cutoff_d.is_file():
        write_keep(cutoff_d, disc_keep, body_sep="\t")
    if force or not sids.is_file():
        write_iids(sids, disc_keep)

    print(
        f"  {label:7s}  Kinship>={thr:<8}  pairs={len(sub):6d}  "
        f"samples={len(full):5d}  discovery={len(disc_keep):5d}"
    )
    print(f"           {cutoff.name}")
    print(f"           {cutoff_d.name}")
    print(f"           {sids.name}")
PY

# -----------------------------------------------------------------------------
# Step 3: subset the full GCTA dense GRM to the keep lists
#   king_cutoff_over4p5.txt            → king_over4p5_gcta.grm.*
#   king_cutoff_over4p5_discovery.txt  → king_over4p5_gcta_discovery.grm.*
#   king_cutoff_over4_discovery.txt    → king_over4_gcta_discovery.grm.*
# Sample order in .grm.id follows the keep file (header skipped).
# -----------------------------------------------------------------------------
if [[ "$SKIP_GRM" == "1" ]]; then
  echo "SKIP_GRM=1: not writing .grm.bin"
  exit 0
fi

need "${GRM_PREFIX}.grm.id"
need "${GRM_PREFIX}.grm.bin"
need "${GRM_PREFIX}.grm.N.bin"

filter_grm() {
  local keep="$1"
  local out="$2"
  if [[ "$FORCE" != "1" && -f "${out}.grm.bin" ]]; then
    echo "  (${out}.grm.bin exists, skip)"
    return
  fi
  python3 -u "$PYTHON_FILTER" "$GRM_PREFIX" "$keep" -o "$out"
}

echo "=== Step 3: filter dense GCTA GRM ==="
echo "--- over4p5 full (3,318) ---"
filter_grm "$OUT_DIR/king_cutoff_over4p5.txt" "$OUT_DIR/king_over4p5_gcta"

echo "--- over4p5 discovery (2,158) ---"
filter_grm "$OUT_DIR/king_cutoff_over4p5_discovery.txt" "$OUT_DIR/king_over4p5_gcta_discovery"

echo "--- over4 discovery (1,470) ---"
filter_grm "$OUT_DIR/king_cutoff_over4_discovery.txt" "$OUT_DIR/king_over4_gcta_discovery"

echo ""
echo "Done."
echo "  keep : $OUT_DIR/king_cutoff_over4p5_discovery.txt"
echo "  GRM  : $OUT_DIR/king_over4p5_gcta_discovery.grm.bin"
echo "  GRM  : $OUT_DIR/king_over4_gcta_discovery.grm.bin"
