#!/usr/bin/env bash
# =============================================================================
# Build deterministic kinship-selected test-cohort lists + GCTA GRMs
# =============================================================================
#
# Relatedness source:
#   default         KING .kin0  (Kinship ≥ 0.022 for over4p5)
#   USE_GRM_REL=1   dense GCTA GRM (A ≥ 0.044 for over4p5; ≈ 2× KING)
#
# The final cohort includes every eligible participant in a pair at/above the cutoff.
# It is deterministic and never a random fixed-size subset.
#
# Usage
#   KIN0=... GRM_PREFIX=... FORCE=1 bash build_king_over4p5_discovery.sh
#   USE_GRM_REL=1 GRM_PREFIX=... FORCE=1 bash build_king_over4p5_discovery.sh
#   SKIP_GRM=1 ...
#
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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
USE_GRM_REL="${USE_GRM_REL:-0}"
GRM_PREFIX="${GRM_PREFIX:-$GCTA_DIR/ukb_all}"
PYTHON_FILTER="${PYTHON_FILTER:-$SCRIPT_DIR/filter_dense_grm_by_ids.py}"
KING_EXE="${KING_EXE:-$KING_DIR/king_bin/king}"
FORCE="${FORCE:-0}"
SKIP_GRM="${SKIP_GRM:-0}"

need() { [[ -f "$1" ]] || { echo "Missing: $1" >&2; exit 1; }; }

if [[ -n "${DISC_TARGET_N:-}" || -n "${DISC_ID:-}" ]]; then
  echo "DISC_TARGET_N and DISC_ID are no longer supported: the cohort is deterministic from the full GRM ID list." >&2
  exit 2
fi

need "$PYTHON_FILTER"
need "${GRM_PREFIX}.grm.id"

if [[ "$USE_GRM_REL" == "1" ]]; then
  echo "=== Relatedness from GCTA GRM (USE_GRM_REL=1; over4p5 A≥0.044) ==="
  need "${GRM_PREFIX}.grm.bin"
  REL_MODE="grm"
  REL_SRC="$GRM_PREFIX"
else
  if [[ ! -f "$KIN0" ]]; then
    echo "=== KING --kinship (kin0 missing) ==="
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
  REL_MODE="kin0"
  REL_SRC="$KIN0"
fi

echo "=== Cutoff lists + deterministic test-cohort filter ==="
python3 - "$REL_MODE" "$REL_SRC" "$OUT_DIR" "$FORCE" "${GRM_PREFIX}.grm.id" <<'PY'
import sys
from pathlib import Path

import numpy as np
import pandas as pd

rel_mode = sys.argv[1]
rel_src = Path(sys.argv[2])
out_dir = Path(sys.argv[3])
force = sys.argv[4] == "1"
cohort_id_path = Path(sys.argv[5])

THRESHOLDS_KIN0 = {"over5": 0.015625, "over4p5": 0.022, "over4": 0.03125}
THRESHOLDS_GRM = {"over5": 0.03125, "over4p5": 0.044, "over4": 0.0625}

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

def related_from_kin0(kin0_path: Path):
    kin = pd.read_csv(
        kin0_path, sep=r"\s+",
        dtype={"FID1": str, "ID1": str, "FID2": str, "ID2": str},
    )
    for c in ("FID1", "ID1", "FID2", "ID2"):
        kin[c] = kin[c].str.strip()
    print(f"  kin0 pairs: {len(kin):,}  Kinship [{kin['Kinship'].min():.4f}, {kin['Kinship'].max():.4f}]")
    out = {}
    for label, thr in THRESHOLDS_KIN0.items():
        sub = kin.loc[kin["Kinship"] >= thr]
        ids = set(zip(sub["FID1"], sub["ID1"])) | set(zip(sub["FID2"], sub["ID2"]))
        out[label] = {
            "pairs": len(sub),
            "full": sorted(ids, key=lambda t: sort_key(t[0])),
            "thr": thr,
            "thr_name": "Kinship",
        }
    return out

def related_from_grm(grm_prefix: Path):
    grm_id = Path(str(grm_prefix) + ".grm.id")
    grm_bin = Path(str(grm_prefix) + ".grm.bin")
    df = load_fid_iid(grm_id)
    n = len(df)
    nval = n * (n + 1) // 2
    print(f"  scanning GRM {grm_bin.name}  n={n:,}  lower-tri={nval:,} …")
    data = np.fromfile(grm_bin, dtype=np.float32, count=nval)
    if len(data) != nval:
        raise SystemExit(f"GRM size mismatch: got {len(data)}, expected {nval}")

    fids = df["FID"].tolist()
    iids = df["IID"].tolist()
    ids_by = {lab: set() for lab in THRESHOLDS_GRM}
    pairs_by = {lab: 0 for lab in THRESHOLDS_GRM}
    pos = 0
    for i in range(n):
        for j in range(i + 1):
            if i != j:
                a = float(data[pos])
                for lab, thr in THRESHOLDS_GRM.items():
                    if a >= thr:
                        ids_by[lab].add((fids[i], iids[i]))
                        ids_by[lab].add((fids[j], iids[j]))
                        pairs_by[lab] += 1
            pos += 1

    out = {}
    for label, thr in THRESHOLDS_GRM.items():
        out[label] = {
            "pairs": pairs_by[label],
            "full": sorted(ids_by[label], key=lambda t: sort_key(t[0])),
            "thr": thr,
            "thr_name": "GRM_A",
        }
        print(f"  GRM {label}: A>={thr}  pairs={pairs_by[label]:,}  samples={len(ids_by[label]):,}")
    return out

if rel_mode == "kin0":
    related = related_from_kin0(rel_src)
elif rel_mode == "grm":
    related = related_from_grm(rel_src)
else:
    raise SystemExit(f"Unknown REL_MODE={rel_mode}")

thresholds = THRESHOLDS_KIN0 if rel_mode == "kin0" else THRESHOLDS_GRM
cohort_df = load_fid_iid(cohort_id_path)
cohort = list(zip(cohort_df["FID"], cohort_df["IID"]))
print(f"  eligible full cohort: {len(cohort):,}")

for label in thresholds:
    full = related[label]["full"]
    n_pairs = related[label]["pairs"]
    thr_v = related[label]["thr"]
    thr_name = related[label]["thr_name"]
    cohort_set = set(cohort)
    test_keep = [p for p in full if p in cohort_set]
    missing_from_cohort = len(full) - len(test_keep)
    if missing_from_cohort:
        print(f"  WARNING: {missing_from_cohort} related IDs are absent from the full cohort ID list")

    cutoff = out_dir / f"king_cutoff_{label}.txt"
    cutoff_test = out_dir / f"test_cohort_kinship_ge_{label}.txt"
    sids = out_dir / f"sample_ids_test_cohort_kinship_ge_{label}.txt"

    if force or not cutoff.is_file():
        write_keep(cutoff, full, body_sep=" ")
    if force or not cutoff_test.is_file():
        write_keep(cutoff_test, test_keep, body_sep="\t")
    if force or not sids.is_file():
        write_iids(sids, test_keep)

    print(
        f"  {label:7s}  {thr_name}>={thr_v:<8}  pairs={n_pairs:6d}  "
        f"related={len(full):5d}  final_keep={len(test_keep):5d}"
    )
    print(f"           {cutoff.name}")
    print(f"           {cutoff_test.name}")
    print(f"           {sids.name}")
PY

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

echo "=== Subset dense GCTA GRM ==="
echo "--- over4p5 deterministic test cohort ---"
filter_grm "$OUT_DIR/test_cohort_kinship_ge_over4p5.txt" "$OUT_DIR/test_cohort_over4p5_gcta"

echo ""
echo "Done."
echo "  keep : $OUT_DIR/test_cohort_kinship_ge_over4p5.txt"
echo "  GRM  : $OUT_DIR/test_cohort_over4p5_gcta.grm.bin"
