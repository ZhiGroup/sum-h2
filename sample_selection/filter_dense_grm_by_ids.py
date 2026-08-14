#!/usr/bin/env python3
"""
Filter a GCTA dense GRM (.grm.id, .grm.bin, .grm.N.bin) to a subset of samples
given by a keep file (FID IID, same format as .grm.id). Writes the subset in
GCTA binary format (lower triangle) under the output prefix.

Alternatively, use --gen-kin-lists to generate ID lists of samples that have
at least one pair with kinship above a threshold (over 5: >= 0.015625;
over 4.5: >= 0.022; over 4: >= 0.03125).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Kinship thresholds for "over 5", "over 4.5", and "over 4" ID lists
KIN_OVER5_LOWER = 0.015625   # no upper bound
KIN_OVER4P5_LOWER = 0.022    # no upper bound
KIN_OVER4_LOWER = 0.03125   # no upper bound


def get_grm_paths(grm_prefix: Path):
    grm_prefix = Path(grm_prefix)
    grm_id_path = grm_prefix.with_suffix(grm_prefix.suffix + ".grm.id")
    if grm_prefix.suffix:
        grm_bin_path = grm_prefix.with_suffix(grm_prefix.suffix + ".grm.bin")
        grm_N_path = grm_prefix.with_suffix(grm_prefix.suffix + ".grm.N.bin")
    else:
        grm_bin_path = grm_prefix.parent / (grm_prefix.name + ".grm.bin")
        grm_N_path = grm_prefix.parent / (grm_prefix.name + ".grm.N.bin")
    if not grm_id_path.is_file():
        grm_id_path = grm_prefix.parent / (grm_prefix.name + ".grm.id")
        grm_bin_path = grm_prefix.parent / (grm_prefix.name + ".grm.bin")
        grm_N_path = grm_prefix.parent / (grm_prefix.name + ".grm.N.bin")
    return grm_id_path, grm_bin_path, grm_N_path


def gen_kin_id_lists(grm_prefix: Path, out_prefix: Path):
    """Generate ID lists for over 5, over 4.5, and over 4 kinship thresholds."""
    grm_id_path, grm_bin_path, grm_N_path = get_grm_paths(grm_prefix)
    for p in (grm_id_path, grm_bin_path, grm_N_path):
        if not p.is_file():
            print(f"Error: missing {p}", file=sys.stderr)
            sys.exit(1)

    df_full = pd.read_csv(grm_id_path, sep=r"\s+", header=None, names=["FID", "IID"])
    n_full = len(df_full)
    nval_full = n_full * (n_full + 1) // 2

    with open(grm_bin_path, "rb") as f:
        data = np.frombuffer(f.read(), dtype=np.float32, count=nval_full)
    if len(data) != nval_full:
        print(f"Error: .grm.bin has {len(data)} elements, expected {nval_full}", file=sys.stderr)
        sys.exit(1)

    def pos_lt(i: int, j: int) -> int:
        if i < j:
            i, j = j, i
        return i * (i + 1) // 2 + j

    ids_over5 = set()
    ids_over4p5 = set()
    ids_over4 = set()
    for i in range(n_full):
        for j in range(i):
            k = data[pos_lt(i, j)]
            if k >= KIN_OVER5_LOWER:
                fid_i, iid_i = str(df_full.iloc[i]["FID"]).strip(), str(df_full.iloc[i]["IID"]).strip()
                fid_j, iid_j = str(df_full.iloc[j]["FID"]).strip(), str(df_full.iloc[j]["IID"]).strip()
                ids_over5.add((fid_i, iid_i))
                ids_over5.add((fid_j, iid_j))
            if k >= KIN_OVER4P5_LOWER:
                fid_i, iid_i = str(df_full.iloc[i]["FID"]).strip(), str(df_full.iloc[i]["IID"]).strip()
                fid_j, iid_j = str(df_full.iloc[j]["FID"]).strip(), str(df_full.iloc[j]["IID"]).strip()
                ids_over4p5.add((fid_i, iid_i))
                ids_over4p5.add((fid_j, iid_j))
            if k >= KIN_OVER4_LOWER:
                fid_i, iid_i = str(df_full.iloc[i]["FID"]).strip(), str(df_full.iloc[i]["IID"]).strip()
                fid_j, iid_j = str(df_full.iloc[j]["FID"]).strip(), str(df_full.iloc[j]["IID"]).strip()
                ids_over4.add((fid_i, iid_i))
                ids_over4.add((fid_j, iid_j))

    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    base = out_prefix.parent / out_prefix.name

    out_over5 = base.with_name(base.name + "_ids_over5.txt")
    with open(out_over5, "w") as f:
        for (fid, iid) in sorted(ids_over5):
            f.write(f"{fid}\t{iid}\n")
    print(f"Wrote {out_over5} ({len(ids_over5)} IDs, kinship >= {KIN_OVER5_LOWER})")

    out_over4p5 = base.with_name(base.name + "_kin_over4p5.txt")
    with open(out_over4p5, "w") as f:
        for (fid, iid) in sorted(ids_over4p5):
            f.write(f"{fid}\t{iid}\n")
    print(f"Wrote {out_over4p5} ({len(ids_over4p5)} IDs, kinship >= {KIN_OVER4P5_LOWER})")

    out_over4 = base.with_name(base.name + "_ids_over4.txt")
    with open(out_over4, "w") as f:
        for (fid, iid) in sorted(ids_over4):
            f.write(f"{fid}\t{iid}\n")
    print(f"Wrote {out_over4} ({len(ids_over4)} IDs, kinship >= {KIN_OVER4_LOWER})")
    print("Done.")


def main():
    parser = argparse.ArgumentParser(description="Filter dense GRM to samples in keep file or generate kinship ID lists")
    parser.add_argument("grm_prefix", type=Path, help="Full GRM prefix (e.g. .../ukb_all)")
    parser.add_argument("keep_file", type=Path, nargs="?", default=None, help="Keep file: FID IID (optional if --gen-kin-lists)")
    parser.add_argument("-o", "--out", type=Path, required=True, help="Output GRM prefix or base for ID list filenames")
    parser.add_argument("--gen-kin-lists", action="store_true", help="Generate ID lists for over 5 (>=0.015625), over 4.5 (>=0.022), and over 4 (>=0.03125) only")
    args = parser.parse_args()

    grm_prefix = Path(args.grm_prefix)
    out_prefix = Path(args.out)

    if args.gen_kin_lists:
        gen_kin_id_lists(grm_prefix, out_prefix)
        return

    keep_path = args.keep_file
    if keep_path is None:
        print("Error: keep_file required when not using --gen-kin-lists", file=sys.stderr)
        sys.exit(1)
    keep_path = Path(keep_path)

    grm_id_path, grm_bin_path, grm_N_path = get_grm_paths(grm_prefix)
    for p in (grm_id_path, grm_bin_path, grm_N_path):
        if not p.is_file():
            print(f"Error: missing {p}", file=sys.stderr)
            sys.exit(1)
    if not keep_path.is_file():
        print(f"Error: missing {keep_path}", file=sys.stderr)
        sys.exit(1)

    # Full GRM sample list
    df_full = pd.read_csv(grm_id_path, sep=r"\s+", header=None, names=["FID", "IID"])
    n_full = len(df_full)
    id_to_idx = {}
    for idx, row in df_full.iterrows():
        key = (str(row["FID"]).strip(), str(row["IID"]).strip())
        id_to_idx[key] = idx

    # Keep list (preserve order)
    df_keep = pd.read_csv(keep_path, sep=r"\s+", header=None, names=["FID", "IID"])
    if len(df_keep) and str(df_keep.iloc[0]["FID"]).upper() == "FID":
        df_keep = df_keep.iloc[1:].reset_index(drop=True)
    df_keep["FID"] = df_keep["FID"].astype(str).str.strip()
    df_keep["IID"] = df_keep["IID"].astype(str).str.strip()

    keep_ix = []
    missing = []
    for _, row in df_keep.iterrows():
        key = (row["FID"], row["IID"])
        if key in id_to_idx:
            keep_ix.append(id_to_idx[key])
        else:
            missing.append(f"{row['FID']} {row['IID']}")
    if missing:
        print(f"Warning: {len(missing)} IDs in keep file not in GRM (first: {missing[0]})", file=sys.stderr)
    n_keep = len(keep_ix)
    print(f"Keeping {n_keep} samples (full GRM has {n_full})")

    # Load full lower triangle (1D)
    nval_full = n_full * (n_full + 1) // 2
    with open(grm_bin_path, "rb") as f:
        data = np.frombuffer(f.read(), dtype=np.float32, count=nval_full)
    if len(data) != nval_full:
        print(f"Error: .grm.bin has {len(data)} elements, expected {nval_full}", file=sys.stderr)
        sys.exit(1)

    with open(grm_N_path, "rb") as f:
        N_data = np.frombuffer(f.read(), dtype=np.int32, count=nval_full)
    if len(N_data) != nval_full:
        print(f"Error: .grm.N.bin has {len(N_data)} elements, expected {nval_full}", file=sys.stderr)
        sys.exit(1)

    def pos_lt(i: int, j: int) -> int:
        if i < j:
            i, j = j, i
        return i * (i + 1) // 2 + j

    # Write output .grm.id (same order as keep_ix: kept samples only)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    id_out = out_prefix.parent / (out_prefix.name + ".grm.id")
    kept_df = df_keep[df_keep.apply(lambda r: (r["FID"], r["IID"]) in id_to_idx, axis=1)]
    with open(id_out, "w") as f:
        for _, row in kept_df.iterrows():
            f.write(f"{row['FID']}\t{row['IID']}\n")
    print(f"Wrote {id_out}")

    # Extract lower triangle for kept samples and write .grm.bin and .grm.N.bin
    bin_out = out_prefix.parent / (out_prefix.name + ".grm.bin")
    N_out = out_prefix.parent / (out_prefix.name + ".grm.N.bin")
    nval_keep = n_keep * (n_keep + 1) // 2
    out_vals = np.zeros(nval_keep, dtype=np.float32)
    out_N = np.zeros(nval_keep, dtype=np.int32)
    idx_out = 0
    for a in range(n_keep):
        for b in range(a + 1):
            i, j = keep_ix[a], keep_ix[b]
            p = pos_lt(i, j)
            out_vals[idx_out] = data[p]
            out_N[idx_out] = N_data[p]
            idx_out += 1

    with open(bin_out, "wb") as f:
        f.write(out_vals.tobytes())
    with open(N_out, "wb") as f:
        f.write(out_N.tobytes())
    print(f"Wrote {bin_out} ({nval_keep} float32)")
    print(f"Wrote {N_out}")
    print("Done.")


if __name__ == "__main__":
    main()
