#!/usr/bin/env python3
"""
Run GCTA HEreg+REML pipeline on all 22 phenotype models and save results.

Reads model list and loci counts from results.json["kinship_sensitivity"]["rows"].
If a model has already been run (arena_manifest.json exists in hereg_runs/), it is
skipped. After all models finish, updates hereg_runs/hereg_summary.json.

Run plot_scatter.py separately to regenerate scatter figures from results.json.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE           = Path(__file__).resolve().parent
RESULTS_JSON   = HERE / "results.json"
HEREG_RUNS_DIR = HERE / "hereg_runs"          # created locally; not in repo
HEREG_SCRIPT   = HERE.parent / "run_idp_pipeline_king_hereg.py"
PYTHON         = sys.executable


# ---------------------------------------------------------------------------
# Load model list from results.json
# ---------------------------------------------------------------------------

_store = json.loads(RESULTS_JSON.read_text())
loci_rows = _store["kinship_sensitivity"]["rows"]
models = [
    (r["raw_label"], r["original_feature_dir"], float(r["minp_loci_count"]))
    for r in loci_rows
    if not r.get("source_row_error")
    and r.get("original_feature_dir")
    and r.get("minp_loci_count") is not None
]
print(f"[pipeline] {len(models)} models loaded from results.json")


# ---------------------------------------------------------------------------
# Pipeline runner
# ---------------------------------------------------------------------------

def _safe(label: str) -> str:
    return re.sub(r"[/\\]", "_", label)


def run_hereg(label: str, feature_dir: str, mode: str) -> dict:
    """
    mode 'pca': --h2-method both (HEreg + REML, writes EVR JSON)
    mode 'qr':  --h2-method hereg only
    Returns manifest dict or {"error": ...}.
    """
    out_root = HEREG_RUNS_DIR / _safe(label) / mode
    manifest_path = out_root / "arena_manifest.json"

    if manifest_path.is_file():
        print(f"[pipeline] skip (cached): {label}/{mode}")
        return json.loads(manifest_path.read_text())

    out_root.mkdir(parents=True, exist_ok=True)
    cmd = [PYTHON, str(HEREG_SCRIPT),
           "--phenotype_csv", feature_dir,
           "--output_root", str(out_root)]
    if mode == "qr":
        cmd += ["--dim-reduction", "qr", "--h2-method", "hereg"]
    else:
        cmd += ["--h2-method", "both"]

    print(f"[pipeline] running {label}/{mode} ...")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[pipeline] ERROR {label}/{mode}:\n{res.stderr[-400:]}", flush=True)
        return {"error": res.stderr[-300:]}

    if manifest_path.is_file():
        return json.loads(manifest_path.read_text())
    return {"error": "manifest not written"}


# ---------------------------------------------------------------------------
# Run all models (parallel, 4 workers)
# ---------------------------------------------------------------------------

_summary_path = HEREG_RUNS_DIR / "hereg_summary.json"
results: dict[str, dict] = (
    json.loads(_summary_path.read_text()) if _summary_path.is_file()
    else _store.get("h2", {})
)

tasks = [(label, feat_dir, loci, mode)
         for label, feat_dir, loci in models
         for mode in ("pca", "qr")]


def _run_task(t):
    label, feat_dir, loci, mode = t
    manifest = run_hereg(label, feat_dir, mode)
    return label, loci, mode, manifest


with ThreadPoolExecutor(max_workers=4) as ex:
    futures = {ex.submit(_run_task, t): t for t in tasks}
    for fut in as_completed(futures):
        label, loci, mode, manifest = fut.result()
        if label not in results:
            results[label] = {"minp_loci_count": loci}
        results[label][mode] = manifest
        print(f"[pipeline] done: {label}/{mode}  "
              f"hereg={manifest.get('sum_h2')}  reml={manifest.get('reml_sum_h2')}",
              flush=True)

HEREG_RUNS_DIR.mkdir(parents=True, exist_ok=True)
_summary_path.write_text(json.dumps(results, indent=2))
print(f"[pipeline] Summary → {_summary_path}")
print("[pipeline] Run plot_scatter.py to regenerate scatter figures.")
