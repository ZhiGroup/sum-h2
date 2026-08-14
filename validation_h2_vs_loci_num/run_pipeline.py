#!/usr/bin/env python3
"""
Run GCTA HEreg+REML pipeline on phenotype models and save results.

Model list comes from results.json["kinship_sensitivity"]["rows"], restricted to
labels that have JAGWAS lmm_only loci in results.json["loci"].

If a model has already been run (arena_manifest.json exists in hereg_runs/), it is
skipped. After all models finish, updates hereg_runs/hereg_summary.json.

Uses PCA 128PCs only.
Run plot_scatter.py separately to regenerate the JAGWAS scatter figure.
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
HEREG_RUNS_DIR = HERE / "hereg_runs"
HEREG_SCRIPT   = HERE.parent / "run_sum_h2.py"
PYTHON         = sys.executable


_store = json.loads(RESULTS_JSON.read_text())
loci_by_label = _store.get("loci", {})
loci_rows = _store["kinship_sensitivity"]["rows"]
models = []
for r in loci_rows:
    if r.get("source_row_error") or not r.get("original_feature_dir"):
        continue
    label = r["raw_label"]
    jag = (loci_by_label.get(label) or {}).get("jagwas")
    if jag is None:
        continue
    models.append((label, r["original_feature_dir"], float(jag)))
print(f"[pipeline] {len(models)} models with JAGWAS loci loaded from results.json")


def _safe(label: str) -> str:
    return re.sub(r"[/\\]", "_", label)


def run_hereg(label: str, feature_dir: str) -> dict:
    """PCA 128PCs with --h2-method both (HEreg + REML)."""
    out_root = HEREG_RUNS_DIR / _safe(label) / "pca"
    manifest_path = out_root / "arena_manifest.json"

    if manifest_path.is_file():
        print(f"[pipeline] skip (cached): {label}/pca")
        return json.loads(manifest_path.read_text())

    out_root.mkdir(parents=True, exist_ok=True)
    cmd = [PYTHON, str(HEREG_SCRIPT),
           "--phenotype_csv", feature_dir,
           "--output_root", str(out_root),
           "--dim-reduction", "pca",
           "--n_pca", "128",
           "--h2-method", "both"]

    print(f"[pipeline] running {label}/pca ...")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[pipeline] ERROR {label}/pca:\n{res.stderr[-400:]}", flush=True)
        return {"error": res.stderr[-300:]}

    if manifest_path.is_file():
        return json.loads(manifest_path.read_text())
    return {"error": "manifest not written"}


_summary_path = HEREG_RUNS_DIR / "hereg_summary.json"
results: dict[str, dict] = (
    json.loads(_summary_path.read_text()) if _summary_path.is_file()
    else _store.get("h2", {})
)

tasks = [(label, feat_dir, loci) for label, feat_dir, loci in models]


def _run_task(t):
    label, feat_dir, loci = t
    manifest = run_hereg(label, feat_dir)
    return label, loci, manifest


with ThreadPoolExecutor(max_workers=4) as ex:
    futures = {ex.submit(_run_task, t): t for t in tasks}
    for fut in as_completed(futures):
        label, loci, manifest = fut.result()
        if label not in results:
            results[label] = {}
        results[label]["jagwas_loci_count"] = loci
        results[label]["pca"] = manifest
        print(f"[pipeline] done: {label}/pca  "
              f"hereg={manifest.get('sum_h2')}  reml={manifest.get('reml_sum_h2')}",
              flush=True)

HEREG_RUNS_DIR.mkdir(parents=True, exist_ok=True)
_summary_path.write_text(json.dumps(results, indent=2))
print(f"[pipeline] Summary → {_summary_path}")
print("[pipeline] Run plot_scatter.py to regenerate the JAGWAS scatter figure.")
