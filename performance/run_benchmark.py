#!/usr/bin/env python3
"""
End-to-end wall-clock benchmarks for IDP heritability methods.

Two modes (--mode):
  trait_scaling   p ∈ {100, 1000, 10000, 100000}, fixed n ≈ 2158 subjects
  sample_scaling  n ∈ {n/3, n, 3n}, fixed p = 1000 traits

Seven methods:
  HE pipeline – PCA 128PCs       (HE_pipeline_PCA_capped_e2e)
  HE pipeline – PCA EVR0.8       (HE_EVR08_e2e)        --only-evr08 [--use-hereg]
  HE pipeline – QR               (QR_HE_e2e)           --only-qr-he [--use-hereg]
  REML pipeline                  (REML_PCA_capped_e2e)
  tr(P⁻¹G)                      (tr_PinvG_e2e)
  SVD(P⁻¹G)                     (SVD_PinvG_e2e)
  blockwise mvGREML              (blockwise_mvGREML_trG_trP_e2e)

All results are written to / merged into results.json in the same directory.

Internal subprocess worker mode (do not call directly):
  run_benchmark.py --_worker <kind> <json-payload>
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor as _TPE
from functools import partial as _partial
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
_AGENT_ROOT = _HERE.parent.parent          # …/performance/ → …/repo/ → …/AGENT/
for _p in (str(_AGENT_ROOT), str(_HERE.parent)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

PCA_N_COMPONENTS = 128

_DEFAULT_GCTA_BIN = "/data4012/zxie3/gcta/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
_DEFAULT_CCOVAR   = "/data484_4/txia2/gwas_practice/T1_ccovar_discovery"
_DEFAULT_QCOVAR   = "/data484_4/txia2/gwas_practice/T1_qcovar_discovery"
DEFAULT_COVAR_MERGED = Path("/data484_4/txia2/gwas_practice/T1_covar_discovery")
RESULTS_JSON = _HERE / "results.json"


# ---------------------------------------------------------------------------
# Lazy imports (from agent root — not bundled in repo)
# ---------------------------------------------------------------------------

def _imports():
    from model.heritability_regression import he_regression_h2, load_grm_dense_from_prefix
    from run_idp_pipeline_king_standard import KinPreset, _ensure_sample_ids, _kin_preset, _norm_iid
    from run_idp_pipeline_king_hereg import pca_only as _pca_only
    return (he_regression_h2, load_grm_dense_from_prefix,
            KinPreset, _ensure_sample_ids, _kin_preset, _norm_iid, _pca_only)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _peak_rss_mb() -> float:
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 3)


def _default_he_workers() -> int:
    return min(8, max(2, os.cpu_count() or 4))


def _checkpoint(t0: float, checkpoints: list, pct: float, label: str) -> None:
    checkpoints.append({"pct": round(pct, 1), "label": label,
                        "elapsed_sec": round(time.perf_counter() - t0, 3)})


def load_covariates_merged(path: str) -> "pd.DataFrame":
    _, _, _, _, _, _norm_iid, _ = _imports()
    df = pd.read_csv(path, sep=r"\s+")
    df["IID_norm"] = df["IID"].astype(str).apply(lambda x: _norm_iid(x))
    id_cols = {"FID", "IID", "IID_norm"}
    cov_cols = [c for c in df.columns if c not in id_cols]
    return df[["IID_norm"] + cov_cols]


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def load_grm_ids(grm_prefix: str):
    _, load_grm_dense_from_prefix, *_ = _imports()
    return load_grm_dense_from_prefix(grm_prefix)


def load_Y_from_feature_dir(feature_dir: Path, df_id: "pd.DataFrame",
                             sample_ids: set, p: int) -> np.ndarray:
    _, _, _, _, _, _norm_iid, _ = _imports()
    n = len(df_id)
    iid_series = df_id["IID"].map(_norm_iid)
    id_to_row: dict[str, int] = dict(zip(iid_series, range(n)))
    Y = np.full((n, p), np.nan, dtype=np.float64)

    def _read_col(j: int) -> None:
        fp = feature_dir / f"Feature_{j}.csv"
        if not fp.is_file():
            fp = feature_dir / f"Feature_{j}"
        if not fp.is_file():
            raise FileNotFoundError(feature_dir / f"Feature_{j}[.csv]")
        df = pd.read_csv(fp, sep=r"\s+")
        pheno_cols = [c for c in df.columns if c not in ("FID", "IID")]
        col = pheno_cols[0]
        iid_norm = df["IID"].map(_norm_iid)
        keep = iid_norm.isin(sample_ids) & iid_norm.isin(id_to_row)
        if keep.any():
            row_idx = iid_norm[keep].map(id_to_row).to_numpy(dtype=int)
            Y[row_idx, j] = pd.to_numeric(df.loc[keep, col], errors="coerce").to_numpy(dtype=np.float64)

    with _TPE(max_workers=min(16, max(2, os.cpu_count() or 4))) as ex:
        list(ex.map(_read_col, range(p)))
    nan_mask = np.isnan(Y)
    if nan_mask.any():
        col_means = np.nanmean(Y, axis=0)
        col_means[np.isnan(col_means)] = 0.0
        Y[np.where(nan_mask)] = col_means[np.where(nan_mask)[1]]
    return Y


def align_Y_to_covariates(Y: np.ndarray, df_id: "pd.DataFrame", sample_ids: set):
    _, _, _, _, _, _norm_iid, _ = _imports()
    covar_df = load_covariates_merged(str(DEFAULT_COVAR_MERGED))
    covar_df = covar_df[covar_df["IID_norm"].isin(sample_ids)].dropna(axis=0, how="any")
    iid_series = df_id["IID"].map(_norm_iid)
    id_to_row: dict[str, int] = dict(zip(iid_series, range(len(iid_series))))
    aligned_iids = [iid for iid in covar_df["IID_norm"].astype(str).tolist() if iid in id_to_row]
    if not aligned_iids:
        raise ValueError("No overlap between synthetic phenotype IIDs and covariates.")
    covar_df = covar_df.set_index("IID_norm").loc[aligned_iids].reset_index()
    indices = [id_to_row[iid] for iid in aligned_iids]
    return Y[indices, :].astype(np.float64, copy=False), aligned_iids, covar_df, indices


def standardize_columns(Y: np.ndarray) -> np.ndarray:
    out = Y.astype(np.float64, copy=True)
    out -= out.mean(axis=0)
    stds = out.std(axis=0)
    stds[stds == 0] = 1.0
    out /= stds
    return out


def residualize_raw_traits(Y: np.ndarray, covar_df: "pd.DataFrame") -> np.ndarray:
    cov_cols = [c for c in covar_df.columns if c != "IID_norm"]
    X_cov = covar_df[cov_cols].values.astype(np.float64)
    X_design = np.hstack([np.ones((X_cov.shape[0], 1), dtype=np.float64), X_cov])
    beta, *_ = np.linalg.lstsq(X_design, Y.astype(np.float64, copy=False), rcond=None)
    return Y - X_design @ beta


def build_G_P_trait(Y: np.ndarray, B: np.ndarray):
    n, p = Y.shape
    BY = B @ Y
    G = Y.T @ BY
    P = (Y.T @ Y) / float(n - 1)
    ridge = float(max(1e-12, np.trace(P) / max(p, 1) * 1e-6))
    P = P + ridge * np.eye(p, dtype=np.float64)
    return G, P


# ---------------------------------------------------------------------------
# PKP helpers
# ---------------------------------------------------------------------------

def _project_grm_pkp(GRM: np.ndarray, X_design: np.ndarray) -> np.ndarray:
    Q, _ = np.linalg.qr(X_design, mode="reduced")
    QtK = Q.T @ GRM
    PK = GRM - Q @ QtK
    return PK - (PK @ Q) @ Q.T


def _build_x_design_for_pkp(covar_df: "pd.DataFrame", eids_sub: list) -> np.ndarray:
    id_col = "IID_norm"
    sub_df = covar_df.set_index(id_col).reindex([str(e) for e in eids_sub]).reset_index()
    cov_cols = [c for c in sub_df.columns if c != id_col]
    X_cov = sub_df[cov_cols].values.astype(np.float64)
    return np.hstack([np.ones((len(eids_sub), 1)), X_cov])


# ---------------------------------------------------------------------------
# HE regression helpers
# ---------------------------------------------------------------------------

_he_fork_scores_e2e = None
_he_fork_grm_e2e = None


def _he_fork_worker_e2e(j: int):
    he_regression_h2, *_ = _imports()
    return he_regression_h2(_he_fork_scores_e2e[:, j], _he_fork_grm_e2e)


def _run_he_parallel(scores: np.ndarray, GRM: np.ndarray, n_he: int, workers: int) -> list:
    import multiprocessing as _mp
    global _he_fork_scores_e2e, _he_fork_grm_e2e
    w = _default_he_workers() if workers < 0 else workers
    if w > 1:
        _he_fork_scores_e2e = scores
        _he_fork_grm_e2e = GRM
        with _mp.get_context("fork").Pool(processes=w) as pool:
            raw = pool.map(_he_fork_worker_e2e, range(n_he))
        _he_fork_scores_e2e = _he_fork_grm_e2e = None
    else:
        he_regression_h2, *_ = _imports()
        raw = [he_regression_h2(scores[:, j], GRM) for j in range(n_he)]
    return [float(v) for v in raw if v is not None and np.isfinite(v)]


def _run_gcta_hereg_feature(i: int, *, features_dir: Path, hereg_out_dir: Path,
                              grm_prefix: str, gcta_bin: str, gcta_threads: int,
                              ccovar_path: str, qcovar_path: str) -> Optional[float]:
    cmd = [gcta_bin, "--HEreg", "--grm", grm_prefix,
           "--pheno", str(features_dir / f"Feature_{i}.csv"),
           "--covar", ccovar_path, "--qcovar", qcovar_path,
           "--thread-num", str(gcta_threads),
           "--out", str(hereg_out_dir / f"Feature_{i}")]
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out_file = hereg_out_dir / f"Feature_{i}.HEreg"
    if not out_file.is_file():
        return None
    for line in out_file.read_text().splitlines():
        parts = line.strip().split()
        if len(parts) >= 2 and parts[0] == "V(G)/Vp":
            try:
                return float(parts[1])
            except (ValueError, IndexError):
                return None
    return None


def _run_hereg_parallel(scores: np.ndarray, eids_sub: list, n_feats: int,
                         grm_prefix: str, gcta_bin: str, hereg_parallel_jobs: int,
                         gcta_threads: int, ccovar_path: str, qcovar_path: str) -> list:
    n_jobs = _default_he_workers() if hereg_parallel_jobs < 1 else max(1, hereg_parallel_jobs)
    threads_per = max(1, gcta_threads // n_jobs)
    tmp_dir = Path(tempfile.mkdtemp(prefix="hereg_e2e_", dir=str(_HERE)))
    feats_dir = tmp_dir / "features"
    hout_dir = tmp_dir / "hereg_out"
    feats_dir.mkdir(parents=True)
    hout_dir.mkdir(parents=True)
    _, _, _, _, _, _norm_iid, _ = _imports()
    try:
        def _write_feat(i: int) -> None:
            out_df = pd.DataFrame({"FID": eids_sub, "IID": eids_sub, str(i): scores[:, i]})
            out_df.to_csv(feats_dir / f"Feature_{i}.csv", sep=" ", index=False)
        with _TPE(max_workers=min(8, n_feats)) as ex:
            list(ex.map(_write_feat, range(n_feats)))
        worker = _partial(_run_gcta_hereg_feature,
                          features_dir=feats_dir, hereg_out_dir=hout_dir,
                          grm_prefix=grm_prefix, gcta_bin=gcta_bin,
                          gcta_threads=threads_per, ccovar_path=ccovar_path,
                          qcovar_path=qcovar_path)
        with _TPE(max_workers=n_jobs) as ex:
            h2_results = list(ex.map(worker, range(n_feats)))
    finally:
        import shutil
        shutil.rmtree(tmp_dir, ignore_errors=True)
    return [float(v) for v in h2_results if v is not None and np.isfinite(v)]


# ---------------------------------------------------------------------------
# Benchmark functions
# ---------------------------------------------------------------------------

def bench_he_pca_capped_e2e(*, preset, feature_dir: Path, p: int,
                              n_pca: int = PCA_N_COMPONENTS,
                              he_parallel_workers: int = -1,
                              use_pkp: bool = False, use_hereg: bool = False,
                              gcta_bin: str = _DEFAULT_GCTA_BIN,
                              hereg_parallel_jobs: int = -1, gcta_threads: int = 8,
                              ccovar_path: str = _DEFAULT_CCOVAR,
                              qcovar_path: str = _DEFAULT_QCOVAR,
                              verbose_progress: bool = True) -> dict:
    _, _, _, _ensure_sample_ids, _, _norm_iid, _pca_only = _imports()
    if not DEFAULT_COVAR_MERGED.is_file():
        raise FileNotFoundError(DEFAULT_COVAR_MERGED)
    sample_ids_path = Path(_ensure_sample_ids(preset, []))
    sample_ids = {_norm_iid(l) for l in open(sample_ids_path) if l.strip()} - {""}
    t0 = time.perf_counter(); r0 = resource.getrusage(resource.RUSAGE_SELF)
    checkpoints: list = []
    GRM, df_id = load_grm_ids(preset.grm_prefix)
    _checkpoint(t0, checkpoints, 15.0, "load_grm")
    if verbose_progress: print(f"[HE p={p}] 15% load_grm  elapsed={checkpoints[-1]['elapsed_sec']}s", flush=True)
    Y = load_Y_from_feature_dir(feature_dir, df_id, sample_ids, p)
    _checkpoint(t0, checkpoints, 35.0, "load_feature_csv_to_Y")
    if verbose_progress: print(f"[HE p={p}] 35% load_Y  elapsed={checkpoints[-1]['elapsed_sec']}s", flush=True)
    Ym, eids, covar_df, grm_indices = align_Y_to_covariates(Y, df_id, sample_ids)
    _checkpoint(t0, checkpoints, 50.0, "covariate_merge_ready")
    pc_scores, eids_sub, evr = _pca_only(Ym, eids, n_components=min(int(n_pca), PCA_N_COMPONENTS))
    _checkpoint(t0, checkpoints, 78.0, "pca_done")
    if verbose_progress: print(f"[HE p={p}] 78% PCA done  elapsed={checkpoints[-1]['elapsed_sec']}s", flush=True)
    eids_to_grm = {eids[k]: grm_indices[k] for k in range(len(eids))}
    indices = [eids_to_grm[e] for e in eids_sub if e in eids_to_grm]
    GRM_sub = GRM[np.ix_(indices, indices)]
    if use_pkp:
        GRM_sub = _project_grm_pkp(GRM_sub, _build_x_design_for_pkp(covar_df, eids_sub))
    n_he = min(int(n_pca), int(p), pc_scores.shape[1])
    if use_hereg:
        h2_list = _run_hereg_parallel(pc_scores, eids_sub, n_he, preset.grm_prefix,
                                       gcta_bin, hereg_parallel_jobs, gcta_threads,
                                       ccovar_path, qcovar_path)
    else:
        h2_list = _run_he_parallel(pc_scores, GRM_sub, n_he, he_parallel_workers)
    _checkpoint(t0, checkpoints, 100.0, "he_on_pcs_done")
    if verbose_progress: print(f"[HE p={p}] 100% done  wall={time.perf_counter()-t0:.1f}s", flush=True)
    wall = time.perf_counter() - t0; r1 = resource.getrusage(resource.RUSAGE_SELF)
    return {"wall_sec": round(wall, 6), "cpu_user_sec": round(r1.ru_utime-r0.ru_utime, 6),
            "cpu_sys_sec": round(r1.ru_stime-r0.ru_stime, 6), "max_rss_mb": _peak_rss_mb(),
            "n_pca_requested": min(int(n_pca), PCA_N_COMPONENTS), "n_he_regressions": n_he,
            "mean_h2_sum": float(np.sum(h2_list)) if h2_list else None, "n_h2": len(h2_list),
            "pc_evr_head": [float(x) for x in evr[:5]], "timing_checkpoints_pct": checkpoints}


def bench_he_evr08_e2e(*, preset, feature_dir: Path, p: int, evr_threshold: float = 0.8,
                        he_parallel_workers: int = -1, use_pkp: bool = False,
                        use_hereg: bool = False, gcta_bin: str = _DEFAULT_GCTA_BIN,
                        hereg_parallel_jobs: int = -1, gcta_threads: int = 8,
                        ccovar_path: str = _DEFAULT_CCOVAR, qcovar_path: str = _DEFAULT_QCOVAR,
                        verbose_progress: bool = True) -> dict:
    _, _, _, _ensure_sample_ids, _, _norm_iid, _pca_only = _imports()
    if not DEFAULT_COVAR_MERGED.is_file():
        raise FileNotFoundError(DEFAULT_COVAR_MERGED)
    sample_ids_path = Path(_ensure_sample_ids(preset, []))
    sample_ids = {_norm_iid(l) for l in open(sample_ids_path) if l.strip()} - {""}
    t0 = time.perf_counter(); r0 = resource.getrusage(resource.RUSAGE_SELF)
    checkpoints: list = []
    GRM, df_id = load_grm_ids(preset.grm_prefix)
    _checkpoint(t0, checkpoints, 15.0, "load_grm")
    if verbose_progress: print(f"[HE-EVR08 p={p}] 15% load_grm  elapsed={checkpoints[-1]['elapsed_sec']}s", flush=True)
    Y = load_Y_from_feature_dir(feature_dir, df_id, sample_ids, p)
    _checkpoint(t0, checkpoints, 35.0, "load_feature_csv_to_Y")
    Ym, eids, covar_df, grm_indices = align_Y_to_covariates(Y, df_id, sample_ids)
    _checkpoint(t0, checkpoints, 50.0, "covariate_merge_ready")
    n_full = min(Ym.shape[0] - 1, p)
    pc_scores, eids_sub, evr = _pca_only(Ym, eids, n_components=n_full)
    cumevr = np.cumsum(evr)
    k_evr = min(int(np.searchsorted(cumevr, evr_threshold)) + 1, pc_scores.shape[1])
    _checkpoint(t0, checkpoints, 78.0, "pca_done")
    if verbose_progress:
        print(f"[HE-EVR08 p={p}] 78% PCA  n_full={n_full}  k_evr={k_evr}  "
              f"cumEVR={float(cumevr[k_evr-1]):.3f}  elapsed={checkpoints[-1]['elapsed_sec']}s", flush=True)
    eids_to_grm = {eids[k]: grm_indices[k] for k in range(len(eids))}
    indices = [eids_to_grm[e] for e in eids_sub if e in eids_to_grm]
    GRM_sub = GRM[np.ix_(indices, indices)]
    if use_pkp:
        GRM_sub = _project_grm_pkp(GRM_sub, _build_x_design_for_pkp(covar_df, eids_sub))
    if use_hereg:
        h2_list = _run_hereg_parallel(pc_scores, eids_sub, k_evr, preset.grm_prefix,
                                       gcta_bin, hereg_parallel_jobs, gcta_threads,
                                       ccovar_path, qcovar_path)
    else:
        h2_list = _run_he_parallel(pc_scores, GRM_sub, k_evr, he_parallel_workers)
    _checkpoint(t0, checkpoints, 100.0, "he_on_pcs_done")
    if verbose_progress: print(f"[HE-EVR08 p={p}] 100% done  wall={time.perf_counter()-t0:.1f}s", flush=True)
    wall = time.perf_counter() - t0; r1 = resource.getrusage(resource.RUSAGE_SELF)
    return {"wall_sec": round(wall, 6), "cpu_user_sec": round(r1.ru_utime-r0.ru_utime, 6),
            "cpu_sys_sec": round(r1.ru_stime-r0.ru_stime, 6), "max_rss_mb": _peak_rss_mb(),
            "n_pca_full": n_full, "n_he_regressions": k_evr, "evr_threshold": evr_threshold,
            "cum_evr_at_cutoff": float(cumevr[k_evr-1]) if len(cumevr) >= k_evr else None,
            "sum_h2": float(np.sum(h2_list)) if h2_list else None, "n_h2": len(h2_list),
            "timing_checkpoints_pct": checkpoints}


def bench_qr_he_e2e(*, preset, feature_dir: Path, p: int,
                     he_parallel_workers: int = -1, use_pkp: bool = False,
                     use_hereg: bool = False, gcta_bin: str = _DEFAULT_GCTA_BIN,
                     hereg_parallel_jobs: int = -1, gcta_threads: int = 8,
                     ccovar_path: str = _DEFAULT_CCOVAR, qcovar_path: str = _DEFAULT_QCOVAR,
                     verbose_progress: bool = True) -> dict:
    _, _, _, _ensure_sample_ids, _, _norm_iid, _ = _imports()
    if not DEFAULT_COVAR_MERGED.is_file():
        raise FileNotFoundError(DEFAULT_COVAR_MERGED)
    sample_ids_path = Path(_ensure_sample_ids(preset, []))
    sample_ids = {_norm_iid(l) for l in open(sample_ids_path) if l.strip()} - {""}
    t0 = time.perf_counter(); r0 = resource.getrusage(resource.RUSAGE_SELF)
    checkpoints: list = []
    GRM, df_id = load_grm_ids(preset.grm_prefix)
    _checkpoint(t0, checkpoints, 15.0, "load_grm")
    Y = load_Y_from_feature_dir(feature_dir, df_id, sample_ids, p)
    _checkpoint(t0, checkpoints, 35.0, "load_feature_csv_to_Y")
    Ym, eids, covar_df, grm_indices = align_Y_to_covariates(Y, df_id, sample_ids)
    _checkpoint(t0, checkpoints, 50.0, "covariate_merge_ready")
    Q, _ = np.linalg.qr(Ym, mode="reduced")
    k = Q.shape[1]
    _checkpoint(t0, checkpoints, 78.0, "qr_done")
    if verbose_progress: print(f"[QR-HE p={p}] 78% QR done  k={k}  elapsed={checkpoints[-1]['elapsed_sec']}s", flush=True)
    GRM_sub = GRM[np.ix_(grm_indices, grm_indices)]
    if use_pkp:
        GRM_sub = _project_grm_pkp(GRM_sub, _build_x_design_for_pkp(covar_df, eids))
    if use_hereg:
        h2_list = _run_hereg_parallel(Q, eids, k, preset.grm_prefix,
                                       gcta_bin, hereg_parallel_jobs, gcta_threads,
                                       ccovar_path, qcovar_path)
    else:
        h2_list = _run_he_parallel(Q, GRM_sub, k, he_parallel_workers)
    _checkpoint(t0, checkpoints, 100.0, "he_on_qr_done")
    if verbose_progress: print(f"[QR-HE p={p}] 100% done  wall={time.perf_counter()-t0:.1f}s", flush=True)
    wall = time.perf_counter() - t0; r1 = resource.getrusage(resource.RUSAGE_SELF)
    return {"wall_sec": round(wall, 6), "cpu_user_sec": round(r1.ru_utime-r0.ru_utime, 6),
            "cpu_sys_sec": round(r1.ru_stime-r0.ru_stime, 6), "max_rss_mb": _peak_rss_mb(),
            "n_qr_components": k, "n_he_regressions": k,
            "sum_h2": float(np.sum(h2_list)) if h2_list else None, "n_h2": len(h2_list),
            "timing_checkpoints_pct": checkpoints}


def bench_dense_ops_e2e(*, preset, feature_dir: Path, p: int, kind: str,
                         n_pca: int = PCA_N_COMPONENTS,
                         verbose_progress: bool = True) -> dict:
    _, _, _, _ensure_sample_ids, _, _norm_iid, _ = _imports()
    sample_ids_path = Path(_ensure_sample_ids(preset, []))
    sample_ids = {_norm_iid(l) for l in open(sample_ids_path) if l.strip()} - {""}
    t0 = time.perf_counter(); r0 = resource.getrusage(resource.RUSAGE_SELF)
    checkpoints: list = []; tag = f"dense-{kind} p={p}"
    GRM, df_id = load_grm_ids(preset.grm_prefix)
    _checkpoint(t0, checkpoints, 20.0, "load_grm")
    Y_raw = load_Y_from_feature_dir(feature_dir, df_id, sample_ids, p)
    _checkpoint(t0, checkpoints, 45.0, "load_feature_csv")
    Ym, _eids, covar_df, indices = align_Y_to_covariates(Y_raw, df_id, sample_ids)
    Y_resid = residualize_raw_traits(Ym, covar_df)
    del Y_raw, Ym
    _checkpoint(t0, checkpoints, 55.0, "residualize_raw_traits")
    Y = standardize_columns(Y_resid); del Y_resid
    GRM_sub = GRM[np.ix_(indices, indices)]
    n = Y.shape[0]
    B = GRM_sub - np.eye(n, dtype=np.float64)
    G, P = build_G_P_trait(Y, B)
    _checkpoint(t0, checkpoints, 70.0, "build_G_P")
    extra: dict = {}
    if kind == "trace":
        L = np.linalg.cholesky(P); tmp = np.linalg.solve(L, G); tmp = np.linalg.solve(L.T, tmp)
        tr = float(np.trace(tmp)); extra["trace_PinvG"] = tr; extra["trace_over_p"] = tr / p
    elif kind == "svd":
        L = np.linalg.cholesky(P); tmp = np.linalg.solve(L, G); tmp = np.linalg.solve(L.T, tmp)
        _u, s, _vh = np.linalg.svd(tmp, full_matrices=False)
        extra["sv_top1"] = float(s[0]) if len(s) else float("nan")
        extra["sv_top10_sum"] = float(np.sum(s[:10])) if len(s) else 0.0
    elif kind == "mvgreml":
        tr_g = float(np.trace(G)); tr_p = float(np.trace(P))
        extra.update({"trG": tr_g, "trP": tr_p, "trG_over_trP": tr_g / tr_p if tr_p else float("nan")})
    else:
        raise ValueError(kind)
    _checkpoint(t0, checkpoints, 100.0, f"solve_{kind}")
    if verbose_progress: print(f"[{tag}] 100% done  wall={time.perf_counter()-t0:.1f}s", flush=True)
    wall = time.perf_counter() - t0; r1 = resource.getrusage(resource.RUSAGE_SELF)
    return {"wall_sec": round(wall, 6), "cpu_user_sec": round(r1.ru_utime-r0.ru_utime, 6),
            "cpu_sys_sec": round(r1.ru_stime-r0.ru_stime, 6), "max_rss_mb": _peak_rss_mb(),
            "n_samples": int(n), "residualized_raw_traits": True,
            "timing_checkpoints_pct": checkpoints, **extra}


def bench_blockwise_mvgreml_e2e(*, preset, feature_dir: Path, p: int,
                                  block_size: int = 1024,
                                  verbose_progress: bool = True) -> dict:
    _, _, _, _ensure_sample_ids, _, _norm_iid, _ = _imports()
    sample_ids_path = Path(_ensure_sample_ids(preset, []))
    sample_ids = {_norm_iid(l) for l in open(sample_ids_path) if l.strip()} - {""}
    t0 = time.perf_counter(); r0 = resource.getrusage(resource.RUSAGE_SELF)
    checkpoints: list = []; tag = f"blockwise-mvgreml p={p}"
    GRM, df_id = load_grm_ids(preset.grm_prefix)
    _checkpoint(t0, checkpoints, 20.0, "load_grm")
    Y_raw = load_Y_from_feature_dir(feature_dir, df_id, sample_ids, p)
    _checkpoint(t0, checkpoints, 45.0, "load_feature_csv")
    Ym, _eids, covar_df, indices = align_Y_to_covariates(Y_raw, df_id, sample_ids)
    Y_resid = residualize_raw_traits(Ym, covar_df); del Y_raw, Ym
    _checkpoint(t0, checkpoints, 55.0, "residualize_raw_traits")
    Y = standardize_columns(Y_resid); del Y_resid
    GRM_sub = GRM[np.ix_(indices, indices)]; del GRM
    n = Y.shape[0]; B = GRM_sub - np.eye(n, dtype=np.float64); del GRM_sub
    _checkpoint(t0, checkpoints, 62.0, "build_B")
    trace_g = 0.0; trace_p_raw = 0.0; n_blocks = 0
    for start in range(0, p, block_size):
        Yb = Y[:, start:min(p, start + block_size)]
        trace_g += float(np.sum(Yb * (B @ Yb)))
        trace_p_raw += float(np.sum(Yb * Yb) / float(n - 1))
        n_blocks += 1
    ridge = float(max(1e-12, trace_p_raw / p * 1e-6))
    trace_p = trace_p_raw + p * ridge
    _checkpoint(t0, checkpoints, 100.0, "blockwise_traces_done")
    if verbose_progress: print(f"[{tag}] 100% done  wall={time.perf_counter()-t0:.1f}s", flush=True)
    wall = time.perf_counter() - t0; r1 = resource.getrusage(resource.RUSAGE_SELF)
    return {"wall_sec": round(wall, 6), "cpu_user_sec": round(r1.ru_utime-r0.ru_utime, 6),
            "cpu_sys_sec": round(r1.ru_stime-r0.ru_stime, 6), "max_rss_mb": _peak_rss_mb(),
            "n_samples": int(n), "block_size": block_size, "n_blocks": n_blocks,
            "residualized_raw_traits": True, "trG": trace_g, "trP": trace_p,
            "trG_over_trP": trace_g / trace_p if trace_p else float("nan"),
            "timing_checkpoints_pct": checkpoints}


def bench_reml_pca_capped_e2e(*, preset, feature_dir: Path, p: int,
                                gcta_bin: str, gcta_threads: int, n_pca: int,
                                tmp_root: Path, verbose_progress: bool = True) -> dict:
    from run_idp_pipeline_king_reml import _read_h2_from_hsq, _run_reml_for_feature
    _, _, _, _ensure_sample_ids, _, _norm_iid, _pca_only = _imports()
    if not DEFAULT_COVAR_MERGED.is_file():
        raise FileNotFoundError(DEFAULT_COVAR_MERGED)
    sample_ids_path = Path(_ensure_sample_ids(preset, []))
    sample_ids = {_norm_iid(l) for l in open(sample_ids_path) if l.strip()} - {""}
    run_dir = tmp_root / f"reml_e2e_p{p}"; run_dir.mkdir(parents=True, exist_ok=True)
    features_dir = run_dir / "features"; reml_out_dir = run_dir / "reml_out"
    checkpoints: list = []; t0 = time.perf_counter(); r0 = resource.getrusage(resource.RUSAGE_SELF)
    GRM, df_id = load_grm_ids(preset.grm_prefix)
    _checkpoint(t0, checkpoints, 15.0, "load_grm")
    Y = load_Y_from_feature_dir(feature_dir, df_id, sample_ids, p)
    _checkpoint(t0, checkpoints, 35.0, "load_feature_csv_to_Y")
    Ym, eids, _covar_df, _ = align_Y_to_covariates(Y, df_id, sample_ids)
    _checkpoint(t0, checkpoints, 45.0, "covariate_alignment_ready")
    pc_scores, eids_sub, evr = _pca_only(Ym, eids, n_components=min(int(n_pca), PCA_N_COMPONENTS))
    _checkpoint(t0, checkpoints, 75.0, "pca_done")
    features_dir.mkdir(parents=True, exist_ok=True)
    n_save = pc_scores.shape[1]
    sfw = min(8, max(2, (os.cpu_count() or 4) // 2))
    def _write_f(i): pd.DataFrame({"FID": eids_sub, "IID": eids_sub, str(i): pc_scores[:, i]}).to_csv(features_dir / f"Feature_{i}.csv", sep=" ", index=False)
    with _TPE(max_workers=sfw) as ex: list(ex.map(_write_f, range(n_save)))
    with open(features_dir / "pca_explained_variance_ratio.json", "w") as f: json.dump(evr, f)
    _checkpoint(t0, checkpoints, 82.0, "pc_features_saved")
    reml_out_dir.mkdir(parents=True, exist_ok=True)
    max_features = min(int(n_pca), int(p), n_save)
    n_jobs = min(8, max_features); threads_per = max(1, gcta_threads // n_jobs)
    def _run_one(i):
        code, log = _run_reml_for_feature(gcta_bin=gcta_bin, grm_prefix=preset.grm_prefix,
            pheno_file=features_dir / f"Feature_{i}.csv",
            out_prefix=reml_out_dir / f"Feature_{i}", threads=threads_per)
        open(reml_out_dir / f"Feature_{i}.log", "w").write(log)
        h2 = _read_h2_from_hsq(reml_out_dir / f"Feature_{i}.hsq")
        ev = float(evr[i]) if i < len(evr) else 0.0
        return i, code, h2, ev
    rows_dict: dict = {}; returncodes: list = []; log_tail = ""
    with _TPE(max_workers=n_jobs) as ex:
        for i, code, h2, ev in ex.map(_run_one, range(max_features)):
            returncodes.append(int(code)); rows_dict[i] = {"feature_index": i, "return_code": code, "h2": h2, "explained_variance_ratio": ev}
            if verbose_progress: print(f"[REML p={p}] Feature_{i}: code={code}, h2={h2}", flush=True)
    reml_df = pd.DataFrame([rows_dict[i] for i in range(max_features)])
    reml_df.to_csv(run_dir / "reml_h2_summary.csv", index=False)
    valid = reml_df["h2"].notna() if len(reml_df) else pd.Series(dtype=bool)
    mean_h2 = None
    if len(reml_df) and valid.any():
        w = reml_df.loc[valid, "explained_variance_ratio"].astype(float).values
        h2v = reml_df.loc[valid, "h2"].astype(float).values
        w_sum = float(w.sum())
        if w_sum > 0: mean_h2 = float((h2v * (w / w_sum)).sum())
    _checkpoint(t0, checkpoints, 100.0, "gcta_reml_done")
    wall = time.perf_counter() - t0; r1 = resource.getrusage(resource.RUSAGE_SELF)
    ok = all(c == 0 for c in returncodes)
    return {"wall_sec": round(wall, 6), "cpu_user_sec": round(r1.ru_utime-r0.ru_utime, 6),
            "cpu_sys_sec": round(r1.ru_stime-r0.ru_stime, 6), "max_rss_mb": _peak_rss_mb(),
            "returncode": 0 if ok else 1, "ok": ok, "n_reml_features": max_features,
            "mean_h2_reml": mean_h2, "pc_evr_head": [float(x) for x in evr[:5]],
            "timing_checkpoints_pct": checkpoints}


# ---------------------------------------------------------------------------
# Subprocess isolation helper + worker dispatch
# ---------------------------------------------------------------------------

def _run_bench_isolated(kind: str, payload: dict, *, label: str) -> dict:
    proc = subprocess.run(
        [sys.executable, str(_HERE / "run_benchmark.py"), "--_worker", kind, json.dumps(payload)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=str(_HERE))
    if proc.returncode != 0:
        raise RuntimeError(f"{label}: worker exit {proc.returncode}\n{proc.stderr[-4000:]}")
    line = proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else ""
    try:
        return json.loads(line)
    except json.JSONDecodeError as e:
        raise RuntimeError(f"{label}: bad JSON: {e}\nstdout={proc.stdout[-2000:]!r}") from e


def _worker_main(kind: str, payload: dict) -> None:
    _, _, _, _, _kin_preset, _norm_iid, _ = _imports()
    preset = _kin_preset(cohort=payload["cohort"], kin=payload["kin"])
    feature_dir = Path(payload["feature_dir"])
    p = int(payload["p"])
    verbose = bool(payload.get("verbose", False))
    use_hereg = bool(payload.get("use_hereg", False))

    if kind == "HE":
        out = bench_he_pca_capped_e2e(preset=preset, feature_dir=feature_dir, p=p,
            n_pca=int(payload.get("n_pca", PCA_N_COMPONENTS)),
            he_parallel_workers=int(payload.get("he_parallel_workers", -1)),
            use_hereg=use_hereg,
            gcta_bin=payload.get("gcta_bin", _DEFAULT_GCTA_BIN),
            hereg_parallel_jobs=int(payload.get("hereg_parallel_jobs", -1)),
            gcta_threads=int(payload.get("gcta_threads", 8)),
            ccovar_path=payload.get("ccovar_path", _DEFAULT_CCOVAR),
            qcovar_path=payload.get("qcovar_path", _DEFAULT_QCOVAR),
            verbose_progress=verbose)
    elif kind == "he_evr08":
        out = bench_he_evr08_e2e(preset=preset, feature_dir=feature_dir, p=p,
            he_parallel_workers=int(payload.get("he_parallel_workers", -1)),
            use_hereg=use_hereg,
            gcta_bin=payload.get("gcta_bin", _DEFAULT_GCTA_BIN),
            hereg_parallel_jobs=int(payload.get("hereg_parallel_jobs", -1)),
            gcta_threads=int(payload.get("gcta_threads", 8)),
            ccovar_path=payload.get("ccovar_path", _DEFAULT_CCOVAR),
            qcovar_path=payload.get("qcovar_path", _DEFAULT_QCOVAR),
            verbose_progress=verbose)
    elif kind == "qr_he":
        out = bench_qr_he_e2e(preset=preset, feature_dir=feature_dir, p=p,
            he_parallel_workers=int(payload.get("he_parallel_workers", -1)),
            use_hereg=use_hereg,
            gcta_bin=payload.get("gcta_bin", _DEFAULT_GCTA_BIN),
            hereg_parallel_jobs=int(payload.get("hereg_parallel_jobs", -1)),
            gcta_threads=int(payload.get("gcta_threads", 8)),
            ccovar_path=payload.get("ccovar_path", _DEFAULT_CCOVAR),
            qcovar_path=payload.get("qcovar_path", _DEFAULT_QCOVAR),
            verbose_progress=verbose)
    elif kind in ("dense_trace", "dense_svd", "dense_mvgreml"):
        out = bench_dense_ops_e2e(preset=preset, feature_dir=feature_dir, p=p,
            kind=payload["dense_kind"],
            n_pca=int(payload.get("n_pca", PCA_N_COMPONENTS)), verbose_progress=verbose)
    elif kind == "blockwise_mvgreml":
        out = bench_blockwise_mvgreml_e2e(preset=preset, feature_dir=feature_dir, p=p,
            block_size=int(payload.get("block_size", 1024)), verbose_progress=verbose)
    elif kind == "REML":
        out = bench_reml_pca_capped_e2e(preset=preset, feature_dir=feature_dir, p=p,
            gcta_bin=payload.get("gcta_bin", _DEFAULT_GCTA_BIN),
            gcta_threads=int(payload.get("gcta_threads", 8)),
            n_pca=int(payload.get("n_pca", PCA_N_COMPONENTS)),
            tmp_root=Path(payload["tmp_root"]), verbose_progress=verbose)
    else:
        raise SystemExit(f"unknown worker kind {kind!r}")
    print(json.dumps(out), flush=True)


# ---------------------------------------------------------------------------
# JSON merge helpers
# ---------------------------------------------------------------------------

def _load_results() -> dict:
    if RESULTS_JSON.is_file():
        return json.loads(RESULTS_JSON.read_text())
    return {"trait_scaling": {"p_values": [], "rows": []},
            "sample_scaling": {"p_fixed": 1000, "rows": []}}


def _save_results(data: dict) -> None:
    RESULTS_JSON.write_text(json.dumps(data, indent=2))


def _merge_trait_rows(existing: list, new_rows: list) -> list:
    by_p: dict = {int(r["p"]): dict(r) for r in existing}
    for r in new_rows:
        p = int(r["p"])
        if p not in by_p:
            by_p[p] = dict(r)
        else:
            for k, v in r.items():
                if k == "p": continue
                if isinstance(v, dict) and v.get("skipped"): continue
                by_p[p][k] = v
    return sorted(by_p.values(), key=lambda x: int(x["p"]))


# ---------------------------------------------------------------------------
# Trait scaling main
# ---------------------------------------------------------------------------

def run_trait_scaling(args) -> None:
    _, _, _, _, _kin_preset, _, _ = _imports()
    traits_root = Path(args.traits_root)
    p_values = [int(x.strip()) for x in args.p_list.split(",") if x.strip()]
    preset = _kin_preset(cohort="discovery", kin="over4p5")
    cohort_s, kin_s = "discovery", "over4p5"
    vp = not args.quiet
    skip_dense = int(args.skip_dense_above_p)
    isolate = not args.no_isolate
    use_hereg = args.use_hereg

    rows_out: list = []
    tmp_root = Path(tempfile.mkdtemp(prefix="e2e_bench_", dir=str(_HERE)))

    try:
        for p in p_values:
            fdir = traits_root / f"traits_{p}"
            if not fdir.is_dir():
                raise FileNotFoundError(fdir)
            if vp:
                print(f"\n=== E2E benchmark p={p} ===", flush=True)

            base_payload = {"cohort": cohort_s, "kin": kin_s, "feature_dir": str(fdir),
                            "p": p, "verbose": vp, "use_hereg": use_hereg,
                            "gcta_bin": args.gcta_bin, "gcta_threads": args.gcta_threads,
                            "hereg_parallel_jobs": args.hereg_parallel_jobs,
                            "ccovar_path": _DEFAULT_CCOVAR, "qcovar_path": _DEFAULT_QCOVAR,
                            "he_parallel_workers": args.he_parallel_workers}

            if args.only_evr08:
                row_evr = (_run_bench_isolated("he_evr08", base_payload, label=f"EVR08 p={p}")
                           if isolate else bench_he_evr08_e2e(preset=preset, feature_dir=fdir, p=p,
                               use_hereg=use_hereg, verbose_progress=vp))
                rows_out.append({"p": p, "feature_dir": str(fdir), "HE_EVR08_e2e": row_evr})
                continue

            if args.only_qr_he:
                row_qr = (_run_bench_isolated("qr_he", base_payload, label=f"QR-HE p={p}")
                          if isolate else bench_qr_he_e2e(preset=preset, feature_dir=fdir, p=p,
                              use_hereg=use_hereg, verbose_progress=vp))
                rows_out.append({"p": p, "feature_dir": str(fdir), "QR_HE_e2e": row_qr})
                continue

            row_he = (_run_bench_isolated("HE", {**base_payload, "n_pca": args.n_pca}, label=f"HE p={p}")
                      if isolate else bench_he_pca_capped_e2e(preset=preset, feature_dir=fdir, p=p,
                          n_pca=args.n_pca, use_hereg=use_hereg, verbose_progress=vp))
            row: dict = {"p": p, "feature_dir": str(fdir), "HE_pipeline_PCA_capped_e2e": row_he}

            def _skipped(why): return {"skipped": True, "reason": why}
            if p <= skip_dense:
                for kind, key, wkind in (("trace", "tr_PinvG_e2e", "dense_trace"),
                                          ("svd", "SVD_PinvG_e2e", "dense_svd")):
                    row[key] = (_run_bench_isolated(wkind, {**base_payload, "dense_kind": kind, "n_pca": args.n_pca}, label=f"dense {kind} p={p}")
                                if isolate else bench_dense_ops_e2e(preset=preset, feature_dir=fdir, p=p, kind=kind, verbose_progress=vp))
            else:
                why = f"p={p} > skip_dense_above_p={skip_dense}"
                row["tr_PinvG_e2e"] = _skipped(why); row["SVD_PinvG_e2e"] = _skipped(why)

            row["blockwise_mvGREML_trG_trP_e2e"] = (
                _run_bench_isolated("blockwise_mvgreml", base_payload, label=f"blockwise p={p}")
                if isolate else bench_blockwise_mvgreml_e2e(preset=preset, feature_dir=fdir, p=p, verbose_progress=vp))

            if args.skip_reml:
                row["REML_PCA_capped_e2e"] = {"skipped": True}
            else:
                row["REML_PCA_capped_e2e"] = (
                    _run_bench_isolated("REML", {**base_payload, "n_pca": args.n_pca, "tmp_root": str(tmp_root)}, label=f"REML p={p}")
                    if isolate else bench_reml_pca_capped_e2e(preset=preset, feature_dir=fdir, p=p,
                        gcta_bin=args.gcta_bin, gcta_threads=args.gcta_threads, n_pca=args.n_pca, tmp_root=tmp_root, verbose_progress=vp))
            rows_out.append(row)
    finally:
        pass  # keep tmp_root for debugging

    store = _load_results()
    store["trait_scaling"]["rows"] = _merge_trait_rows(store["trait_scaling"].get("rows", []), rows_out)
    store["trait_scaling"]["p_values"] = sorted({int(r["p"]) for r in store["trait_scaling"]["rows"]})
    _save_results(store)
    print(f"[trait_scaling] Updated {RESULTS_JSON}", flush=True)


# ---------------------------------------------------------------------------
# Sample scaling main
# ---------------------------------------------------------------------------

_SAMPLE_METHODS = ["HE_pipeline_PCA_capped_e2e", "HE_EVR08_e2e", "QR_HE_e2e",
                   "REML_PCA_capped_e2e", "tr_PinvG_e2e", "SVD_PinvG_e2e",
                   "blockwise_mvGREML_trG_trP_e2e"]


def run_sample_scaling(args) -> None:
    _, _, KinPreset, _, _kin_preset, _norm_iid, _ = _imports()
    AGENT_ROOT = _HERE.parent.parent
    OVER5_PREFIX = AGENT_ROOT / "grm_subset_king" / "king_over5_gcta_discovery"
    OVER45_IDS  = AGENT_ROOT / "grm_subset_king" / "sample_ids_king_over4p5_discovery.txt"
    OVER5_IDS   = AGENT_ROOT / "grm_subset_king" / "sample_ids_king_over5_discovery.txt"
    SAMPLE_DATA_DIR = _HERE / "sample_scaling_data"
    SAMPLE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    TRAITS_DIR = SAMPLE_DATA_DIR / "synthetic_over5_traits_1000"
    GRM_DIR = SAMPLE_DATA_DIR / "subset_grms"
    GRM_DIR.mkdir(parents=True, exist_ok=True)
    P_FIXED = 1000; N_PCA = 128; GCTA_BIN = args.gcta_bin

    def read_ids(path: Path) -> list:
        return [_norm_iid(x.strip()) for x in path.read_text().splitlines() if x.strip()]

    def write_lower_triangle(path: Path, matrix: np.ndarray) -> None:
        vals = []
        for i in range(matrix.shape[0]):
            vals.append(matrix[i, :i+1].astype(np.float32, copy=False))
        np.concatenate(vals).astype(np.float32, copy=False).tofile(path)

    def subset_grm(master_prefix: Path, ids: list, out_prefix: Path) -> None:
        id_out = Path(str(out_prefix) + ".grm.id")
        bin_out = Path(str(out_prefix) + ".grm.bin")
        n_out = Path(str(out_prefix) + ".grm.N.bin")
        if id_out.exists() and bin_out.exists() and n_out.exists():
            return
        df_id = pd.read_csv(str(master_prefix) + ".grm.id", sep=r"\s+", header=None, names=["FID", "IID"])
        master_ids = [_norm_iid(x) for x in df_id["IID"].tolist()]
        idx_map = {iid: i for i, iid in enumerate(master_ids)}
        indices = [idx_map[iid] for iid in ids if iid in idx_map]
        n = len(master_ids); nval = n * (n + 1) // 2
        data = np.fromfile(str(master_prefix) + ".grm.bin", dtype=np.float32, count=nval)
        counts = np.fromfile(str(master_prefix) + ".grm.N.bin", dtype=np.float32, count=nval)
        full = np.zeros((n, n), dtype=np.float32); full_n = np.zeros((n, n), dtype=np.float32)
        k = 0
        for i in range(n):
            span = i + 1
            full[i, :span] = data[k:k+span]; full[:span, i] = data[k:k+span]
            full_n[i, :span] = counts[k:k+span]; full_n[:span, i] = counts[k:k+span]
            k += span
        sub = full[np.ix_(indices, indices)]; sub_n = full_n[np.ix_(indices, indices)]
        df_id.iloc[indices].to_csv(id_out, sep="\t", header=False, index=False)
        write_lower_triangle(bin_out, sub); write_lower_triangle(n_out, sub_n)

    def prepare_traits(ids: list) -> None:
        TRAITS_DIR.mkdir(parents=True, exist_ok=True)
        if (TRAITS_DIR / f"Feature_{P_FIXED-1}.csv").exists():
            return
        fid = [int(x) if str(x).isdigit() else x for x in ids]
        for j in range(P_FIXED):
            rng = np.random.default_rng(1000003 + j)
            pd.DataFrame({"FID": fid, "IID": fid, str(j): rng.standard_normal(len(ids))}).to_csv(
                TRAITS_DIR / f"Feature_{j}.csv", sep=" ", index=False)

    def subset_traits(label: str, ids: list) -> Path:
        out_dir = SAMPLE_DATA_DIR / f"traits_{label}"
        if (out_dir / f"Feature_{P_FIXED-1}.csv").exists():
            return out_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        id_set = set(ids)
        for j in range(P_FIXED):
            df = pd.read_csv(TRAITS_DIR / f"Feature_{j}.csv", sep=r"\s+")
            iid_norm = df["IID"].map(_norm_iid)
            df.loc[iid_norm.isin(id_set)].to_csv(out_dir / f"Feature_{j}.csv", sep=" ", index=False)
        return out_dir

    over45 = read_ids(OVER45_IDS); over5 = read_ids(OVER5_IDS)
    extras = [iid for iid in over5 if iid not in set(over45)]
    sample_sets = [
        ("n_small_1x3", over45[:max(1, len(over45)//3)]),
        ("n_mid_current", over45),
        ("n_large_3x", over45 + extras[:max(0, len(over45)*3 - len(over45))]),
    ]
    prepare_traits(read_ids(OVER5_IDS))

    store = _load_results()
    by_label = {r["label"]: r for r in store["sample_scaling"].get("rows", [])}

    for label, ids in sample_sets:
        sample_path = SAMPLE_DATA_DIR / f"sample_ids_{label}.txt"
        sample_path.write_text("\n".join(ids) + "\n")
        grm_prefix = GRM_DIR / f"king_over5_discovery_{label}"
        subset_grm(OVER5_PREFIX, ids, grm_prefix)
        feature_dir = subset_traits(label, ids)
        preset = KinPreset(label=label, keep_file=sample_path, grm_prefix=str(grm_prefix), sample_ids_file=sample_path)
        row = by_label.get(label, {"label": label, "n_target": len(ids),
                                    "sample_ids_file": str(sample_path), "grm_prefix": str(grm_prefix)})
        row["feature_dir"] = str(feature_dir)
        sample_ids_set = {_norm_iid(x) for x in ids}

        for method in _SAMPLE_METHODS:
            if method in row and not row[method].get("failed"):
                continue
            print(f"\n=== sample scaling {label} n={len(ids)} method={method} ===", flush=True)
            try:
                if method == "QR_HE_e2e":
                    row[method] = bench_qr_he_e2e(preset=preset, feature_dir=feature_dir, p=P_FIXED)
                elif method == "HE_pipeline_PCA_capped_e2e":
                    row[method] = bench_he_pca_capped_e2e(preset=preset, feature_dir=feature_dir, p=P_FIXED, n_pca=N_PCA)
                elif method == "HE_EVR08_e2e":
                    row[method] = bench_he_evr08_e2e(preset=preset, feature_dir=feature_dir, p=P_FIXED)
                elif method == "REML_PCA_capped_e2e":
                    tmp_root = SAMPLE_DATA_DIR / "tmp_reml"
                    row[method] = bench_reml_pca_capped_e2e(preset=preset, feature_dir=feature_dir,
                        p=P_FIXED, gcta_bin=GCTA_BIN, gcta_threads=8, n_pca=N_PCA, tmp_root=tmp_root)
                elif method == "tr_PinvG_e2e":
                    row[method] = bench_dense_ops_e2e(preset=preset, feature_dir=feature_dir, p=P_FIXED, kind="trace")
                elif method == "SVD_PinvG_e2e":
                    row[method] = bench_dense_ops_e2e(preset=preset, feature_dir=feature_dir, p=P_FIXED, kind="svd")
                elif method == "blockwise_mvGREML_trG_trP_e2e":
                    row[method] = bench_blockwise_mvgreml_e2e(preset=preset, feature_dir=feature_dir, p=P_FIXED)
            except Exception as e:
                row[method] = {"failed": True, "error": repr(e)}
            row["n_observed"] = int(row.get(method, {}).get("n_samples", len(ids)))
        by_label[label] = row
        store["sample_scaling"]["rows"] = sorted(by_label.values(), key=lambda x: x.get("n_target", 0))
        _save_results(store)

    print(f"[sample_scaling] Updated {RESULTS_JSON}", flush=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    # Internal worker dispatch (subprocess isolation mode)
    if len(sys.argv) >= 4 and sys.argv[1] == "--_worker":
        _worker_main(sys.argv[2], json.loads(sys.argv[3]))
        return

    ap = argparse.ArgumentParser(description="IDP heritability scaling benchmarks.")
    ap.add_argument("--mode", choices=["trait_scaling", "sample_scaling"],
                    default="trait_scaling")

    # Trait scaling options
    ap.add_argument("--traits-root", type=str,
                    default=str(_HERE.parent / "synthetic_phenotypes"))
    ap.add_argument("--p-list", type=str, default="100,1000,10000")
    ap.add_argument("--only-evr08", action="store_true",
                    help="Run only HE pipeline with EVR=0.8 dynamic PCA cutoff.")
    ap.add_argument("--only-qr-he", action="store_true",
                    help="Run only QR+HE benchmark.")
    ap.add_argument("--use-hereg", action="store_true",
                    help="Use GCTA --HEreg (with --covar/--qcovar) for EVR0.8 and QR variants.")
    ap.add_argument("--skip-reml", action="store_true")
    ap.add_argument("--skip-dense-above-p", type=int, default=10000)
    ap.add_argument("--no-isolate", action="store_true",
                    help="Run all benchmarks in one process (faster, but RSS is a shared high-water mark).")
    ap.add_argument("--he-parallel-workers", type=int, default=-1)
    ap.add_argument("--hereg-parallel-jobs", type=int, default=8)
    ap.add_argument("--n-pca", type=int, default=128)

    # Shared options
    ap.add_argument("--gcta-bin", type=str, default=_DEFAULT_GCTA_BIN)
    ap.add_argument("--gcta-threads", type=int, default=8)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    if args.mode == "trait_scaling":
        run_trait_scaling(args)
    else:
        run_sample_scaling(args)


if __name__ == "__main__":
    main()
