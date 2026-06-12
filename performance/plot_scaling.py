#!/usr/bin/env python3
"""
Generate scaling plots from results.json.

Writes two PNG files in the same directory:
  IDP_synthetic_trait_scaling_time_resource_v8_hereg.png
  IDP_synthetic_sample_size_scaling_time_resource_v2_hereg.png

Reads all data from results.json — no external file dependencies.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = Path(__file__).resolve().parent
RESULTS_JSON = _HERE / "results.json"
OUT_TRAIT  = _HERE / "IDP_synthetic_trait_scaling_time_resource_v8_hereg.png"
OUT_SAMPLE = _HERE / "IDP_synthetic_sample_size_scaling_time_resource_v2_hereg.png"

# ---------------------------------------------------------------------------
# Series definitions
# ---------------------------------------------------------------------------

TRAIT_KEY_MAP: list[tuple[str, str]] = [
    ("HE_pipeline_PCA_capped_e2e",    "HE pipeline - PCA 128PCs (ours)"),
    ("HE_EVR08_e2e",                  "HE pipeline - PCA EVR0.8 (ours)"),
    ("QR_HE_e2e",                     "HE pipeline - QR (ours)"),
    ("REML_PCA_capped_e2e",           "REML pipeline"),
    ("tr_PinvG_e2e",                  "tr(P^-1G)"),
    ("SVD_PinvG_e2e",                 "SVD(P^-1G)"),
    ("blockwise_mvGREML_trG_trP_e2e", "blockwise_multivariate_GREML (trG/trP)"),
]

SAMPLE_KEY_MAP: list[tuple[str, str]] = [
    ("HE_pipeline_PCA_capped_e2e",    "HE pipeline - PCA 128PCs (ours)"),
    ("HE_EVR08_e2e",                  "HE pipeline - PCA EVR0.8 (ours)"),
    ("QR_HE_e2e",                     "HE pipeline - QR (ours)"),
    ("REML_PCA_capped_e2e",           "REML pipeline"),
    ("tr_PinvG_e2e",                  "tr(P^-1G)"),
    ("SVD_PinvG_e2e",                 "SVD(P^-1G)"),
    ("blockwise_mvGREML_trG_trP_e2e", "blockwise_multivariate_GREML (trG/trP)"),
]

COLOR_MAP: dict[str, str] = {
    "HE pipeline - PCA 128PCs (ours)":           "red",
    "HE pipeline - PCA EVR0.8 (ours)":           "deeppink",
    "HE pipeline - QR (ours)":                   "darkorange",
    "REML pipeline":                              "gold",
    "tr(P^-1G)":                                 "green",
    "SVD(P^-1G)":                                "blue",
    "blockwise_multivariate_GREML (trG/trP)":    "purple",
    "fastGWA + FUMA (estimated)":                "dimgray",
}

MARKERS = ["o", "X", "P", "s", "^", "D", "v"]

# Estimated fastGWA + FUMA times: linear scaling in p (each trait is independent GWAS).
# At p≈128 phenotypes: ~5 h measured; extrapolated linearly.
FASTGWA_TRAIT_P   = [100,   1_000,    10_000,     100_000]
FASTGWA_TRAIT_SEC = [18_000, 180_000, 1_800_000, 18_000_000]   # 5h, 50h, 500h, 5000h

# Sample scaling (p=1 000 fixed). fastGWA runtime is dominated by SNP count (m~8M),
# not sample count — approximately constant across the n range tested here.
FASTGWA_SAMPLE_N   = [719,     2_158,   6_474]
FASTGWA_SAMPLE_SEC = [180_000, 180_000, 180_000]                # ~50h flat


# ---------------------------------------------------------------------------
# Extrapolation helpers
# ---------------------------------------------------------------------------

def _loglog_quadratic(p: np.ndarray, t: np.ndarray, p_new: float) -> float:
    logp = np.log10(p.astype(float))
    logt = np.log10(np.maximum(t.astype(float), 1e-12))
    coef = np.polyfit(logp, logt, 2)
    return float(10 ** float(np.polyval(coef, np.log10(float(p_new)))))


def _append_100k_estimates(p_arr: list, series_time: dict, series_mb: dict,
                            dense_est: dict) -> tuple[list, dict, dict, dict]:
    tr_kernel = float(dense_est["estimate_100k"]["total_wall_sec_est"])
    mb_dense   = float(dense_est["estimate_100k"]["memory"]["approx_with_solve_temporaries_GB"]) * 1024.0
    ratio_svd_tr = 181719.06869811536 / 24270.171238115352
    ps = np.asarray(p_arr, dtype=float)
    p_new = 100_000.0
    new_time: dict = {}; new_mb: dict = {}; methods: dict = {}
    for _, label in TRAIT_KEY_MAP:
        t_m = np.asarray(series_time[label], dtype=float)
        mb_m = np.asarray(series_mb[label], dtype=float)
        if label in ("HE pipeline - PCA 128PCs (ours)", "HE pipeline - PCA EVR0.8 (ours)",
                     "HE pipeline - QR (ours)", "REML pipeline",
                     "blockwise_multivariate_GREML (trG/trP)"):
            est = _loglog_quadratic(ps, t_m, p_new)
            mb_est = _loglog_quadratic(ps, mb_m, p_new)
            if label == "blockwise_multivariate_GREML (trG/trP)":
                mb_est = min(mb_est, mb_dense * 1.05)
            new_time[label] = list(series_time[label]) + [est]
            new_mb[label]   = list(series_mb[label]) + [mb_est]
            methods[label] = "log10–log10 quadratic fit"
        elif label == "tr(P^-1G)":
            est = max(_loglog_quadratic(ps, t_m, p_new), tr_kernel)
            new_time[label] = list(series_time[label]) + [est]
            new_mb[label]   = list(series_mb[label]) + [mb_dense]
            methods[label] = f"max(quadratic, dense O(p³) kernel est={tr_kernel:.1f}s)"
        elif label == "SVD(P^-1G)":
            est = max(_loglog_quadratic(ps, t_m, p_new), tr_kernel * ratio_svd_tr)
            new_time[label] = list(series_time[label]) + [est]
            new_mb[label]   = list(series_mb[label]) + [mb_dense]
            methods[label] = "max(quadratic, tr_kernel×SVD/tr ratio)"
        else:
            new_time[label] = list(series_time[label])
            new_mb[label]   = list(series_mb[label])
    return p_arr + [100_000], new_time, new_mb, methods


# ---------------------------------------------------------------------------
# Trait scaling plot
# ---------------------------------------------------------------------------

def _positive(y: np.ndarray) -> np.ndarray:
    out = np.asarray(y, dtype=float).copy()
    finite = np.isfinite(out)
    out[finite] = np.maximum(out[finite], 1e-300)
    return out


def _plot_segments(ax, x: np.ndarray, y_raw: np.ndarray, extrap: np.ndarray,
                   color, marker: str, label: str) -> None:
    shown = False
    for i in range(len(x) - 1):
        if not (np.isfinite(y_raw[i]) and np.isfinite(y_raw[i+1])): continue
        ls = "--" if extrap[i+1] else "-"
        ax.loglog(x[i:i+2], _positive(y_raw[i:i+2]), linestyle=ls, color=color,
                  linewidth=2, label=label if not shown else None)
        shown = True
    for i in range(len(x)):
        if not np.isfinite(y_raw[i]): continue
        yp = float(_positive(np.array([y_raw[i]]))[0])
        ax.plot(x[i], yp, marker=marker, color=color, linestyle="",
                markersize=8 if extrap[i] else 7,
                markerfacecolor="white" if extrap[i] else color,
                markeredgewidth=1.4 if extrap[i] else 1.0, zorder=4)


def plot_trait_scaling(store: dict) -> None:
    ts = store["trait_scaling"]
    rows = sorted(ts["rows"], key=lambda r: int(r["p"]))
    p_values = [int(r["p"]) for r in rows]
    dense_est = store.get("dense_trace_extrapolation", {})

    series_time: dict = {}; series_mb: dict = {}; extrap: dict = {}
    for _, label in TRAIT_KEY_MAP:
        series_time[label] = []; series_mb[label] = []; extrap[label] = []
    for r in rows:
        for key, label in TRAIT_KEY_MAP:
            block = r.get(key, {})
            if isinstance(block, dict) and block.get("skipped"):
                series_time[label].append(float("nan"))
                series_mb[label].append(float("nan"))
                extrap[label].append(False)
            else:
                series_time[label].append(float(block.get("wall_sec", float("nan"))))
                series_mb[label].append(float(block.get("max_rss_mb", float("nan"))))
                extrap[label].append(bool(block.get("extrapolated", False)))

    # Fill skipped tr/SVD at 100k using physics-based estimates from dense_trace_extrapolation
    if dense_est and 100_000 in p_values:
        idx = p_values.index(100_000)
        p3 = np.asarray(p_values[:idx], dtype=float)
        tr_kernel = float(dense_est["estimate_100k"]["total_wall_sec_est"])
        mb_dense = float(dense_est["estimate_100k"]["memory"]["approx_with_solve_temporaries_GB"]) * 1024.0
        ratio_svd_tr = 181719.06869811536 / 24270.171238115352
        for label, is_svd in [("tr(P^-1G)", False), ("SVD(P^-1G)", True)]:
            if not np.isfinite(series_time[label][idx]):
                t_m = np.asarray(series_time[label][:idx], dtype=float)
                floor = tr_kernel * ratio_svd_tr if is_svd else tr_kernel
                est = max(_loglog_quadratic(p3, t_m, 100_000.0), floor)
                series_time[label][idx] = est
                series_mb[label][idx] = mb_dense
                extrap[label][idx] = True

    if p_values == [100, 1000, 10_000] and dense_est:
        p_values, series_time, series_mb, _ = _append_100k_estimates(
            p_values, series_time, series_mb, dense_est)
        for _, label in TRAIT_KEY_MAP:
            extrap[label].append(True)

    x = np.asarray(p_values, dtype=float)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.0, 5.6))
    rng = np.random.default_rng(42)
    for i, (_, label) in enumerate(TRAIT_KEY_MAP):
        color = COLOR_MAP[label]; marker = MARKERS[i % len(MARKERS)]
        _plot_segments(ax1, x, np.asarray(series_time[label]), np.asarray(extrap[label]), color, marker, label)
        mb_jitter = np.asarray(series_mb[label], dtype=float).copy()
        for k in range(len(mb_jitter)):
            if np.isfinite(mb_jitter[k]):
                mb_jitter[k] = mb_jitter[k] * (1.0 + 1e-4 * i) + rng.uniform(0, 0.02)
        _plot_segments(ax2, x, mb_jitter, np.asarray(extrap[label]), color, marker, label)

    # fastGWA + FUMA estimated reference line (time only; memory N/A)
    fgwa_x = np.asarray(FASTGWA_TRAIT_P, dtype=float)
    fgwa_t = np.asarray(FASTGWA_TRAIT_SEC, dtype=float)
    fgwa_color = COLOR_MAP["fastGWA + FUMA (estimated)"]
    ax1.loglog(fgwa_x, fgwa_t, linestyle="--", linewidth=2, color=fgwa_color,
               marker="*", markersize=9, label="fastGWA + FUMA (estimated)")
    # Add hour annotations at each point
    for px, pt in zip(FASTGWA_TRAIT_P, FASTGWA_TRAIT_SEC):
        h = pt / 3600
        label_str = f"{h:.0f}h" if h >= 1 else f"{pt:.0f}s"
        ax1.annotate(label_str, (float(px), float(pt)),
                     textcoords="offset points", xytext=(5, 4),
                     fontsize=7, color=fgwa_color)

    for ax, title, ylabel in [(ax1, "Time", "Wall time (s), log scale"),
                               (ax2, "Memory", "Peak RSS (MiB, process), log scale")]:
        ax.set_title(title); ax.set_xlabel("Trait count p (synthetic raw features)"); ax.set_ylabel(ylabel)
        ax.legend(loc="upper left", fontsize=8)
        ax.set_xticks(list(x)); ax.set_xticklabels([str(int(t)) for t in x])
    fig.tight_layout()
    fig.savefig(OUT_TRAIT, dpi=160); plt.close(fig)
    print(f"[plot] Saved → {OUT_TRAIT}")


# ---------------------------------------------------------------------------
# Sample scaling plot
# ---------------------------------------------------------------------------

def plot_sample_scaling(store: dict) -> None:
    ss = store["sample_scaling"]
    rows = sorted(ss["rows"], key=lambda r: r.get("n_target", 0))
    if not rows:
        print("[plot] sample_scaling: no rows — skipping"); return

    x = np.array([r.get("n_observed", r.get("n_target", 0)) for r in rows], dtype=float)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5.6))
    for i, (key, label) in enumerate(SAMPLE_KEY_MAP):
        if any(key not in r or (isinstance(r.get(key), dict) and r[key].get("failed")) for r in rows):
            continue
        yt = [float(r[key]["wall_sec"]) for r in rows]
        ym = [float(r[key]["max_rss_mb"]) for r in rows]
        c = COLOR_MAP[label]; m = MARKERS[i % len(MARKERS)]
        ax1.loglog(x, yt, marker=m, linewidth=2, color=c, label=label)
        ax2.loglog(x, ym, marker=m, linewidth=2, color=c, label=label)

    # fastGWA + FUMA estimated reference line (time only; memory N/A)
    fgwa_x = np.asarray(FASTGWA_SAMPLE_N, dtype=float)
    fgwa_t = np.asarray(FASTGWA_SAMPLE_SEC, dtype=float)
    fgwa_color = COLOR_MAP["fastGWA + FUMA (estimated)"]
    ax1.loglog(fgwa_x, fgwa_t, linestyle="--", linewidth=2, color=fgwa_color,
               marker="*", markersize=9, label="fastGWA + FUMA (estimated)")
    for nx, nt in zip(FASTGWA_SAMPLE_N, FASTGWA_SAMPLE_SEC):
        h = nt / 3600
        ax1.annotate(f"{h:.0f}h", (float(nx), float(nt)),
                     textcoords="offset points", xytext=(5, 4),
                     fontsize=7, color=fgwa_color)

    for ax, title, ylabel in [(ax1, "Time", "Wall time (s), log scale"),
                               (ax2, "Memory", "Peak RSS (MiB), log scale")]:
        ax.set_title(title)
        ax.set_xlabel("Sample count n (p=1000 fixed)"); ax.set_ylabel(ylabel)
        ax.set_xticks(list(x)); ax.set_xticklabels([str(int(v)) for v in x])
        ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_SAMPLE, dpi=160); plt.close(fig)
    print(f"[plot] Saved → {OUT_SAMPLE}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    store = json.loads(RESULTS_JSON.read_text())
    plot_trait_scaling(store)
    plot_sample_scaling(store)
