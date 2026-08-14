#!/usr/bin/env python3
"""
Generate scaling plots from results.json.

Writes:
  IDP_synthetic_trait_scaling_time_resource_v8_hereg.png
  IDP_synthetic_sample_size_scaling_time_resource_v2_hereg.png

Methods plotted: HE PCA 128PCs, REML, tr(P⁻¹G), SVD(P⁻¹G), blockwise mvGREML.
EVR0.8, QR, and whitened kernel are not included.
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

TRAIT_KEY_MAP: list[tuple[str, str]] = [
    ("HE_pipeline_PCA_capped_e2e",    "HE pipeline - PCA 128PCs (ours)"),
    ("REML_PCA_capped_e2e",           "REML pipeline"),
    ("tr_PinvG_e2e",                  "tr(P^-1G)"),
    ("SVD_PinvG_e2e",                 "SVD(P^-1G)"),
    ("blockwise_mvGREML_trG_trP_e2e", "blockwise_multivariate_GREML (trG/trP)"),
]

SAMPLE_KEY_MAP: list[tuple[str, str]] = list(TRAIT_KEY_MAP)

COLOR_MAP: dict[str, str] = {
    "HE pipeline - PCA 128PCs (ours)":           "red",
    "REML pipeline":                              "gold",
    "tr(P^-1G)":                                 "green",
    "SVD(P^-1G)":                                "blue",
    "blockwise_multivariate_GREML (trG/trP)":    "purple",
    "fastGWA + FUMA (estimated)":                "dimgray",
}

MARKER_MAP: dict[str, str] = {
    "HE pipeline - PCA 128PCs (ours)":           "o",
    "REML pipeline":                              "s",
    "tr(P^-1G)":                                 "^",
    "SVD(P^-1G)":                                "D",
    "blockwise_multivariate_GREML (trG/trP)":    "v",
}

FASTGWA_TRAIT_P   = [100,   1_000,    10_000,     100_000]
FASTGWA_TRAIT_SEC = [18_000, 180_000, 1_800_000, 18_000_000]
# Match observed n in results.json sample_scaling rows (approx. 719 / 2158 / 6474)
FASTGWA_SAMPLE_N   = [717,     2_153,   6_458]
FASTGWA_SAMPLE_SEC = [180_000, 180_000, 180_000]


def _loglog_quadratic(p: np.ndarray, t: np.ndarray, p_new: float) -> float:
    logp = np.log10(p.astype(float))
    logt = np.log10(np.maximum(t.astype(float), 1e-12))
    coef = np.polyfit(logp, logt, 2)
    return float(10 ** float(np.polyval(coef, np.log10(float(p_new)))))


def _estimate_100k_for_label(label: str, p_arr: list, t_m: np.ndarray, mb_m: np.ndarray,
                             tr_kernel: float, mb_dense: float, ratio_svd_tr: float
                             ) -> tuple[float, float, str] | None:
    """Return (time_est, mb_est, method_note) or None if this series has no estimate."""
    ps = np.asarray(p_arr, dtype=float)
    p_new = 100_000.0
    finite_t = np.isfinite(t_m)
    finite_m = np.isfinite(mb_m)
    if label in ("HE pipeline - PCA 128PCs (ours)", "REML pipeline",
                 "blockwise_multivariate_GREML (trG/trP)"):
        if finite_t.sum() < 2:
            return None
        est = _loglog_quadratic(ps[finite_t], t_m[finite_t], p_new)
        mb_est = _loglog_quadratic(ps[finite_m], mb_m[finite_m], p_new) if finite_m.sum() >= 2 else mb_dense
        if label == "blockwise_multivariate_GREML (trG/trP)":
            mb_est = min(mb_est, mb_dense * 1.05)
        return est, mb_est, "log10–log10 quadratic fit"
    if label == "tr(P^-1G)":
        if finite_t.sum() < 2:
            return None
        est = max(_loglog_quadratic(ps[finite_t], t_m[finite_t], p_new), tr_kernel)
        return est, mb_dense, f"max(quadratic, dense O(p³) kernel est={tr_kernel:.1f}s)"
    if label == "SVD(P^-1G)":
        if finite_t.sum() < 2:
            return None
        est = max(_loglog_quadratic(ps[finite_t], t_m[finite_t], p_new), tr_kernel * ratio_svd_tr)
        return est, mb_dense, "max(quadratic, tr_kernel×SVD/tr ratio)"
    return None


def _append_100k_estimates(p_arr: list, series_time: dict, series_mb: dict,
                            extrap: dict, dense_est: dict
                            ) -> tuple[list, dict, dict, dict, dict]:
    """Fill or append p=100k. Never leave a NaN gap before the estimate (that breaks the line)."""
    tr_kernel = float(dense_est["estimate_100k"]["total_wall_sec_est"])
    mb_dense = float(dense_est["estimate_100k"]["memory"]["approx_with_solve_temporaries_GB"]) * 1024.0
    ratio_svd_tr = 181719.06869811536 / 24270.171238115352

    p_arr = list(p_arr)
    new_time = {k: list(v) for k, v in series_time.items()}
    new_mb = {k: list(v) for k, v in series_mb.items()}
    new_ex = {k: list(v) for k, v in extrap.items()}
    methods: dict = {}

    has_100k = 100_000 in p_arr
    idx_100k = p_arr.index(100_000) if has_100k else None
    # Fit using measured points only (exclude existing 100k slot if present)
    fit_ps = [p for i, p in enumerate(p_arr) if not (has_100k and i == idx_100k)]

    for _, label in TRAIT_KEY_MAP:
        t_m = np.asarray(new_time[label], dtype=float)
        mb_m = np.asarray(new_mb[label], dtype=float)
        if has_100k:
            # Keep measured 100k; only replace NaN/skipped with an estimate.
            if np.isfinite(t_m[idx_100k]):
                new_ex[label][idx_100k] = False
                methods[label] = "measured"
                continue
            fit_t = np.asarray([t_m[i] for i, p in enumerate(p_arr) if i != idx_100k], dtype=float)
            fit_m = np.asarray([mb_m[i] for i, p in enumerate(p_arr) if i != idx_100k], dtype=float)
            est = _estimate_100k_for_label(label, fit_ps, fit_t, fit_m, tr_kernel, mb_dense, ratio_svd_tr)
            if est is None:
                continue
            t_est, mb_est, note = est
            new_time[label][idx_100k] = t_est
            new_mb[label][idx_100k] = mb_est
            new_ex[label][idx_100k] = True
            methods[label] = note
        else:
            est = _estimate_100k_for_label(label, p_arr, t_m, mb_m, tr_kernel, mb_dense, ratio_svd_tr)
            if est is None:
                continue
            t_est, mb_est, note = est
            new_time[label].append(t_est)
            new_mb[label].append(mb_est)
            new_ex[label].append(True)
            methods[label] = note

    if not has_100k:
        # All series that received an estimate were appended; pad any that did not
        for _, label in TRAIT_KEY_MAP:
            if len(new_time[label]) == len(p_arr):
                new_time[label].append(float("nan"))
                new_mb[label].append(float("nan"))
                new_ex[label].append(False)
        p_arr = p_arr + [100_000]

    return p_arr, new_time, new_mb, methods, new_ex


def _positive(y: np.ndarray) -> np.ndarray:
    out = np.asarray(y, dtype=float).copy()
    finite = np.isfinite(out)
    out[finite] = np.maximum(out[finite], 1e-300)
    return out


def _plot_segments(ax, x: np.ndarray, y_raw: np.ndarray, extrap: np.ndarray,
                   color, marker: str, label: str) -> None:
    shown = False
    for i in range(len(x) - 1):
        if not (np.isfinite(y_raw[i]) and np.isfinite(y_raw[i+1])):
            continue
        ls = "--" if extrap[i+1] else "-"
        ax.loglog(x[i:i+2], _positive(y_raw[i:i+2]), linestyle=ls, color=color,
                  linewidth=2, label=label if not shown else None)
        shown = True
    for i in range(len(x)):
        if not np.isfinite(y_raw[i]):
            continue
        yp = float(_positive(np.array([y_raw[i]]))[0])
        ax.plot(x[i], yp, marker=marker, color=color, linestyle="",
                markersize=8 if extrap[i] else 7,
                markerfacecolor="white" if extrap[i] else color,
                markeredgewidth=1.4 if extrap[i] else 1.0, zorder=4)


def plot_trait_scaling(store: dict, *, out_path: Path | None = None) -> None:
    out_path = out_path or OUT_TRAIT
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
            if not isinstance(block, dict) or block.get("skipped") or block.get("failed") or not block:
                series_time[label].append(float("nan"))
                series_mb[label].append(float("nan"))
                extrap[label].append(False)
            else:
                series_time[label].append(float(block.get("wall_sec", float("nan"))))
                series_mb[label].append(float(block.get("max_rss_mb", float("nan"))))
                extrap[label].append(bool(block.get("extrapolated", False)))

    if dense_est:
        p_plot, series_time, series_mb, _, extrap = _append_100k_estimates(
            p_values, series_time, series_mb, extrap, dense_est)
    else:
        p_plot = p_values

    x = np.asarray(p_plot, dtype=float)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    for _, label in TRAIT_KEY_MAP:
        yt = np.asarray(series_time[label], dtype=float)
        ym = np.asarray(series_mb[label], dtype=float)
        if len(yt) != len(x) or not np.any(np.isfinite(yt)):
            continue
        ex = np.asarray(extrap[label], dtype=bool)
        if len(ex) != len(x):
            ex = np.resize(ex, len(x))
        c = COLOR_MAP[label]; m = MARKER_MAP[label]
        _plot_segments(ax1, x, yt, ex, c, m, label)
        _plot_segments(ax2, x, ym, ex, c, m, label)

    fgwa_x = np.asarray(FASTGWA_TRAIT_P, dtype=float)
    fgwa_t = np.asarray(FASTGWA_TRAIT_SEC, dtype=float)
    fgwa_color = COLOR_MAP["fastGWA + FUMA (estimated)"]
    ax1.loglog(fgwa_x, fgwa_t, linestyle="--", linewidth=2, color=fgwa_color,
               marker="*", markersize=9, label="fastGWA + FUMA (estimated)")
    for px, pt in zip(FASTGWA_TRAIT_P, FASTGWA_TRAIT_SEC):
        h = pt / 3600
        ax1.annotate(f"{h:.0f}h" if h < 100 else f"{h:.0f}h",
                     (float(px), float(pt)), textcoords="offset points",
                     xytext=(5, 4), fontsize=7, color=fgwa_color)

    for ax, title, ylabel in [(ax1, "Time", "Wall time (s), log scale"),
                               (ax2, "Memory", "Peak RSS (MiB), log scale")]:
        ax.set_title(title)
        ax.set_xlabel("Trait count p (n≈2158 fixed)"); ax.set_ylabel(ylabel)
        ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160); plt.close(fig)
    print(f"[plot] Saved → {out_path}")


def plot_sample_scaling(store: dict, *, out_path: Path | None = None) -> None:
    out_path = out_path or OUT_SAMPLE
    rows = sorted(store["sample_scaling"]["rows"], key=lambda r: int(r.get("n_target", 0)))
    x = np.asarray([int(r.get("n_observed", r.get("n_target", 0))) for r in rows], dtype=float)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    for key, label in SAMPLE_KEY_MAP:
        yt, ym = [], []
        for r in rows:
            block = r.get(key, {})
            if (not isinstance(block, dict) or block.get("failed") or block.get("skipped")
                    or block.get("wall_sec") is None):
                yt.append(float("nan"))
                ym.append(float("nan"))
            else:
                yt.append(float(block["wall_sec"]))
                ym.append(float(block.get("max_rss_mb", float("nan"))))
        yt_a = np.asarray(yt, dtype=float)
        ym_a = np.asarray(ym, dtype=float)
        if not np.any(np.isfinite(yt_a)):
            continue
        ex = np.zeros(len(x), dtype=bool)
        c = COLOR_MAP[label]; m = MARKER_MAP[label]
        _plot_segments(ax1, x, yt_a, ex, c, m, label)
        _plot_segments(ax2, x, ym_a, ex, c, m, label)

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
    fig.savefig(out_path, dpi=160); plt.close(fig)
    print(f"[plot] Saved → {out_path}")


if __name__ == "__main__":
    store = json.loads(RESULTS_JSON.read_text())
    plot_trait_scaling(store)
    plot_sample_scaling(store)
