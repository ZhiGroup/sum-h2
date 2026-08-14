#!/usr/bin/env python3
"""
Generate scaling plots from results.json.

Writes:
  IDP_synthetic_trait_scaling_time_resource_v8_hereg.png
  IDP_synthetic_sample_size_scaling_time_resource_v2_hereg.png

Methods plotted: HE PCA 128PCs, REML, tr(P⁻¹G), SVD(P⁻¹G), blockwise mvGREML
(plus whitened kernel when present). EVR0.8 and QR are not included.
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
    ("HE_whitened_kernel_e2e",        "HE pipeline - whitened kernel (ours)"),
]

SAMPLE_KEY_MAP: list[tuple[str, str]] = list(TRAIT_KEY_MAP)

COLOR_MAP: dict[str, str] = {
    "HE pipeline - PCA 128PCs (ours)":           "red",
    "REML pipeline":                              "gold",
    "tr(P^-1G)":                                 "green",
    "SVD(P^-1G)":                                "blue",
    "blockwise_multivariate_GREML (trG/trP)":    "purple",
    "HE pipeline - whitened kernel (ours)":      "black",
    "fastGWA + FUMA (estimated)":                "dimgray",
}

MARKER_MAP: dict[str, str] = {
    "HE pipeline - PCA 128PCs (ours)":           "o",
    "REML pipeline":                              "s",
    "tr(P^-1G)":                                 "^",
    "SVD(P^-1G)":                                "D",
    "blockwise_multivariate_GREML (trG/trP)":    "v",
    "HE pipeline - whitened kernel (ours)":      "o",
}

FASTGWA_TRAIT_P   = [100,   1_000,    10_000,     100_000]
FASTGWA_TRAIT_SEC = [18_000, 180_000, 1_800_000, 18_000_000]
FASTGWA_SAMPLE_N   = [719,     2_158,   6_474]
FASTGWA_SAMPLE_SEC = [180_000, 180_000, 180_000]


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
        if label in ("HE pipeline - PCA 128PCs (ours)", "REML pipeline",
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
        elif label == "HE pipeline - whitened kernel (ours)":
            # may be all-nan if stripped from some rows
            if not np.any(np.isfinite(t_m)):
                new_time[label] = list(series_time[label])
                new_mb[label] = list(series_mb[label])
                continue
            est = max(_loglog_quadratic(ps[np.isfinite(t_m)], t_m[np.isfinite(t_m)], p_new), tr_kernel)
            new_time[label] = list(series_time[label]) + [est]
            new_mb[label]   = list(series_mb[label]) + [mb_dense]
            methods[label] = f"max(quadratic, dense O(p³) kernel est={tr_kernel:.1f}s)"
        else:
            new_time[label] = list(series_time[label])
            new_mb[label]   = list(series_mb[label])
    return p_arr + [100_000], new_time, new_mb, methods


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
        p_plot, series_time, series_mb, _ = _append_100k_estimates(
            p_values, series_time, series_mb, dense_est)
        for label in series_time:
            extrap[label] = list(extrap.get(label, [])) + [True]
    else:
        p_plot = p_values

    x = np.asarray(p_plot, dtype=float)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    for _, label in TRAIT_KEY_MAP:
        yt = np.asarray(series_time[label], dtype=float)
        ym = np.asarray(series_mb[label], dtype=float)
        if not np.any(np.isfinite(yt)):
            continue
        ex = np.asarray(extrap[label], dtype=bool)
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
        if any(key not in r or (isinstance(r.get(key), dict) and (r[key].get("failed") or r[key].get("skipped")))
               for r in rows):
            # still plot if some rows have data
            if not any(isinstance(r.get(key), dict) and r[key].get("wall_sec") is not None for r in rows):
                continue
        yt, ym = [], []
        ok = True
        for r in rows:
            block = r.get(key, {})
            if not isinstance(block, dict) or block.get("failed") or block.get("skipped") or "wall_sec" not in block:
                ok = False
                break
            yt.append(float(block["wall_sec"]))
            ym.append(float(block["max_rss_mb"]))
        if not ok:
            continue
        c = COLOR_MAP[label]; m = MARKER_MAP[label]
        ax1.loglog(x, yt, marker=m, linewidth=2, color=c, label=label)
        ax2.loglog(x, ym, marker=m, linewidth=2, color=c, label=label)

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
