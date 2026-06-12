#!/usr/bin/env python3
"""
Generate 5-panel scatter plots of heritability h² vs GWAS loci count.

Reads all data from results.json. Writes:
  scatter_5panel_hereg_fuma.png   — x-axis: minP loci count (22 models, FUMA-clumped)
  scatter_5panel_hereg_jagwas.png — x-axis: JAGWAS lmm_only loci count (21 models, FUMA-clumped)

5 panels per figure:
  1. HE pipeline – PCA 128PCs    (sum h²)
  2. HE pipeline – PCA EVR0.8   (mean h² up to EVR≥0.8 cutoff)
  3. HE pipeline – QR            (mean h² = sum_h2 / n_dims)
  4. REML pipeline – PCA 128PCs  (reml_sum_h2)
  5. blockwise mvGREML           (tr(G)/tr(P))
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE           = Path(__file__).resolve().parent
RESULTS_JSON   = HERE / "results.json"
HEREG_RUNS_DIR = HERE / "hereg_runs"    # local cache written by run_pipeline.py
OUT_FUMA       = HERE / "scatter_5panel_hereg_fuma.png"
OUT_JAGWAS     = HERE / "scatter_5panel_hereg_jagwas.png"

EVR_TARGET = 0.8
FIT_LINE_COLOR = "#7eb8e6"
DISPLAY_RANDOM_X_SHIFT = 1.35
FORCE_BOTTOM_RIGHT = frozenset({"ViT_T1_Random", "ViT_T1_PCA"})

PALETTE = {
    "ViT": "#ff7f0e", "CNN": "#2ca02c", "MoCo": "#9467bd",
    "Random": "#d62728", "Surf": "#7f7f7f", "Other": "#7f7f7f",
    "IDP": "#000000", "IDP_T1": "#000000", "IDP_T2": "#000000",
}
LEGEND_FAMILIES = [
    ("ViT", "Vol_ViT"), ("CNN", "Vol_CNN"), ("MoCo", "Vol_MoCo"),
    ("Random", "Vol_Random"), ("Surf", "Surf"), ("IDP", "IDP"),
]


# ---------------------------------------------------------------------------
# Label helpers
# ---------------------------------------------------------------------------

def family(raw: str) -> str:
    l = raw.lower()
    if l == "idp_t1": return "IDP_T1"
    if l == "idp_t2": return "IDP_T2"
    if "t1_idp_pheno_discovery" in l: return "IDP"
    if "z_graph" in l: return "Surf"
    if "mocov2" in l:  return "MoCo"
    if "cnn" in l:     return "CNN"
    if "random" in l:  return "Random"
    if "vit" in l:     return "ViT"
    return "Other"


def clean_label(raw: str) -> str:
    l = raw.lower().replace("\\", "/")
    if l == "idp_t1" or l.endswith("/idp/t1"): return "IDP_T1"
    if l == "idp_t2" or l.endswith("/idp/t2"): return "IDP_T2"
    if "t1_idp_pheno_discovery" in l: return "IDP"
    if "z_graph_fsaverage4" in l or l.startswith("z_graph_fsaverage"):
        return "surf_mesh_mask75" if ("mask0p75" in l or "mask75" in l) else "surf_mesh_mask25"
    if "graphunet" in l: return "z_GraphUNet"
    if "residual_distill" in l or "fusion_residual" in l: return "z_FusionDistill"
    if "vit_t1" in l and "384" in l and "vit_t2" not in l and "mae" not in l: return "ViT_T1_PCA"
    if "vit_t2" in l and "384" in l: return "ViT_T2_PCA"
    if "mae_nomask" in l: return "ViT_T1_MAE"
    if l.rstrip("/").endswith("cnn_adni"): return "CNN_T1_ADNI"
    if l.rstrip("/").endswith("cnn_replicate1"): return "CNN_T1_replicate1"
    if "mocov2" in l:
        return "MoCo_T2" if ("t2_5std" in l or ("t2" in l and "5std" in l)) else "MoCo_T1"
    if "vit_t1_fixed_adni" in l: return "ViT_T1_ADNI"
    if "vit_t1_fixed_random" in l: return "ViT_T1_Random"
    m = re.search(r"cv_fold(\d+)", l)
    if m: return f"ViT_T1_cv_fold{m.group(1)}"
    if "vit_t1_fixed_replicate2" in l: return "ViT_T1_replicate2"
    if "vit_t1_fixed_replicate" in l: return "ViT_T1_replicate"
    s = raw
    s = s.replace("contrast_learning_mocov2", "MoCo").replace("partial_fit", "T1")
    s = s.replace("fixed_", "").replace("fixed", "")
    s = s.replace("vit_", "ViT_").replace("vit", "ViT")
    s = s.replace("cnn_", "CNN_").replace("cnn", "CNN")
    s = s.replace("random", "Random")
    s = re.sub(r"(?i)\bt1\b", "T1", s)
    s = re.sub(r"(?i)\bt2\b", "T2", s)
    if "T1" not in s and "T2" not in s: s = f"{s}_T1"
    return s.replace("__", "_").strip("_")


# ---------------------------------------------------------------------------
# Load data from results.json
# ---------------------------------------------------------------------------

_store = json.loads(RESULTS_JSON.read_text())

# Use local hereg_runs cache if available (from run_pipeline.py), else fall back
_summary_path = HEREG_RUNS_DIR / "hereg_summary.json"
results: dict[str, dict] = (
    json.loads(_summary_path.read_text()) if _summary_path.is_file()
    else _store["h2"]
)

# loci counts from results.json["loci"]
MINP_LOCI: dict[str, float] = {
    k: float(v["minp"])
    for k, v in _store["loci"].items()
    if v.get("minp") is not None
}
JAGWAS_LOCI: dict[str, float] = {
    k: float(v["jagwas"])
    for k, v in _store["loci"].items()
    if v.get("jagwas") is not None
}

# Blockwise tr(G)/tr(P) — now stored in h2[label]["blockwise"]
BLOCKWISE_BY_LABEL: dict[str, dict] = {}
BLOCKWISE_JAGWAS_BY_LABEL: dict[str, dict] = {}
for _raw, _entry in _store["h2"].items():
    _bw = _entry.get("blockwise", {})
    _trv_raw = _bw.get("trG_over_trP")
    if _trv_raw is None:
        continue
    try:
        _trv = float(_trv_raw)
        if not np.isfinite(_trv):
            continue
    except (TypeError, ValueError):
        continue
    _pt = {"y": _trv, "label": clean_label(_raw.lower()), "fam": family(_raw.lower())}
    if (_lf := MINP_LOCI.get(_raw)) is not None:
        BLOCKWISE_BY_LABEL[_raw] = {**_pt, "x": _lf}
    if (_lj := JAGWAS_LOCI.get(_raw)) is not None:
        BLOCKWISE_JAGWAS_BY_LABEL[_raw] = {**_pt, "x": _lj}


# ---------------------------------------------------------------------------
# Metric extractors
# ---------------------------------------------------------------------------

def _finite(v):
    try:
        f = float(v)
        return f if np.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def _evr80_mean(label: str) -> float | None:
    run_dir = HEREG_RUNS_DIR / re.sub(r"[/\\]", "_", label) / "pca"
    csv_path = run_dir / "hereg_h2_summary.csv"
    evr_path = run_dir / "pca_features" / "pca_explained_variance_ratio.json"
    if not csv_path.is_file() or not evr_path.is_file():
        return None
    df = pd.read_csv(csv_path)
    evr = json.loads(evr_path.read_text())
    if not evr or df.empty:
        return None
    evr_arr = np.array(evr[:len(df)], dtype=np.float64)
    cumevr = np.cumsum(evr_arr)
    k = int(np.searchsorted(cumevr, EVR_TARGET)) + 1
    k = min(k, len(df))
    h2_vals = df["h2"].iloc[:k].dropna().values.astype(np.float64)
    valid = h2_vals[np.isfinite(h2_vals)]
    return float(valid.sum() / k) if len(valid) > 0 else None


def _hereg_pca_sum(label, d):  return d.get("pca", {}).get("sum_h2")
def _hereg_evr80(label, d):
    # Prefer pre-computed value stored in results.json; fall back to local files
    v = d.get("pca", {}).get("mean_h2_evr08")
    if v is not None: return v
    return _evr80_mean(label)
def _hereg_qr_mean(label, d):
    m = d.get("qr", {})
    s = _finite(m.get("sum_h2"))
    n = m.get("n_dims") or m.get("n_phenotype_columns") or 128
    return s / n if s is not None else None
def _reml_pca_sum(label, d):   return d.get("pca", {}).get("reml_sum_h2")


PANELS = [
    ("HE pipeline - PCA 128PCs", "Sum h²",      _hereg_pca_sum),
    ("HE pipeline - PCA EVR0.8", "Mean h²",     _hereg_evr80),
    ("HE pipeline - QR",         "Mean h²",     _hereg_qr_mean),
    ("REML pipeline",            "REML sum h²", _reml_pca_sum),
]


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def build_pts(x_source: dict, metric_fn) -> list[dict]:
    pts = []
    for label, d in results.items():
        xv = x_source.get(label)
        if xv is None:
            continue
        yv = _finite(metric_fn(label, d))
        if yv is None:
            continue
        pts.append({"x": float(xv), "y": yv,
                    "label": clean_label(label), "fam": family(label)})
    return pts


def draw_panel(ax, pts, *, y_label, title_base, x_label):
    if not pts:
        ax.set_title(f"{title_base}\n(no data)")
        return
    x = np.array([p["x"] for p in pts])
    y = np.array([p["y"] for p in pts])

    if len(x) >= 2 and np.std(x) > 0 and np.std(y) > 0:
        m, b = np.polyfit(x, y, 1)
        r_pearson = float(np.corrcoef(x, y)[0, 1])
    else:
        m = b = float("nan"); r_pearson = None

    x_span = float(np.ptp(x)) or 1.0
    y_span = float(np.ptp(y)) or 1.0
    x_pad = max(2.0, 0.07 * x_span)
    y_pad = max(0.02 * (abs(y).max() or 1.0), 0.08 * y_span, 1e-6)
    ax.set_xlim(x.min() - x_pad, x.max() + x_pad + DISPLAY_RANDOM_X_SHIFT)
    ax.set_ylim(y.min() - y_pad, y.max() + y_pad)

    for fam, legend_name in LEGEND_FAMILIES:
        if fam == "Surf":
            sub = [p for p in pts if p["fam"] in ("Surf", "Other")]
        elif fam == "IDP":
            sub = [p for p in pts if p["fam"] in ("IDP", "IDP_T1", "IDP_T2")]
        else:
            sub = [p for p in pts if p["fam"] == fam]
        if not sub: continue
        xs_plot = [p["x"] + (DISPLAY_RANDOM_X_SHIFT if fam == "Random" else 0.0) for p in sub]
        ec = "none" if fam == "IDP" else "white"
        ax.scatter(xs_plot, [p["y"] for p in sub], s=42, alpha=0.92,
                   color=PALETTE.get(fam, "#7f7f7f"), edgecolors=ec, linewidths=0.5,
                   label=legend_name, zorder=3)

    xlim = ax.get_xlim()
    xs_line = np.linspace(xlim[0], xlim[1], 200)
    if np.isfinite(m) and np.isfinite(b):
        ax.plot(xs_line, m * xs_line + b, linewidth=1.8, color=FIT_LINE_COLOR,
                label=f"y={m:.4g}x+{b:.4g}", zorder=2)

    x_median = float(np.median(x))
    for p in pts:
        x_ann = p["x"] + (DISPLAY_RANDOM_X_SHIFT if p["fam"] == "Random" else 0.0)
        put_right = (p["x"] <= x_median) or (p["label"] in FORCE_BOTTOM_RIGHT)
        ax.annotate(p["label"], (x_ann, p["y"]),
                    textcoords="offset points",
                    xytext=(5, 0) if put_right else (-5, 0),
                    ha="left" if put_right else "right",
                    va="center", fontsize=6.5, color="#1f1f1f",
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.78),
                    zorder=4, clip_on=True)

    ax.set_xlabel(x_label, fontsize=8)
    ax.set_ylabel(y_label, fontsize=12)
    ax.tick_params(labelsize=7)
    r_str = f"Pearson r = {r_pearson:.3f}" if r_pearson is not None else "Pearson r = NA"
    ax.set_title(f"{title_base}\n{r_str}", fontsize=13)
    ax.grid(alpha=0.2, zorder=1)
    handles, labels_leg = ax.get_legend_handles_labels()
    ax.legend(handles, labels_leg, fontsize=6.5, loc="best", framealpha=0.88)


def make_figure(x_source: dict, x_label: str, out_path: Path,
                bw_pts: list[dict]) -> None:
    panel_w, panel_h = 8.8 * (2/3), 5.9 * (2/3)
    fig, axes = plt.subplots(3, 2, figsize=(panel_w * 2, panel_h * 3))
    for i, (title_base, y_label, metric_fn) in enumerate(PANELS):
        row, col = divmod(i, 2)
        pts = build_pts(x_source, metric_fn)
        draw_panel(axes[row, col], pts, y_label=y_label,
                   title_base=title_base, x_label=x_label)
    draw_panel(axes[2, 0], bw_pts, y_label="tr(G)/tr(P)",
               title_base="blockwise_multivariate_GREML (trG/trP)", x_label=x_label)
    axes[2, 1].set_visible(False)
    fig.tight_layout(pad=1.2)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] Saved → {out_path}")


make_figure(MINP_LOCI,   "minP loci count (22 models)",
            OUT_FUMA,   list(BLOCKWISE_BY_LABEL.values()))
make_figure(JAGWAS_LOCI, "JAGWAS lmm_only loci count (21 models)",
            OUT_JAGWAS, list(BLOCKWISE_JAGWAS_BY_LABEL.values()))
