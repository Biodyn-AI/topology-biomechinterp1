from __future__ import annotations

import csv
import json
import math
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from umap import UMAP


PROJECT_ROOT = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/"
    "biodyn-work/subproject_41_claude_topology_hypothesis_screening_autoloop"
)
ITER_0009 = PROJECT_ROOT / "iterations/iter_0009"
ITER_0010 = PROJECT_ROOT / "iterations/iter_0010"
FIG_DIR = PROJECT_ROOT / "reports/figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

EMB_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/"
    "subproject_38_geometric_residual_stream_interpretability/"
    "implementation/outputs/cycle1_main/layer_gene_embeddings.npy"
)
EDGE_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/"
    "subproject_38_geometric_residual_stream_interpretability/"
    "implementation/outputs/cycle1_main/cycle1_edge_dataset.tsv"
)
GENE2GO_PKL = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/"
    "single_cell_mechinterp/data/perturb/gene2go_all.pkl"
)
TRRUST_PATH = Path(
    "/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work/"
    "single_cell_mechinterp/external/networks/trrust_human.tsv"
)


def _safe_float(value: str) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except Exception:
        return None


def _plot_spectral_trajectory() -> Path:
    spectral_path = ITER_0010 / "h03_spectral_profile.csv"
    rows = list(csv.DictReader(spectral_path.open()))

    layers = [int(r["layer"]) for r in rows]
    sv1 = [float(r["SV1_var_frac"]) for r in rows]
    sv2 = [float(r["SV2_var_frac"]) for r in rows]
    sv3 = [float(r["SV3_var_frac"]) for r in rows]
    ratio = [float(r["SV1_SV2_ratio"]) for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    axes[0].plot(layers, sv1, marker="o", label="SV1")
    axes[0].plot(layers, sv2, marker="o", label="SV2")
    axes[0].plot(layers, sv3, marker="o", label="SV3")
    axes[0].set_title("Variance Fractions Across Layers")
    axes[0].set_xlabel("Layer")
    axes[0].set_ylabel("Fraction of variance")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    axes[1].plot(layers, ratio, marker="o", color="tab:purple")
    axes[1].set_title("SV1/SV2 Magnitude Ratio")
    axes[1].set_xlabel("Layer")
    axes[1].set_ylabel("SV1/SV2")
    axes[1].grid(alpha=0.3)

    fig.suptitle("Claude Iteration 0010: Spectral Geometry")
    fig.tight_layout()

    out = FIG_DIR / "claude_spectral_trajectory_iter0010.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def _heatmap(
    matrix: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
    title: str,
    cbar_label: str,
    out_path: Path,
    cmap: str,
    center_zero: bool = False,
) -> Path:
    fig, ax = plt.subplots(figsize=(12, 5.2))

    if center_zero:
        vmax = np.nanmax(np.abs(matrix))
        vmin = -vmax
    else:
        vmin = np.nanmin(matrix)
        vmax = np.nanmax(matrix)

    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xlabel("Layer")
    ax.set_ylabel("Compartment")
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)

    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def _build_iter0009_sv1_matrices() -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
    csv_path = ITER_0009 / "h02_layer_compartment_map.csv"
    rows = list(csv.DictReader(csv_path.open()))

    compartments = sorted({r["compartment"] for r in rows})
    layers = sorted({int(r["layer"]) for r in rows})

    signed_log2_or = np.full((len(compartments), len(layers)), np.nan, dtype=float)
    neglog10_emp = np.full((len(compartments), len(layers)), np.nan, dtype=float)

    comp_idx = {c: i for i, c in enumerate(compartments)}
    layer_idx = {l: i for i, l in enumerate(layers)}

    for r in rows:
        c = r["compartment"]
        l = int(r["layer"])
        or_val = _safe_float(r.get("OR", ""))
        emp = _safe_float(r.get("emp_p", ""))
        pole = (r.get("pole") or "").strip().lower()

        if or_val is not None and or_val > 0:
            sign = 1.0 if pole == "top" else -1.0
            signed_log2_or[comp_idx[c], layer_idx[l]] = sign * math.log2(or_val)
        if emp is not None and emp > 0:
            neglog10_emp[comp_idx[c], layer_idx[l]] = -math.log10(emp)
        elif emp == 0:
            # Lower-bound of null resolution (iter_0009 used N=200 for this scan)
            neglog10_emp[comp_idx[c], layer_idx[l]] = -math.log10(1.0 / 200.0)

    return signed_log2_or, neglog10_emp, compartments, [str(l) for l in layers]


def _build_iter0010_sv2_matrices() -> tuple[np.ndarray, np.ndarray, list[str], list[str]]:
    json_path = ITER_0010 / "h01_sv2sv3_layer_compartment_scan.json"
    payload = json.loads(json_path.read_text())

    compartments = [
        "mitochondrion",
        "ER_lumen",
        "extracellular_vesicle",
        "cytoskeleton",
        "plasma_membrane",
        "nucleus",
        "secreted",
        "cytoplasm",
        "lysosome",
    ]
    layers = list(range(12))

    signed_log2_or = np.full((len(compartments), len(layers)), np.nan, dtype=float)
    neglog10_emp = np.full((len(compartments), len(layers)), np.nan, dtype=float)

    for i, comp in enumerate(compartments):
        for j, layer in enumerate(layers):
            cell = payload["SV2"][str(layer)][comp]
            if cell.get("skip"):
                continue
            or_val = float(cell["or"])
            emp = float(cell["emp_p"])
            pole = str(cell.get("pole", "")).lower()

            if or_val > 0:
                sign = 1.0 if pole == "top" else -1.0
                signed_log2_or[i, j] = sign * math.log2(or_val)

            if emp > 0:
                neglog10_emp[i, j] = -math.log10(emp)
            else:
                # N=200 shuffles here, so minimum resolvable p is 1/200
                neglog10_emp[i, j] = -math.log10(1.0 / 200.0)

    return signed_log2_or, neglog10_emp, compartments, [str(l) for l in layers]


def _plot_trrust_null_vs_obs() -> Path:
    result_path = ITER_0010 / "h02_trrust_copole_result.json"
    null_sv1_path = ITER_0010 / "h02_trrust_copole_null_sv1.npy"
    null_sv2_path = ITER_0010 / "h02_trrust_copole_null_sv2.npy"

    result = json.loads(result_path.read_text())
    null_sv1 = np.load(null_sv1_path)
    null_sv2 = np.load(null_sv2_path)

    obs1 = float(result["SV1"]["obs_copole_rate"])
    obs2 = float(result["SV2"]["obs_copole_rate"])

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)

    axes[0].hist(null_sv1, bins=25, alpha=0.8, color="tab:blue")
    axes[0].axvline(obs1, color="red", linestyle="--", linewidth=2)
    axes[0].set_title("TRRUST Co-pole Null (SV1)")
    axes[0].set_xlabel("Co-pole rate")
    axes[0].set_ylabel("Count")

    axes[1].hist(null_sv2, bins=25, alpha=0.8, color="tab:green")
    axes[1].axvline(obs2, color="red", linestyle="--", linewidth=2)
    axes[1].set_title("TRRUST Co-pole Null (SV2)")
    axes[1].set_xlabel("Co-pole rate")

    fig.suptitle(
        "Claude Iteration 0010: TF-Target Co-pole Signal\n"
        f"Observed SV1={obs1:.4f}, SV2={obs2:.4f}"
    )
    fig.tight_layout()

    out = FIG_DIR / "claude_trrust_copole_null_vs_obs_iter0010.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def _plot_sv1_stability_matrix() -> Path:
    mat_path = ITER_0009 / "h03_sv1_spearman_mat.npy"
    mat = np.load(mat_path)
    labels = [str(i) for i in range(mat.shape[0])]

    fig, ax = plt.subplots(figsize=(6, 5.2))
    im = ax.imshow(mat, cmap="viridis", vmin=0, vmax=1)
    ax.set_title("SV1 Cross-layer Spearman Stability")
    ax.set_xlabel("Layer")
    ax.set_ylabel("Layer")
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Spearman r")

    fig.tight_layout()
    out = FIG_DIR / "claude_sv1_spearman_matrix_iter0009.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def _load_named_gene_index_map() -> tuple[list[str], dict[str, int]]:
    gene2idx: dict[str, int] = {}
    with EDGE_PATH.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            source = row["source"].strip().upper()
            target = row["target"].strip().upper()
            source_idx = int(row["source_idx"])
            target_idx = int(row["target_idx"])
            gene2idx[source] = source_idx
            gene2idx[target] = target_idx
    genes = sorted(gene2idx.keys())
    return genes, gene2idx


def _load_go_term_genes(
    gene_set: set[str], term_ids: list[str]
) -> set[str]:
    with GENE2GO_PKL.open("rb") as f:
        gene2go_raw = pickle.load(f)

    term_set = set(term_ids)
    selected: set[str] = set()
    for raw_gene, terms in gene2go_raw.items():
        gene = raw_gene.strip().upper()
        if gene not in gene_set:
            continue
        for term in terms:
            if isinstance(term, dict):
                go_id = term.get("GO_ID", term.get("go_id", ""))
            else:
                go_id = str(term)
            if go_id in term_set:
                selected.add(gene)
                break
    return selected


def _load_trrust_sets(gene_set: set[str]) -> tuple[set[str], set[str], set[str]]:
    tf_set: set[str] = set()
    target_set: set[str] = set()
    with TRRUST_PATH.open() as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            tf = parts[0].strip().upper()
            tgt = parts[1].strip().upper()
            if tf in gene_set and tgt in gene_set and tf != tgt:
                tf_set.add(tf)
                target_set.add(tgt)
    both = tf_set & target_set
    return tf_set, target_set, both


def _annotate_axis_extremes(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    genes: list[str],
    k: int = 3,
) -> None:
    idx_top_x = np.argsort(x)[-k:]
    idx_bot_x = np.argsort(x)[:k]
    idx_top_y = np.argsort(y)[-k:]
    idx_bot_y = np.argsort(y)[:k]
    idx = set(idx_top_x.tolist() + idx_bot_x.tolist() + idx_top_y.tolist() + idx_bot_y.tolist())

    for i in idx:
        ax.annotate(
            genes[i],
            (x[i], y[i]),
            textcoords="offset points",
            xytext=(4, 4),
            fontsize=7,
            alpha=0.9,
        )


def _plot_layer11_gene_map_sv1_sv2() -> Path:
    genes, gene2idx = _load_named_gene_index_map()
    gene_set = set(genes)
    gene_to_pos = {g: i for i, g in enumerate(genes)}

    emb = np.load(EMB_PATH)
    layer11 = emb[11]
    gene_idx = np.array([gene2idx[g] for g in genes], dtype=int)
    mat = layer11[gene_idx, :]
    mat_centered = mat - mat.mean(axis=0)

    u, s, _ = np.linalg.svd(mat_centered, full_matrices=False)
    sv1 = u[:, 0] * s[0]
    sv2 = u[:, 1] * s[1]

    ev_genes = _load_go_term_genes(gene_set, ["GO:0070062"])
    secretory_genes = _load_go_term_genes(gene_set, ["GO:0005788", "GO:0005615"])
    cytoskeleton_genes = _load_go_term_genes(gene_set, ["GO:0005856"])
    mito_genes = _load_go_term_genes(gene_set, ["GO:0005739"])
    tf_set, target_set, both_set = _load_trrust_sets(gene_set)
    tf_only = tf_set - both_set
    target_only = target_set - both_set

    def idxs(gset: set[str]) -> np.ndarray:
        return np.array([gene_to_pos[g] for g in sorted(gset) if g in gene_to_pos], dtype=int)

    idx_ev = idxs(ev_genes)
    idx_secretory = idxs(secretory_genes)
    idx_cyto = idxs(cytoskeleton_genes)
    idx_mito = idxs(mito_genes)
    idx_tf_only = idxs(tf_only)
    idx_tgt_only = idxs(target_only)
    idx_both = idxs(both_set)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)

    for ax in axes:
        ax.scatter(sv1, sv2, s=18, color="lightgray", alpha=0.55, label="Other genes")
        ax.axhline(0.0, color="black", linewidth=0.7, alpha=0.6)
        ax.axvline(0.0, color="black", linewidth=0.7, alpha=0.6)
        ax.set_xlabel("SV1 coordinate")
        ax.set_ylabel("SV2 coordinate")
        ax.grid(alpha=0.2)

    axes[0].set_title("Layer-11 Gene Geometry: Biological Compartments")
    if idx_ev.size:
        axes[0].scatter(sv1[idx_ev], sv2[idx_ev], s=36, color="#e41a1c", alpha=0.95, label=f"EV (n={idx_ev.size})")
    if idx_secretory.size:
        axes[0].scatter(
            sv1[idx_secretory],
            sv2[idx_secretory],
            s=34,
            color="#377eb8",
            alpha=0.9,
            marker="s",
            label=f"Secretory/ER (n={idx_secretory.size})",
        )
    if idx_cyto.size:
        axes[0].scatter(
            sv1[idx_cyto],
            sv2[idx_cyto],
            s=42,
            color="#ff7f00",
            alpha=0.9,
            marker="^",
            label=f"Cytoskeleton (n={idx_cyto.size})",
        )
    if idx_mito.size:
        axes[0].scatter(
            sv1[idx_mito],
            sv2[idx_mito],
            s=42,
            color="#4daf4a",
            alpha=0.9,
            marker="D",
            label=f"Mitochondrion (n={idx_mito.size})",
        )
    _annotate_axis_extremes(axes[0], sv1, sv2, genes, k=3)
    axes[0].legend(loc="best", fontsize=8, frameon=True)

    axes[1].set_title("Layer-11 Gene Geometry: TRRUST Regulatory Sets")
    if idx_tf_only.size:
        axes[1].scatter(
            sv1[idx_tf_only],
            sv2[idx_tf_only],
            s=45,
            color="#984ea3",
            alpha=0.9,
            marker="^",
            label=f"TF-only (n={idx_tf_only.size})",
        )
    if idx_tgt_only.size:
        axes[1].scatter(
            sv1[idx_tgt_only],
            sv2[idx_tgt_only],
            s=38,
            color="#00a7a7",
            alpha=0.9,
            marker="o",
            label=f"Target-only (n={idx_tgt_only.size})",
        )
    if idx_both.size:
        axes[1].scatter(
            sv1[idx_both],
            sv2[idx_both],
            s=65,
            color="#f781bf",
            alpha=0.95,
            marker="*",
            label=f"TF+Target (n={idx_both.size})",
        )
    _annotate_axis_extremes(axes[1], sv1, sv2, genes, k=3)
    axes[1].legend(loc="best", fontsize=8, frameon=True)

    fig.suptitle("Claude Geometry Readout: Named Gene Manifold (SV1 vs SV2, Layer 11)", fontsize=12)
    fig.tight_layout()

    out = FIG_DIR / "claude_gene_map_sv1_sv2_layer11.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def _plot_manifold_directions_detailed() -> Path:
    genes, gene2idx = _load_named_gene_index_map()
    gene_set = set(genes)
    gene_to_pos = {g: i for i, g in enumerate(genes)}

    emb = np.load(EMB_PATH)
    layer11 = emb[11]
    gene_idx = np.array([gene2idx[g] for g in genes], dtype=int)
    mat = layer11[gene_idx, :]
    mat_centered = mat - mat.mean(axis=0)

    u, s, _ = np.linalg.svd(mat_centered, full_matrices=False)
    sv1 = u[:, 0] * s[0]
    sv2 = u[:, 1] * s[1]
    sv3 = u[:, 2] * s[2]

    ev_genes = _load_go_term_genes(gene_set, ["GO:0070062"])
    secretory_genes = _load_go_term_genes(gene_set, ["GO:0005788", "GO:0005615"])
    cytoskeleton_genes = _load_go_term_genes(gene_set, ["GO:0005856"])
    mito_genes = _load_go_term_genes(gene_set, ["GO:0005739"])

    # Make direction labels biologically stable by orienting axes.
    def _mean_of_set(values: np.ndarray, gset: set[str]) -> float:
        idx = [gene_to_pos[g] for g in gset if g in gene_to_pos]
        if not idx:
            return 0.0
        return float(np.mean(values[idx]))

    if _mean_of_set(sv1, secretory_genes) < 0:
        sv1 = -sv1
    if _mean_of_set(sv2, cytoskeleton_genes) < _mean_of_set(sv2, ev_genes):
        sv2 = -sv2

    reducer = UMAP(
        n_components=2,
        n_neighbors=20,
        min_dist=0.25,
        metric="cosine",
        random_state=42,
    )
    umap_xy = reducer.fit_transform(mat_centered)
    ux = umap_xy[:, 0]
    uy = umap_xy[:, 1]

    xreg = np.column_stack([ux, uy, np.ones(len(ux))])
    beta_sv1, *_ = np.linalg.lstsq(xreg, sv1, rcond=None)
    beta_sv2, *_ = np.linalg.lstsq(xreg, sv2, rcond=None)
    grad_sv1 = beta_sv1[:2]
    grad_sv2 = beta_sv2[:2]

    def idxs(gset: set[str]) -> np.ndarray:
        return np.array([gene_to_pos[g] for g in sorted(gset) if g in gene_to_pos], dtype=int)

    idx_ev = idxs(ev_genes)
    idx_secretory = idxs(secretory_genes)
    idx_cyto = idxs(cytoskeleton_genes)
    idx_mito = idxs(mito_genes)

    categories = np.array(["Other"] * len(genes), dtype=object)
    categories[idx_mito] = "Mitochondrion"
    categories[idx_secretory] = "Secretory/ER"
    categories[idx_ev] = "Extracellular vesicle"
    categories[idx_cyto] = "Cytoskeleton"

    colors = {
        "Other": "#bdbdbd",
        "Secretory/ER": "#377eb8",
        "Extracellular vesicle": "#e41a1c",
        "Cytoskeleton": "#ff7f00",
        "Mitochondrion": "#4daf4a",
    }

    fig = plt.figure(figsize=(18, 6.5))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3, projection="3d")

    # Panel 1: SV1-SV2 with all genes mapped.
    for cat, color in colors.items():
        mask = categories == cat
        ax1.scatter(sv1[mask], sv2[mask], s=30, c=color, alpha=0.85, label=cat)
    ax1.axhline(0.0, color="black", linewidth=0.8, alpha=0.6)
    ax1.axvline(0.0, color="black", linewidth=0.8, alpha=0.6)
    ax1.set_title("Linear Manifold Plane: Genes in (SV1, SV2)")
    ax1.set_xlabel("SV1 coordinate")
    ax1.set_ylabel("SV2 coordinate")
    ax1.grid(alpha=0.2)

    # Draw explicit direction arrows.
    lx = float(np.max(np.abs(sv1)) * 0.85)
    ly = float(np.max(np.abs(sv2)) * 0.85)
    ax1.arrow(0, 0, lx, 0, width=0.0, head_width=0.6, head_length=0.6, color="#2c7fb8", alpha=0.9)
    ax1.arrow(0, 0, -lx, 0, width=0.0, head_width=0.6, head_length=0.6, color="#2c7fb8", alpha=0.9)
    ax1.arrow(0, 0, 0, ly, width=0.0, head_width=0.6, head_length=0.6, color="#d95f0e", alpha=0.9)
    ax1.arrow(0, 0, 0, -ly, width=0.0, head_width=0.6, head_length=0.6, color="#d95f0e", alpha=0.9)
    ax1.text(lx, 0.2, "+SV1 (secretory/ER tendency)", fontsize=8, color="#2c7fb8")
    ax1.text(-lx, -0.6, "-SV1", fontsize=8, color="#2c7fb8", ha="right")
    ax1.text(0.2, ly, "+SV2 (cytoskeleton tendency)", fontsize=8, color="#d95f0e")
    ax1.text(0.2, -ly, "-SV2 (EV tendency)", fontsize=8, color="#d95f0e")

    # Annotate axis-extreme genes.
    idx_ext = set(np.argsort(sv1)[-4:].tolist() + np.argsort(sv1)[:4].tolist() + np.argsort(sv2)[-4:].tolist() + np.argsort(sv2)[:4].tolist())
    for i in idx_ext:
        ax1.annotate(genes[i], (sv1[i], sv2[i]), textcoords="offset points", xytext=(4, 3), fontsize=7)

    ax1.legend(loc="best", fontsize=8, frameon=True)

    # Panel 2: nonlinear manifold (UMAP) + projected direction arrows.
    for cat, color in colors.items():
        mask = categories == cat
        ax2.scatter(ux[mask], uy[mask], s=30, c=color, alpha=0.85, label=cat)
    ax2.set_title("Nonlinear Manifold (UMAP) + Direction Projections")
    ax2.set_xlabel("UMAP-1")
    ax2.set_ylabel("UMAP-2")
    ax2.grid(alpha=0.2)

    center = np.array([float(np.mean(ux)), float(np.mean(uy))])
    scale = 0.22 * max(float(np.ptp(ux)), float(np.ptp(uy)))

    def _draw_grad(ax: plt.Axes, grad: np.ndarray, color: str, label: str, yoff: float) -> None:
        norm = float(np.linalg.norm(grad))
        if norm == 0:
            return
        g = grad / norm
        p1 = center - g * scale
        p2 = center + g * scale
        ax.annotate("", xy=p2, xytext=p1, arrowprops=dict(arrowstyle="->", lw=2.0, color=color))
        ax.text(p2[0], p2[1] + yoff, label, color=color, fontsize=8)

    _draw_grad(ax2, grad_sv1, "#2c7fb8", "increasing SV1", 0.12)
    _draw_grad(ax2, grad_sv2, "#d95f0e", "increasing SV2", -0.12)

    # Panel 3: explicit 3D manifold (SV1, SV2, SV3).
    for cat, color in colors.items():
        mask = categories == cat
        ax3.scatter(sv1[mask], sv2[mask], sv3[mask], s=22, c=color, alpha=0.82)
    ax3.set_title("3D Linear Manifold: (SV1, SV2, SV3)")
    ax3.set_xlabel("SV1")
    ax3.set_ylabel("SV2")
    ax3.set_zlabel("SV3")

    # Add axis direction arrows in 3D.
    mx = float(np.max(np.abs(sv1)) * 0.8)
    my = float(np.max(np.abs(sv2)) * 0.8)
    mz = float(np.max(np.abs(sv3)) * 0.8)
    ax3.quiver(0, 0, 0, mx, 0, 0, color="#2c7fb8", linewidth=2)
    ax3.quiver(0, 0, 0, 0, my, 0, color="#d95f0e", linewidth=2)
    ax3.quiver(0, 0, 0, 0, 0, mz, color="#6a3d9a", linewidth=2)
    ax3.text(mx, 0, 0, "+SV1", color="#2c7fb8", fontsize=8)
    ax3.text(0, my, 0, "+SV2", color="#d95f0e", fontsize=8)
    ax3.text(0, 0, mz, "+SV3", color="#6a3d9a", fontsize=8)

    fig.suptitle("Detailed Gene Manifold with Directions (Layer 11)", fontsize=13)
    fig.tight_layout()

    out = FIG_DIR / "claude_gene_manifold_detailed_layer11.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def main() -> None:
    outputs: list[Path] = []

    outputs.append(_plot_spectral_trajectory())

    sv1_signed, sv1_sig, sv1_comps, layers = _build_iter0009_sv1_matrices()
    outputs.append(
        _heatmap(
            sv1_signed,
            sv1_comps,
            layers,
            "SV1 Compartment Geometry by Layer (iter_0009)",
            "signed log2(OR)  [top=+, bottom=-]",
            FIG_DIR / "claude_sv1_signed_log2or_iter0009.png",
            cmap="coolwarm",
            center_zero=True,
        )
    )
    outputs.append(
        _heatmap(
            sv1_sig,
            sv1_comps,
            layers,
            "SV1 Compartment Significance by Layer (iter_0009)",
            "-log10(empirical p)",
            FIG_DIR / "claude_sv1_significance_iter0009.png",
            cmap="magma",
            center_zero=False,
        )
    )

    sv2_signed, sv2_sig, sv2_comps, layers2 = _build_iter0010_sv2_matrices()
    outputs.append(
        _heatmap(
            sv2_signed,
            sv2_comps,
            layers2,
            "SV2 Compartment Geometry by Layer (iter_0010)",
            "signed log2(OR)  [top=+, bottom=-]",
            FIG_DIR / "claude_sv2_signed_log2or_iter0010.png",
            cmap="coolwarm",
            center_zero=True,
        )
    )
    outputs.append(
        _heatmap(
            sv2_sig,
            sv2_comps,
            layers2,
            "SV2 Compartment Significance by Layer (iter_0010)",
            "-log10(empirical p)",
            FIG_DIR / "claude_sv2_significance_iter0010.png",
            cmap="magma",
            center_zero=False,
        )
    )

    outputs.append(_plot_trrust_null_vs_obs())
    outputs.append(_plot_sv1_stability_matrix())
    outputs.append(_plot_layer11_gene_map_sv1_sv2())
    outputs.append(_plot_manifold_directions_detailed())

    print("Generated figures:")
    for p in outputs:
        print(str(p))


if __name__ == "__main__":
    main()
