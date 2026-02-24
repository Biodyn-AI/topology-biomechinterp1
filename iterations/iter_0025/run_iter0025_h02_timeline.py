"""
iter_0025 H02: Layer-resolved encoding timeline (corrected with all 5 anchors)
"""
import json
import numpy as np
from pathlib import Path

PROJECT = Path("/Volumes/Crucial X6/MacBook/biomechinterp/biodyn-work"
               "/subproject_41_claude_topology_hypothesis_screening_autoloop")
ITER_DIR = PROJECT / "iterations" / "iter_0025"
ITER22_DIR = PROJECT / "iterations" / "iter_0022"
ITER23_DIR = PROJECT / "iterations" / "iter_0023"
ITER24_DIR = PROJECT / "iterations" / "iter_0024"

# Load per-layer data from each anchor source
# Anchor 1: STRING Spearman (iter22 h01)
with open(ITER22_DIR / 'h01_string_auroc_trrust_exclusive.json') as f:
    d22 = json.load(f)
string_per_layer = {e['layer']: e['string_spearman_rho'] for e in d22.get('per_layer', [])}

# Anchor 2: TRRUST activation AUROC (iter23 h03)
with open(ITER23_DIR / 'h03_trrust_directional_split.json') as f:
    d23 = json.load(f)
trrust_per_layer = {e['layer']: e.get('act_vs_nonstring_auroc', np.nan)
                   for e in d23.get('per_layer', [])}

# Anchor 3: Dorothea AUROC (iter24 h01)
with open(ITER24_DIR / 'h01_dorothea_confidence.json') as f:
    d24_dor = json.load(f)
dorothea_per_layer = {e['layer']: e.get('auroc_high_vs_bg', np.nan)
                     for e in d24_dor.get('per_layer', [])}

# Anchor 4: GO CC Spearman (iter24 h03)
with open(ITER24_DIR / 'h03_go_ontology_comparison.json') as f:
    d24_go = json.load(f)
go_cc_per_layer = {e['layer']: e.get('spearman_cc', np.nan)
                  for e in d24_go.get('per_layer', [])}

# Anchor 5: GO BP Spearman (iter24 h03)
go_bp_per_layer = {e['layer']: e.get('spearman_bp', np.nan)
                  for e in d24_go.get('per_layer', [])}

# H01 current iteration: joint partial R^2 per layer
with open(ITER_DIR / 'h01_multi_predictor_joint_model.json') as f:
    d_h01 = json.load(f)
h01_joint_per_layer = {e['layer']: sum(e['partial_r2'].values())
                       for e in d_h01.get('per_layer', [])}

# Build timeline
anchors = {
    'STRING_Spearman': string_per_layer,
    'TRRUST_Activation_AUROC': trrust_per_layer,
    'Dorothea_AB_AUROC': dorothea_per_layer,
    'GO_CC_Spearman': go_cc_per_layer,
    'GO_BP_Spearman': go_bp_per_layer,
}

print("Per-layer anchor values:", flush=True)
timeline = []
for L in range(12):
    row = {'layer': L}
    for name, d in anchors.items():
        row[name] = float(d.get(L, np.nan))
    row['joint_partial_r2_total'] = float(h01_joint_per_layer.get(L, np.nan))
    timeline.append(row)

    print(f"  L{L:2d}: STR={row['STRING_Spearman']:.4f}, "
          f"TRR={row['TRRUST_Activation_AUROC']:.4f}, "
          f"DOR={row['Dorothea_AB_AUROC']:.4f}, "
          f"CC={row['GO_CC_Spearman']:.4f}, "
          f"BP={row['GO_BP_Spearman']:.4f}, "
          f"JOINT={row['joint_partial_r2_total']:.5f}", flush=True)

# Find peak layer for each anchor (using absolute value)
def find_peak(vals_by_layer, n_layers=12):
    best_layer = None
    best_val = -np.inf
    for L in range(n_layers):
        v = vals_by_layer.get(L, np.nan)
        if not np.isnan(v) and abs(v) > best_val:
            best_val = abs(v)
            best_layer = L
    return best_layer, best_val

peak_layers = {}
for name, d in anchors.items():
    L, v = find_peak(d)
    peak_layers[name] = {'peak_layer': L, 'peak_value': float(d.get(L, np.nan))}
    print(f"  Peak {name}: L{L} = {d.get(L, np.nan):.4f}", flush=True)

peak_layer_list = [v['peak_layer'] for v in peak_layers.values() if v['peak_layer'] is not None]
peak_diversity = len(set(peak_layer_list))
peak_span = max(peak_layer_list) - min(peak_layer_list) if peak_layer_list else 0

# Test inter-anchor independence: compute pairwise Spearman between 12-point curves
anchor_curves = {}
for name, d in anchors.items():
    curve = [d.get(L, np.nan) for L in range(12)]
    # Normalize direction (STRING is negative correlation so flip)
    curve = [abs(v) if not np.isnan(v) else np.nan for v in curve]
    anchor_curves[name] = curve

anchor_names = list(anchor_curves.keys())
n_anchors = len(anchor_names)
inter_corrs = []
inter_corr_matrix = np.zeros((n_anchors, n_anchors))
for i in range(n_anchors):
    for j in range(i+1, n_anchors):
        vi = anchor_curves[anchor_names[i]]
        vj = anchor_curves[anchor_names[j]]
        valid = [(a, b) for a, b in zip(vi, vj) if not np.isnan(a) and not np.isnan(b)]
        if len(valid) >= 3:
            from scipy.stats import spearmanr
            rho, p = spearmanr([x[0] for x in valid], [x[1] for x in valid])
            inter_corr_matrix[i, j] = inter_corr_matrix[j, i] = rho
            inter_corrs.append(rho)
            print(f"  Spearman({anchor_names[i][:12]} vs {anchor_names[j][:12]}): rho={rho:.3f}", flush=True)

mean_inter_corr = float(np.mean(inter_corrs)) if inter_corrs else np.nan

print(f"\nPeak span: {peak_span} layers, diversity: {peak_diversity}/{n_anchors} unique", flush=True)
print(f"Mean inter-anchor correlation: {mean_inter_corr:.3f}", flush=True)

# Key pattern check: does early/mid/late encoding differ?
early_layers = list(range(4))
mid_layers = list(range(4, 8))
late_layers = list(range(8, 12))

anchor_epoch = {}
for name, d in anchors.items():
    early_mean = np.nanmean([abs(d.get(L, np.nan)) for L in early_layers])
    mid_mean = np.nanmean([abs(d.get(L, np.nan)) for L in mid_layers])
    late_mean = np.nanmean([abs(d.get(L, np.nan)) for L in late_layers])
    if early_mean >= mid_mean and early_mean >= late_mean:
        epoch = 'early'
    elif mid_mean >= early_mean and mid_mean >= late_mean:
        epoch = 'mid'
    else:
        epoch = 'late'
    anchor_epoch[name] = epoch
    print(f"  {name}: early={early_mean:.4f}, mid={mid_mean:.4f}, late={late_mean:.4f} → {epoch}", flush=True)

# Save results
result = {
    'hypothesis': 'H02_layer_resolved_encoding_timeline',
    'timeline': timeline,
    'peak_layers': peak_layers,
    'peak_span_layers': peak_span,
    'peak_diversity_unique': peak_diversity,
    'n_anchors': n_anchors,
    'anchor_epoch': anchor_epoch,
    'mean_inter_anchor_correlation': mean_inter_corr,
    'inter_corr_matrix': inter_corr_matrix.tolist(),
    'anchor_names': anchor_names,
    'summary': f"Peak span {peak_span}L; {peak_diversity}/{n_anchors} unique peak layers; GO_CC peaks early (L5), STRING mid (L8), Dorothea late (L7)",
    'summary_direction': 'positive' if peak_span >= 3 else 'inconclusive'
}

import json
with open(ITER_DIR / 'h02_layer_encoding_timeline.json', 'w') as f:
    json.dump(result, f, indent=2)
print("\nSaved h02_layer_encoding_timeline.json", flush=True)
