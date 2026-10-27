"""Shared figure style: validated categorical order (light surface), thin marks, recessive axes."""
from pathlib import Path
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

R3 = Path(__file__).resolve().parents[2] / "revision3_experiments"
OUT = Path(__file__).resolve().parents[1] / "figures"
CTX = ["lung", "lung_external", "kidney_B", "immune", "kidney_A"]
LABEL = {"lung": "lung", "lung_external": "lung (ext.)", "kidney_B": "kidney B", "immune": "immune", "kidney_A": "kidney A"}
COLOR = dict(zip(CTX, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]))
STYLE = {"lung": "-", "lung_external": "-", "kidney_B": "-", "immune": "--", "kidney_A": ":"}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"

mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.5, "lines.linewidth": 1.5, "lines.markersize": 4, "legend.frameon": False,
    "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.03, "figure.dpi": 150})


def panel(ax, letter):
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")


def end_label(ax, x, y, text, color, dy=0.0):
    ax.annotate(text, (x, y), xytext=(4, dy), textcoords="offset points", va="center", ha="left",
                fontsize=6.5, color=INK2)
