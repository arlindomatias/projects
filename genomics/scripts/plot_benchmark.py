"""Precision-recall curves from rtg vcfeval ROC files, one line per caller.

Panels: all variants (weighted_roc), SNVs (snp_roc), indels (non_snp_roc).
The dot at the end of each curve is the operating point without a score threshold
(the "None" row of summary.txt); the legend gives its F1.

Called from the Snakefile (`script:`), or by hand:
    python scripts/plot_benchmark.py out.png bcftools=results/eval_bcftools/HG002 \
        GATK=results/eval/HG002 DeepVariant=results/eval_deepvariant/HG002
"""

import gzip
import io
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import CALLER_COLORS, INK_2, apply_style, plt  # noqa: E402

RECALL_MIN = 0.9  # zoom: all callers operate above this recall

PANELS = [
    ("weighted_roc.tsv.gz", "All variants"),
    ("snp_roc.tsv.gz", "SNVs"),
    ("non_snp_roc.tsv.gz", "Indels"),
]


def read_roc(path):
    """Return (table, score field) from a vcfeval ROC file."""
    with gzip.open(path, "rt") as fh:
        lines = fh.readlines()
    score_field = next(
        (l.split(":", 1)[1].strip() for l in lines if l.startswith("#score field")), "?"
    )
    header = next(i for i, l in enumerate(lines) if l.startswith("#score\t"))
    cols = lines[header].lstrip("#").strip().split("\t")
    table = pd.read_csv(io.StringIO("".join(lines[header + 1:])), sep="\t", names=cols)
    return table, score_field


def plot(evals, out_png):
    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), constrained_layout=True)
    score_fields = {}
    for ax, (fname, title) in zip(axes, PANELS):
        visible = []  # precision values inside the zoom window, to set the y range
        for caller, eval_dir in evals.items():
            roc, score_fields[caller] = read_roc(os.path.join(eval_dir, fname))
            color = CALLER_COLORS[caller]
            end = roc.iloc[-1]  # lowest threshold = all calls
            ax.plot(roc["sensitivity"], roc["precision"], color=color, lw=2,
                    label=f"{caller}  (F1 {end['f_measure']:.4f})")
            ax.plot(end["sensitivity"], end["precision"], "o", ms=8, color=color,
                    markeredgecolor="white", markeredgewidth=2, zorder=3)
            visible.append(roc.loc[roc["sensitivity"] >= RECALL_MIN, "precision"])
        ymin = min(v.min() for v in visible if len(v))
        ax.set_xlim(RECALL_MIN, 1.0)
        ax.set_ylim(ymin - (1 - ymin) * 0.05, 1.0 + (1 - ymin) * 0.03)
        ax.set_title(title, loc="left")
        ax.set_xlabel("Recall (sensitivity)")
        ax.set_ylabel("Precision")
        ax.legend(loc="best", fontsize=9)
    scores = ", ".join(f"{c}: {s}" for c, s in score_fields.items())
    fig.suptitle("Precision vs recall against GIAB HG002 (chr20, high-confidence regions)",
                 x=0.01, ha="left", fontsize=12)
    fig.supxlabel(f"Curve = varying score threshold ({scores}); dot = no threshold. "
                  f"Zoomed to recall >= {RECALL_MIN}.", x=0.99, ha="right", fontsize=8,
                  color=INK_2)
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    fig.savefig(out_png, dpi=150)


if __name__ == "__main__":
    if "snakemake" in globals():
        sys.stderr = open(snakemake.log[0], "w")  # noqa: F821
        plot(dict(snakemake.params.evals), snakemake.output[0])  # noqa: F821
    else:
        out, *pairs = sys.argv[1:]
        plot(dict(p.split("=", 1) for p in pairs), out)
