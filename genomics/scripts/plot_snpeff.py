"""SnpEff summary plots: variants per impact class and the most common effects.

Each variant is counted once, by its first ANN entry (SnpEff sorts the most
severe effect first).

Called from the Snakefile (`script:`), or by hand:
    python scripts/plot_snpeff.py results/annot/HG002.snpeff.vcf.gz out.png
"""

import os
import sys
from collections import Counter

from cyvcf2 import VCF

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_style import INK, INK_2, SINGLE_HUE, apply_style, plt  # noqa: E402

IMPACTS = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
CODING = IMPACTS[:3]  # MODIFIER (non-coding) dwarfs the rest; reported in the title
TOP_EFFECTS = 10


def count(vcf_path):
    impacts, effects = Counter(), Counter()
    for var in VCF(vcf_path):
        ann = var.INFO.get("ANN")
        if not ann:
            continue
        fields = ann.split(",", 1)[0].split("|")
        impact, effect = fields[2], fields[1].split("&")[0]
        impacts[impact] += 1
        if impact != "MODIFIER":
            effects[effect] += 1
    return impacts, effects


def hbar(ax, labels, values):
    y = range(len(labels))
    ax.set_axisbelow(True)
    ax.barh(y, values, color=SINGLE_HUE, height=0.6)
    ax.set_yticks(list(y), labels)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.tick_params(axis="y", length=0, labelcolor=INK)
    xmax = max(values)
    ax.set_xlim(0, xmax * 1.15)
    for yi, v in zip(y, values):
        ax.text(v + xmax * 0.01, yi, f"{v:,}", va="center", fontsize=9, color=INK)


def plot(vcf_path, out_png):
    impacts, effects = count(vcf_path)
    apply_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True,
                                   gridspec_kw={"width_ratios": [1, 1.4]})
    hbar(ax1, CODING, [impacts.get(i, 0) for i in CODING])
    ax1.set_title(f"Variants by impact (n = {sum(impacts.values()):,}; "
                  f"MODIFIER {impacts.get('MODIFIER', 0):,} not shown)", loc="left")
    ax1.set_xlabel("Variants")

    top = effects.most_common(TOP_EFFECTS)
    hbar(ax2, [e.replace("_", " ") for e, _ in top], [n for _, n in top])
    ax2.set_title(f"Most common effects, HIGH/MODERATE/LOW (n = {sum(effects.values()):,})",
                  loc="left")
    ax2.set_xlabel("Variants")

    fig.supxlabel("Most severe effect per variant (first SnpEff ANN entry).",
                  x=0.99, ha="right", fontsize=8, color=INK_2)
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    fig.savefig(out_png, dpi=150)
    for i in IMPACTS:
        print(f"{i}\t{impacts.get(i, 0)}", file=sys.stderr)


if __name__ == "__main__":
    if "snakemake" in globals():
        sys.stderr = open(snakemake.log[0], "w")  # noqa: F821
        plot(snakemake.input.vcf, snakemake.output[0])  # noqa: F821
    else:
        plot(sys.argv[1], sys.argv[2])
