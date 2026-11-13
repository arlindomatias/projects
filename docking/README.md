# Taxanes vs. P-glycoprotein — structural bioinformatics and molecular docking

Structure-based docking study of the three clinical taxanes (paclitaxel, docetaxel, cabazitaxel)
against human P-glycoprotein (P-gp / MDR1, `ABCB1`, UniProt P08183), the efflux transporter that
limits exposure of exogenous compounds in the central nervous system at the blood–brain barrier.

The taxanes are compared against two pharmacological anchors with co-structures in the PDB:
a third-generation **inhibitor** (zosuquidar) and a well-characterized **substrate** (vincristine),
plus additional substrate anchors (verapamil, digoxin) and a CNS-penetrant non-substrate negative
control (diazepam).

## Why these structures

| PDB | Bound drug | Role here |
|---|---|---|
| 6QEX | taxol (`TA1`) | primary receptor; the taxane reference pose |
| 7A6F | zosuquidar (`ZQU`) | inhibitor reference / modulator site |
| 7A69 | vincristine (`R1Q`) | substrate reference |
| 7A65 | — | apo control |
| 7A6C / 7A6E | elacridar / tariquidar | optional inhibitors |
| 6C0V | ATP (outward-facing) | conformational discussion only |

Having a deposited **taxol–ABCB1** complex is what makes the question tractable: the docking protocol
is validated by reproducing that pose before any taxane is ranked against another.

## Scope and limits

Docking estimates binding affinity and pose, not transport. It cannot by itself separate a substrate
from an inhibitor, and it does not model the ATP-driven transport cycle. Results are reported as
affinity/pose hypotheses to be read against published transport data.

## Layout

```
config/        all parameters (targets, ligands, docking boxes)
scripts/       numbered CLI steps: 00_fetch, 10_prepare, 20_dock, 30_analyze
workflow/      Snakefile
data/          raw downloads and prepared inputs (gitignored, regenerable)
results/       poses, score tables, figures
docs/          receptor preparation, validation, decision log
notebooks/     exploratory analysis and reporting
```

## Setup

```bash
conda env create -f environment.yml     # or environment.lock.yml for the exact resolved versions
conda activate docking
```

## Running the fetch step

```bash
python scripts/00_fetch_receptors.py    # downloads and verifies the 7 structures
python scripts/00_fetch_ligands.py      # resolves the 8 ligands against PubChem
```

Both steps verify rather than trust: the receptor step fails if a file's resolution or reference
ligand disagrees with `config/targets.yaml`, and the ligand step writes the resolved CIDs and
InChIKeys back into `config/ligands.csv` so no identifier is ever taken from memory.

## Status

Protocol defined (see `CLAUDE.md`) and **step 00 complete**: all seven receptors downloaded and
verified against the deposited mmCIF, all eight ligands resolved on PubChem. See
`docs/receptor_prep.md` for the verification record.

Next: step 10 — receptor cleaning (keep chain A, drop the Fab chains B/C), local 3D generation for
the five ligands PubChem has no 3D record for, protonation at pH 7.4 — then the self-docking
validation gate.
