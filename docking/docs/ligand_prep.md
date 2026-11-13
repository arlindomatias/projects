# Ligand preparation record

Companion to `docs/receptor_prep.md`. Records how each ligand went from a name in
`config/ligands.csv` to a docking-ready PDBQT, and every choice made on the way.

## Step 00 — PubChem resolution (2026-11-13)

`scripts/00_fetch_ligands.py` resolved all eight names against PubChem and wrote the CID and
InChIKey back into `config/ligands.csv`. The InChIKeys below are the resolved values, not typed
ones. Heavy-atom counts returned by PubChem agree with RDKit parsing of the stored isomeric
SMILES for all eight compounds.

| Ligand | CID | Formula | MW | Heavy atoms | Rot. bonds | PubChem 3D? |
|---|---|---|---|---|---|---|
| paclitaxel | 36314 | C47H51NO14 | 853.9 | 62 | 14 | **no** |
| docetaxel | 148124 | C43H53NO14 | 807.9 | 58 | 13 | **no** |
| cabazitaxel | 9854073 | C45H57NO14 | 835.9 | 60 | 15 | **no** |
| zosuquidar | 3036703 | C32H31F2N3O2 | 527.6 | 39 | 6 | yes |
| vincristine | 5978 | C46H56N4O10 | 825.0 | 60 | 10 | **no** |
| verapamil | 2520 | C27H38N2O4 | 454.6 | 33 | 13 | yes |
| digoxin | 2724385 | C41H64O14 | 780.9 | 55 | 7 | **no** |
| diazepam | 3016 | C16H13ClN2O | 284.74 | 20 | 1 | yes |

The size spread (284–854 Da, 20–62 heavy atoms) is exactly why raw Vina scores cannot be compared
across this set — it is a 3x range in heavy-atom count.

### Five ligands need locally generated 3D structures

PubChem publishes no 3D conformer for paclitaxel, docetaxel, cabazitaxel, vincristine or digoxin
(they exceed its conformer-generation limits). Only a 2D record was stored for those, so **step 10
must embed them with RDKit ETKDG** from the isomeric SMILES in
`data/raw/ligands/properties.csv`, then confirm the result against crystallographic `TA1`.

This is not a minor inconvenience for this project: the taxanes carry 10–11 defined stereocentres
each, and a conformer generated from scratch is exactly the thing that must be checked against the
deposited pose before any ranking is believed.

### Undefined stereocentres to resolve, with their nature

RDKit was used to count assigned vs. unassigned stereocentres in the stored SMILES. The three
taxanes and digoxin are fully specified. Three ligands are not, and they are not equivalent cases:

| Ligand | Assigned | Unassigned | What it is | Action for step 10 |
|---|---|---|---|---|
| paclitaxel / docetaxel / cabazitaxel | all | 0 | — | none |
| digoxin | all | 0 | — | none |
| vincristine | 9 | 1 | a **ring nitrogen** (N with three C neighbours) | ignore — an amine nitrogen in a ring inverts freely and is not a fixed stereocentre |
| zosuquidar | 3 | 1 | a **ring carbon** (neighbours N, C, C) | real ambiguity; pick the configuration that matches the deposited `ZQU` in 7A6F and record the choice |
| verapamil | 0 | 1 | the single benzylic quaternary carbon | verapamil is used clinically as a **racemate**; decide explicitly whether to dock one enantiomer or both, and say which |

## Still to be done in step 10

- Embed 3D for the five 2D-only ligands (RDKit ETKDG, multiple conformers, MMFF minimization).
- Verify generated taxane stereochemistry against `TA1` extracted from 6QEX.
- Resolve the zosuquidar centre against `ZQU`; decide verapamil enantiomer handling.
- Protonation at pH 7.4, plus the planned neutral variants for the basic amines
  (zosuquidar, vincristine, verapamil).
- Meeko conversion to PDBQT; record torsion-tree root and the number of active torsions, which for
  the taxanes will be near Vina's practical limit.
