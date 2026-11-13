# Decision log

One line per protocol choice that could reasonably have gone the other way.

| Date | Decision | Rationale |
|---|---|---|
| 2026-11-13 | Primary receptor is 6QEX (ABCB1 + taxol) | Only deposited structure with a taxane bound; lets the protocol be validated by redocking the actual chemotype under study |
| 2026-11-13 | Only inward-facing structures are docked; 6C0V (outward-facing) is used for discussion | The outward-facing cavity is collapsed and is not a meaningful docking site |
| 2026-11-13 | Cholesterol (CLR) excluded from the receptor by default | Keeps the cavity definition simple; a `keep_cholesterol: true` run is the planned sensitivity check |
| 2026-11-13 | Three docking regimes (taxane site, modulator site, whole cavity) | P-gp is polyspecific with a ~6000 A^3 cavity; a single small box would prejudge where the taxanes bind |
| 2026-11-13 | Ranking reported as raw score + ligand efficiency + contact fingerprint | Vina scores scale with ligand size; the set spans 230-920 Da, so raw dG is not comparable across chemotypes |
| 2026-11-13 | Environment pinned to Python 3.10, not 3.11 | 3.10 is what vina + pymol-open-source + snakemake actually co-resolve to on conda-forge/bioconda; `environment.lock.yml` records the exact working set |
| 2026-11-13 | Receptor identity verified from the file, not the config | `00_fetch_receptors.py` compares the mmCIF resolution and chemical components against `targets.yaml` and fails the step on a mismatch, so a wrong PDB ID or ligand code cannot propagate |
| 2026-11-13 | Keep chain A, delete chains B and C, in every receptor | Resolved from the deposited entity descriptions (not assumed): ABCB1 is entity 1 / chain A and the Fab light and heavy chains are B and C in all six nanodisc structures |
| 2026-11-13 | 3D conformers generated locally for the 5 large ligands | PubChem has no 3D record for paclitaxel, docetaxel, cabazitaxel, vincristine or digoxin (too large/flexible), so step 10 must build them with RDKit ETKDG from the stored isomeric SMILES and verify stereochemistry against crystallographic TA1 |
| 2026-11-13 | Undefined stereocentres triaged by nature, not counted | Of three ligands with an unassigned centre, vincristine's is a freely inverting ring nitrogen (ignore), zosuquidar's is a real ring carbon (resolve against deposited ZQU) and verapamil's is the clinical racemate (explicit choice required) — see docs/ligand_prep.md |
