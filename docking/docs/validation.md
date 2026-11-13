# Validation gate

Self-docking (redocking) of the crystallographic ligands must pass before any cross-ligand
comparison is interpreted.

| Receptor | Ligand | Symmetry-corrected RMSD (A) | Top-cluster seeds | Status |
|---|---|---|---|---|
| 6QEX | TA1 (taxol) | pending | pending | pending |
| 7A6F | ZQU (zosuquidar) | pending | pending | pending |
| 7A69 | R1Q (vincristine) | pending | pending | pending |

Pass threshold: RMSD <= 2.5 A for the taxol reference. Failure means revising box, ligand
flexibility or scoring function — not loosening the threshold.

## Cross-docking

Dock TA1 into 7A6F and 7A65 to quantify dependence on the induced-fit snapshot. Results here.
