# Receptor preparation record

For each receptor, record: chain IDs kept vs. deleted (Fab fragments), heteroatoms removed,
missing residues/side chains and whether they were rebuilt, protonation tool and pH, and the
resulting PDBQT path.

## Step 00 — download and verification (2026-11-13)

`scripts/00_fetch_receptors.py` downloaded all seven structures from RCSB (mmCIF + legacy PDB,
both formats available for every entry) and checked each file against `config/targets.yaml`
rather than trusting the config. **All seven passed**: the resolution in
`_em_3d_reconstruction.resolution` matched the configured value exactly, and every declared
reference ligand was found among the file's chemical components.

| PDB | Resolution (Å, from file) | Deposited title | Ref. ligand | Present |
|---|---|---|---|---|
| 6QEX | 3.60 | Nanodisc reconstituted human ABCB1 in complex with UIC2 fab and taxol | TA1 | yes |
| 7A6F | 3.50 | Nanodisc reconstituted human ABCB1 in complex with MRK16 Fab and zosuquidar | ZQU | yes |
| 7A69 | 3.20 | Nanodisc reconstituted human ABCB1 in complex with MRK16 Fab and vincristine | R1Q | yes |
| 7A65 | 3.90 | Nanodisc reconstituted, drug-free human ABCB1 in complex with MRK16 Fab | — | n/a (apo) |
| 7A6C | 3.60 | Nanodisc reconstituted human ABCB1 in complex with MRK16 Fab and elacridar | R0Z | yes |
| 7A6E | 3.60 | Nanodisc reconstituted human ABCB1 in complex with MRK16 Fab and tariquidar | R1H | yes |
| 6C0V | 3.40 | Molecular structure of human P-glycoprotein in the ATP-bound, outward-facing conformation | — | n/a |

### Chain assignment (resolved, no longer a guess)

Taken from the entity table the step writes to `data/raw/pdb/entities.csv`:

| PDB | Entity 1 (keep) | Entity 2 (delete) | Entity 3 (delete) |
|---|---|---|---|
| 6QEX | Multidrug resistance protein 1 — chain **A** | UIC2 Fab light chain — chain B | UIC2 Fab heavy chain — chain C |
| 7A6F | ABCB1 — chain **A** | MRK16 Fab light chain — B | MRK16 Fab heavy chain — C |
| 7A69 | ABCB1 — chain **A** | MRK16 Fab light chain — B | MRK16 Fab heavy chain — C |
| 7A65 | ABCB1 — chain **A** | MRK16 Fab light chain — B | MRK16 Fab heavy chain — C |
| 7A6C | ABCB1 — chain **A** | "If kappa light chain" — B | MRK16 Fab heavy chain — C |
| 7A6E | ABCB1 — chain **A** | MRK16 Fab light chain — B | MRK16 Fab heavy chain — C |
| 6C0V | ABCB1 — chain **A** | — (no Fab in this deposition) | — |

So the rule is uniform: **keep chain A, delete chains B and C.** This is now recorded as
`keep_chains` / `fab_chains` in `config/targets.yaml`.

### Heteroatom inventory (observed, per structure)

- `CLR` (cholesterol) is present in **all six nanodisc structures** — it is the one heteroatom
  whose in/out status actually changes the cavity, which is why it has a dedicated flag.
- `NAG` (glycan) and `3PE` (nanodisc phospholipid) occur **only in 6QEX**.
- `ATP` + `MG` occur **only in 6C0V**, which is not docked.
- No crystallographic waters are deposited in any of these cryo-EM models.

### Still to be done in step 10

- Missing residues / side chains: **not yet assessed.** At 3.2–3.9 Å these models are expected to
  have unmodelled loops (notably the linker between the two halves). Each gap must be listed here
  with the decision taken (rebuilt with PDBFixer vs. left missing); large gaps are not modelled
  silently.
- Hydrogens and protonation at pH 7.4 (`pdb2pqr`/PROPKA), then PDBQT via Meeko.
- Residue numbering must be checked against UniProt P08183 before any cavity residue is quoted.
