#!/usr/bin/env python3
"""Step 00 — resolve the ligand set against PubChem and download 3D conformers.

`config/ligands.csv` is authored with names only. This step resolves each name to
a CID and an InChIKey, writes them back into the config (so no identifier is ever
typed from memory), and stores the isomeric SMILES plus a 3D conformer.

The isomeric SMILES matters here: taxane stereochemistry is complex and the
generated 3D structures must later be checked against crystallographic TA1.

Outputs
  data/raw/ligands/<name>.sdf          3D conformer (2D only if PubChem has no 3D)
  data/raw/ligands/<name>.smi          isomeric SMILES
  data/raw/ligands/properties.csv      CID, InChIKey, formula, MW, heavy atoms, rotatable bonds
  config/ligands.csv                   pubchem_cid and inchikey filled in

Usage
  python scripts/00_fetch_ligands.py
  python scripts/00_fetch_ligands.py --only paclitaxel --force
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
from pathlib import Path

from _common import http_get, repo_root, setup_logging

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# Properties requested from PubChem. The isomeric-SMILES property was renamed
# from `IsomericSMILES` to `SMILES`; try the current name first and fall back.
BASE_PROPERTIES = [
    "MolecularFormula",
    "MolecularWeight",
    "InChIKey",
    "InChI",
    "HeavyAtomCount",
    "RotatableBondCount",
    "HBondDonorCount",
    "HBondAcceptorCount",
    "XLogP",
    "TPSA",
    "Charge",
]
SMILES_PROPERTY_CANDIDATES = ["SMILES", "IsomericSMILES"]


def resolve_cid(name: str, logger) -> int | None:
    url = f"{PUG}/compound/name/{name}/cids/JSON"
    response = http_get(url, logger=logger, allow_404=True)
    if response is None:
        logger.error("%s: PubChem returned no compound for this name", name)
        return None
    cids = response.json().get("IdentifierList", {}).get("CID", [])
    if not cids:
        logger.error("%s: empty CID list", name)
        return None
    if len(cids) > 1:
        logger.warning(
            "%s: name resolves to %d CIDs %s — using the first, verify it is the intended compound",
            name,
            len(cids),
            cids[:5],
        )
    return int(cids[0])


def fetch_properties(cid: int, logger) -> dict:
    """Query the property table, tolerating the SMILES property rename."""
    last_error: Exception | None = None
    for smiles_property in SMILES_PROPERTY_CANDIDATES:
        properties = BASE_PROPERTIES + [smiles_property]
        url = f"{PUG}/compound/cid/{cid}/property/{','.join(properties)}/JSON"
        try:
            response = http_get(url, logger=logger, allow_404=True)
        except RuntimeError as exc:  # property name not accepted by this server
            last_error = exc
            logger.debug("property set with %s rejected for CID %d", smiles_property, cid)
            continue
        if response is None:
            continue
        rows = response.json().get("PropertyTable", {}).get("Properties", [])
        if rows:
            record = dict(rows[0])
            # Normalize whichever SMILES key came back.
            record["isomeric_smiles"] = record.get(smiles_property, "")
            record["smiles_property_used"] = smiles_property
            return record
    raise RuntimeError(f"could not retrieve properties for CID {cid}") from last_error


def sdf_dimensionality(path: Path) -> str:
    """'3d' or '2d' for an existing SDF, decided from the z coordinates.

    Needed so that a re-run reports the same conformer source as the original
    download instead of losing the information to a 'skipped'.
    """
    with path.open() as handle:
        lines = handle.read().splitlines()
    if len(lines) < 4:
        return "unknown"
    try:
        n_atoms = int(lines[3][0:3])
    except ValueError:
        return "unknown"
    for line in lines[4 : 4 + n_atoms]:
        try:
            z = float(line[20:30])
        except ValueError:
            return "unknown"
        if abs(z) > 1e-4:
            return "3d"
    return "2d"


def fetch_structure(cid: int, out_path: Path, force: bool, logger) -> str:
    """Download an SDF conformer. Returns '3d', '2d' or 'unknown'."""
    if out_path.exists() and not force:
        existing = sdf_dimensionality(out_path)
        logger.info("  %s already present (%s), skipping download", out_path.name, existing)
        return existing
    for record_type in ("3d", "2d"):
        url = f"{PUG}/compound/cid/{cid}/SDF?record_type={record_type}"
        response = http_get(url, logger=logger, allow_404=True)
        if response is None:
            logger.warning("  no %s record for CID %d", record_type.upper(), cid)
            continue
        out_path.write_bytes(response.content)
        logger.info("  wrote %s (%s conformer, %d bytes)", out_path.name, record_type, len(response.content))
        return record_type
    raise RuntimeError(f"no SDF record available for CID {cid}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ligands", default="config/ligands.csv")
    parser.add_argument("--outdir", default="data/raw/ligands")
    parser.add_argument("--only", nargs="*", metavar="NAME", help="restrict to these ligand names")
    parser.add_argument("--force", action="store_true", help="re-download structures and overwrite stored identifiers")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    logger = setup_logging("00_fetch_ligands", verbose=args.verbose)
    root = repo_root()
    ligands_csv = root / args.ligands if not Path(args.ligands).is_absolute() else Path(args.ligands)
    out_dir = root / args.outdir if not Path(args.outdir).is_absolute() else Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with ligands_csv.open() as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    if not rows:
        logger.error("%s has no ligand rows", ligands_csv)
        return 2

    selected = {n.lower() for n in args.only} if args.only else None
    properties_table: list[dict] = []
    problems: list[str] = []
    needs_local_3d: list[str] = []

    for row in rows:
        name = row["name"].strip()
        if selected is not None and name.lower() not in selected:
            continue
        logger.info("=== %s (%s) ===", name, row.get("role", ""))

        cid = resolve_cid(name, logger)
        if cid is None:
            problems.append(f"{name}: could not resolve a CID")
            continue

        properties = fetch_properties(cid, logger)
        inchikey = properties.get("InChIKey", "")
        smiles = properties.get("isomeric_smiles", "")
        logger.info("  CID %d  InChIKey %s", cid, inchikey)
        logger.info(
            "  %s  MW %s  heavy atoms %s  rotatable bonds %s",
            properties.get("MolecularFormula", "?"),
            properties.get("MolecularWeight", "?"),
            properties.get("HeavyAtomCount", "?"),
            properties.get("RotatableBondCount", "?"),
        )

        # Never silently change an identifier that is already recorded.
        for column, resolved in (("pubchem_cid", str(cid)), ("inchikey", inchikey)):
            stored = (row.get(column) or "").strip()
            if stored and stored != resolved and not args.force:
                message = (
                    f"{name}: stored {column} {stored!r} disagrees with PubChem {resolved!r}; "
                    "keeping the stored value (use --force to overwrite)"
                )
                logger.error(message)
                problems.append(message)
                continue
            row[column] = resolved

        if smiles:
            (out_dir / f"{name}.smi").write_text(f"{smiles}\t{name}\n")
        else:
            problems.append(f"{name}: PubChem returned no isomeric SMILES")

        try:
            conformer = fetch_structure(cid, out_dir / f"{name}.sdf", args.force, logger)
        except RuntimeError as exc:
            logger.error("  %s", exc)
            problems.append(f"{name}: {exc}")
            conformer = "failed"
        if conformer == "2d":
            # Expected for the large, flexible ligands: PubChem publishes no 3D conformer for
            # them. Not a failure of this step — it is a requirement passed on to step 10, which
            # must embed them with RDKit. Recorded in the `conformer_source` column.
            needs_local_3d.append(name)

        properties_table.append(
            {
                "name": name,
                "class": row.get("class", ""),
                "role": row.get("role", ""),
                "pubchem_cid": cid,
                "inchikey": inchikey,
                "molecular_formula": properties.get("MolecularFormula", ""),
                "molecular_weight": properties.get("MolecularWeight", ""),
                "heavy_atom_count": properties.get("HeavyAtomCount", ""),
                "rotatable_bond_count": properties.get("RotatableBondCount", ""),
                "hbond_donor_count": properties.get("HBondDonorCount", ""),
                "hbond_acceptor_count": properties.get("HBondAcceptorCount", ""),
                "xlogp": properties.get("XLogP", ""),
                "tpsa": properties.get("TPSA", ""),
                "formal_charge": properties.get("Charge", ""),
                "conformer_source": conformer,
                "smiles_property_used": properties.get("smiles_property_used", ""),
                "isomeric_smiles": smiles,
                "inchi": properties.get("InChI", ""),
                "fetched_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            }
        )

    # Write identifiers back into the config, preserving column order.
    with ligands_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    logger.info("identifiers written back to %s", ligands_csv.relative_to(root))

    if properties_table:
        properties_path = out_dir / "properties.csv"
        with properties_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(properties_table[0]))
            writer.writeheader()
            writer.writerows(properties_table)
        logger.info("property table written to %s", properties_path.relative_to(root))

    if needs_local_3d:
        logger.warning(
            "no 3D conformer on PubChem for %d ligand(s): %s", len(needs_local_3d), ", ".join(needs_local_3d)
        )
        logger.warning("  step 10 must embed these with RDKit from the stored isomeric SMILES")

    if problems:
        logger.error("%d issue(s):", len(problems))
        for problem in problems:
            logger.error("  - %s", problem)
        return 1

    logger.info("all %d ligand(s) resolved and downloaded", len(properties_table))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
