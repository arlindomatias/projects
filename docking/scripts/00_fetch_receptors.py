#!/usr/bin/env python3
"""Step 00 — download the ABCB1 receptor structures listed in config/targets.yaml.

Downloads mmCIF (required) and legacy PDB (best effort) from RCSB, then verifies
each file against the config instead of trusting it:

  * resolution reported in the file vs. `resolution` in targets.yaml
  * presence of the declared `reference_ligand` chemical component
  * polymer entities and their chain IDs, so the preparation step can decide
    which chains are Fab fragments rather than guessing

Outputs
  data/raw/pdb/<PDB>.cif, <PDB>.pdb
  data/raw/pdb/manifest.csv   one row per structure (checksum, resolution, checks)
  data/raw/pdb/entities.csv   one row per entity x chain (description, type)

Usage
  python scripts/00_fetch_receptors.py
  python scripts/00_fetch_receptors.py --only 6QEX 7A6F --force
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from _common import http_get, load_yaml, repo_root, setup_logging, sha256_file

RCSB_CIF = "https://files.rcsb.org/download/{pdb}.cif"
RCSB_PDB = "https://files.rcsb.org/download/{pdb}.pdb"

# Tolerance when comparing the resolution in the file with the one in config.
RESOLUTION_TOLERANCE = 0.05


def _as_list(value) -> list[str]:
    """MMCIF2Dict gives a scalar for single-row loops and a list otherwise."""
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(item) for item in value]
    return [str(value)]


def _first_float(values: list[str]) -> float | None:
    for value in values:
        try:
            return float(value)
        except (TypeError, ValueError):
            continue
    return None


def reported_resolution(cif: dict) -> tuple[float | None, str]:
    """Resolution from the file, preferring the cryo-EM reconstruction record."""
    candidates = [
        ("_em_3d_reconstruction.resolution", "em_3d_reconstruction"),
        ("_refine.ls_d_res_high", "refine.ls_d_res_high"),
        ("_reflns.d_resolution_high", "reflns.d_resolution_high"),
    ]
    for key, source in candidates:
        value = _first_float(_as_list(cif.get(key)))
        if value is not None:
            return value, source
    return None, "not_found"


def entity_rows(pdb_id: str, cif: dict) -> list[dict]:
    """One row per entity, with its chain IDs and description."""
    ids = _as_list(cif.get("_entity.id"))
    types = _as_list(cif.get("_entity.type"))
    descriptions = _as_list(cif.get("_entity.pdbx_description"))

    # Chain IDs per polymer entity.
    poly_entity = _as_list(cif.get("_entity_poly.entity_id"))
    poly_strands = _as_list(cif.get("_entity_poly.pdbx_strand_id"))
    chains_by_entity = dict(zip(poly_entity, poly_strands))

    # Component ID per non-polymer entity (ligands, lipids, glycans).
    nonpoly_entity = _as_list(cif.get("_pdbx_entity_nonpoly.entity_id"))
    nonpoly_comp = _as_list(cif.get("_pdbx_entity_nonpoly.comp_id"))
    comp_by_entity = dict(zip(nonpoly_entity, nonpoly_comp))

    rows = []
    for index, entity_id in enumerate(ids):
        rows.append(
            {
                "pdb_id": pdb_id,
                "entity_id": entity_id,
                "entity_type": types[index] if index < len(types) else "",
                "description": descriptions[index] if index < len(descriptions) else "",
                "chains": chains_by_entity.get(entity_id, ""),
                "comp_id": comp_by_entity.get(entity_id, ""),
            }
        )
    return rows


def present_components(cif: dict) -> set[str]:
    """Every chemical component declared in the file (polymer residues included)."""
    return {comp.upper() for comp in _as_list(cif.get("_chem_comp.id"))}


def download(pdb_id: str, out_dir: Path, force: bool, logger) -> dict[str, Path | None]:
    paths: dict[str, Path | None] = {}
    for fmt, url_template, required in (
        ("cif", RCSB_CIF, True),
        ("pdb", RCSB_PDB, False),
    ):
        target = out_dir / f"{pdb_id}.{fmt}"
        if target.exists() and not force:
            logger.info("%s: %s already present, skipping download", pdb_id, target.name)
            paths[fmt] = target
            continue

        url = url_template.format(pdb=pdb_id)
        response = http_get(url, logger=logger, allow_404=not required)
        if response is None:
            logger.warning(
                "%s: no legacy PDB format available at RCSB (expected for very large models)",
                pdb_id,
            )
            paths[fmt] = None
            continue

        target.write_bytes(response.content)
        logger.info("%s: wrote %s (%.1f MB)", pdb_id, target.name, len(response.content) / 1e6)
        paths[fmt] = target
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/targets.yaml")
    parser.add_argument("--outdir", default="data/raw/pdb")
    parser.add_argument("--only", nargs="*", metavar="PDB", help="restrict to these PDB IDs")
    parser.add_argument("--force", action="store_true", help="re-download even if the file exists")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    logger = setup_logging("00_fetch_receptors", verbose=args.verbose)
    root = repo_root()
    out_dir = (root / args.outdir) if not Path(args.outdir).is_absolute() else Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    config = load_yaml(args.config)
    receptors: dict = config.get("receptors", {})
    wanted = [p.upper() for p in args.only] if args.only else list(receptors)
    missing = [p for p in wanted if p not in receptors]
    if missing:
        logger.error("not in %s: %s", args.config, ", ".join(missing))
        return 2

    manifest: list[dict] = []
    entities: list[dict] = []
    problems: list[str] = []

    for pdb_id in wanted:
        spec = receptors[pdb_id]
        logger.info("=== %s (%s) ===", pdb_id, spec.get("role", "?"))
        paths = download(pdb_id, out_dir, args.force, logger)
        cif_path = paths.get("cif")
        if cif_path is None:
            problems.append(f"{pdb_id}: mmCIF download failed")
            continue

        cif = MMCIF2Dict(str(cif_path))
        entities.extend(entity_rows(pdb_id, cif))

        file_resolution, resolution_source = reported_resolution(cif)
        config_resolution = spec.get("resolution")
        resolution_ok = (
            file_resolution is not None
            and config_resolution is not None
            and abs(file_resolution - float(config_resolution)) <= RESOLUTION_TOLERANCE
        )
        if not resolution_ok:
            message = (
                f"{pdb_id}: resolution in file {file_resolution} A ({resolution_source}) "
                f"!= config {config_resolution} A"
            )
            logger.warning(message)
            problems.append(message)
        else:
            logger.info("%s: resolution %.2f A confirmed (%s)", pdb_id, file_resolution, resolution_source)

        components = present_components(cif)
        reference_ligand = spec.get("reference_ligand")
        if reference_ligand:
            reference_ligand = str(reference_ligand).upper()
            ligand_ok = reference_ligand in components
            if ligand_ok:
                logger.info("%s: reference ligand %s present", pdb_id, reference_ligand)
            else:
                message = f"{pdb_id}: reference ligand {reference_ligand} NOT found in the file"
                logger.error(message)
                problems.append(message)
        else:
            ligand_ok = None  # apo / discussion-only structures declare no ligand

        title = "; ".join(_as_list(cif.get("_struct.title"))) or ""
        polymer_chains = sorted(
            {
                chain.strip()
                for row in entities
                if row["pdb_id"] == pdb_id and row["entity_type"] == "polymer"
                for chain in row["chains"].split(",")
                if chain.strip()
            }
        )

        manifest.append(
            {
                "pdb_id": pdb_id,
                "role": spec.get("role", ""),
                "title": title,
                "resolution_file": file_resolution if file_resolution is not None else "",
                "resolution_source": resolution_source,
                "resolution_config": config_resolution if config_resolution is not None else "",
                "resolution_matches_config": resolution_ok,
                "reference_ligand": reference_ligand or "",
                "reference_ligand_present": "" if ligand_ok is None else ligand_ok,
                "polymer_chains": ",".join(polymer_chains),
                "n_polymer_entities": sum(
                    1 for row in entities if row["pdb_id"] == pdb_id and row["entity_type"] == "polymer"
                ),
                "has_cholesterol_CLR": "CLR" in components,
                "cif_sha256": sha256_file(cif_path),
                "pdb_format_available": paths.get("pdb") is not None,
                "downloaded_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            }
        )

    if manifest:
        manifest_path = out_dir / "manifest.csv"
        with manifest_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
            writer.writeheader()
            writer.writerows(manifest)
        logger.info("manifest written to %s", manifest_path.relative_to(root))

    if entities:
        entities_path = out_dir / "entities.csv"
        with entities_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(entities[0]))
            writer.writeheader()
            writer.writerows(entities)
        logger.info("entity table written to %s", entities_path.relative_to(root))

    if problems:
        logger.error("%d check(s) failed:", len(problems))
        for problem in problems:
            logger.error("  - %s", problem)
        return 1

    logger.info("all %d structure(s) downloaded and verified", len(manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
