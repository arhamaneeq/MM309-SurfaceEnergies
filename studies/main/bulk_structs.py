#!/usr/bin/env python3
"""Generate chemically disordered FCC bulk configurations for Fe-Ni-Al.

Run from the repository root:
    python studies/main/bulk_struct.py

This creates one 4x4x4 FCC supercell per alloy/seed at a trial lattice
parameter. The finite-temperature equilibrium lattice parameter is determined
later by bulk_lammps.py; these structures are the starting configurations.
"""

from __future__ import annotations

import json
from pathlib import Path

from src.structures.bulk import generate_bulk, write_bulk, composition_report

ROOT = Path("studies/main/structures/bulk")
SIZE = (4, 4, 4)
TRIAL_LATTICE_CONSTANT = 3.60  # Angstrom; starting value, not assumed equilibrium
SEEDS = (20261024, 20261025, 20261026)

COMPOSITIONS: dict[str, dict[str, float]] = {
    "FeNiAl": {"Fe": 1 / 3, "Ni": 1 / 3, "Al": 1 / 3},
    "Fe20Ni50Al30": {"Fe": 0.20, "Ni": 0.50, "Al": 0.30},
    "Fe50Ni30Al20": {"Fe": 0.50, "Ni": 0.30, "Al": 0.20},
}


def report_to_dict(report) -> dict:
    """Convert the dataclass/nested composition report to JSON-compatible data."""
    if hasattr(report, "__dataclass_fields__"):
        return {key: report_to_dict(value) for key, value in vars(report).items()}
    if isinstance(report, dict):
        return {key: report_to_dict(value) for key, value in report.items()}
    if isinstance(report, (tuple, list)):
        return [report_to_dict(value) for value in report]
    if isinstance(report, (str, int, float, bool)) or report is None:
        return report
    return str(report)


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    manifest: list[dict] = []

    for alloy_name, composition in COMPOSITIONS.items():
        for seed in SEEDS:
            stem = f"{alloy_name}_seed{seed}"
            atoms = generate_bulk(
                composition=composition,
                size=SIZE,
                lattice_constant=TRIAL_LATTICE_CONSTANT,
                seed=seed,
            )
            structure_path = write_bulk(atoms, ROOT / f"{stem}.data")
            metadata = {
                "alloy": alloy_name,
                "composition_target": composition,
                "seed": seed,
                "size_conventional_cells": SIZE,
                "n_atoms": len(atoms),
                "trial_lattice_constant_angstrom": TRIAL_LATTICE_CONSTANT,
                "purpose": "starting configuration for finite-temperature bulk equilibration",
                "structure_file": structure_path.name,
                "composition_report": report_to_dict(composition_report(atoms)),
            }
            metadata_path = ROOT / f"{stem}.json"
            metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")

            manifest.append({
                "alloy": alloy_name,
                "seed": seed,
                "n_atoms": len(atoms),
                "structure": str(structure_path),
                "metadata": str(metadata_path),
            })
            print(f"Wrote {structure_path} ({len(atoms)} atoms)")

    manifest_path = ROOT / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nGenerated {len(manifest)} bulk configurations.")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
