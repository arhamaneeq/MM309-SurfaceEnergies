
import csv
import json
import math
from dataclasses import asdict
from pathlib import Path

from src.structures.bulk import (
    generate_bulk,
    write_bulk,
    composition_report,
)
from src.structures.slabs import (
    generate_slab,
    write_slab,
    slab_report,
)


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

STRUCTURE_ROOT = PROJECT_ROOT / "studies/main/structures"
BULK_DIR = STRUCTURE_ROOT / "bulk_temperature"
SLAB_DIR = STRUCTURE_ROOT / "slabs"

SUMMARY_CSV = (
    PROJECT_ROOT
    / "studies/main/results/bulk_lattice_parameters_summary.csv"
)


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

COMPOSITIONS = {
    "FeNiAl": {
        "Fe": 1 / 3,
        "Ni": 1 / 3,
        "Al": 1 / 3,
    },
    "Fe20Ni50Al30": {
        "Fe": 0.20,
        "Ni": 0.50,
        "Al": 0.30,
    },
    "Fe50Ni30Al20": {
        "Fe": 0.50,
        "Ni": 0.30,
        "Al": 0.20,
    },
}

SEEDS = (20261024, 20261025, 20261026)
MILLER_INDICES = ((1, 1, 1), (1, 1, 0), (1, 0, 0))

BULK_SIZE = (4, 4, 4)
SLAB_LAYERS = 12
VACUUM = 10.0

REQUIRED_COLUMNS = {
    "alloy",
    "temperature_K",
    "mean_a_fcc_A",
}


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def write_metadata(path: Path, metadata: dict) -> None:
    """Write JSON metadata alongside a generated structure."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)


def read_lattice_parameters(path: Path) -> list[dict]:
    """Read and validate alloy-temperature lattice parameters."""
    if not path.is_file():
        raise FileNotFoundError(
            f"Bulk lattice-parameter summary not found: {path}"
        )

    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)

        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")

        missing_columns = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing_columns:
            raise ValueError(
                "Bulk summary is missing required columns: "
                + ", ".join(sorted(missing_columns))
            )

        records = []
        seen = set()

        for line_number, row in enumerate(reader, start=2):
            alloy = row["alloy"].strip()

            if alloy not in COMPOSITIONS:
                raise ValueError(
                    f"Unknown alloy {alloy!r} on CSV line {line_number}"
                )

            try:
                temperature = int(row["temperature_K"])
                lattice_constant = float(row["mean_a_fcc_A"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid temperature or lattice parameter "
                    f"on CSV line {line_number}"
                ) from exc

            if not math.isfinite(lattice_constant) or lattice_constant <= 0:
                raise ValueError(
                    f"Invalid lattice parameter on CSV line {line_number}: "
                    f"{lattice_constant}"
                )

            key = (alloy, temperature)
            if key in seen:
                raise ValueError(
                    f"Duplicate alloy-temperature row: {key}"
                )

            seen.add(key)

            records.append({
                "alloy": alloy,
                "temperature_K": temperature,
                "lattice_constant_A": lattice_constant,
                "source_row": dict(row),
            })

    if not records:
        raise ValueError(f"No data rows found in {path}")

    expected = {
        (alloy, temperature)
        for alloy in COMPOSITIONS
        for temperature in (300, 500, 700)
    }
    missing = sorted(expected - seen)

    if missing:
        print(
            "\nWARNING: The CSV does not contain all expected "
            "alloy-temperature combinations."
        )
        for alloy, temperature in missing:
            print(f"  Missing: {alloy}, {temperature} K")
        print(
            f"  Generating only the {len(records)} combinations "
            "present in the CSV.\n"
        )

    return records


def print_composition(report) -> None:
    """Print the composition report returned by composition_report()."""
    for element, data in asdict(report).items():
        print(
            f"  {element}: {data['count']:3d} atoms "
            f"({100 * data['fraction']:7.3f} at.%)"
        )


# ---------------------------------------------------------------------
# Generation
# ---------------------------------------------------------------------

def main() -> None:
    records = read_lattice_parameters(SUMMARY_CSV)

    BULK_DIR.mkdir(parents=True, exist_ok=True)
    SLAB_DIR.mkdir(parents=True, exist_ok=True)

    n_bulk = 0
    n_slabs = 0

    for record in records:
        alloy = record["alloy"]
        temperature = record["temperature_K"]
        lattice_constant = record["lattice_constant_A"]
        composition = COMPOSITIONS[alloy]

        print("\n" + "=" * 72)
        print(f"Alloy: {alloy}")
        print(f"Target temperature: {temperature} K")
        print(f"Mean bulk lattice parameter: {lattice_constant:.6f} Å")
        print("=" * 72)

        for config_index, seed in enumerate(SEEDS, start=1):
            config_name = f"config{config_index}"
            temperature_name = f"T{temperature}"

            # ---------------------------------------------------------
            # Generate the parent bulk configuration.
            # ---------------------------------------------------------

            bulk_dir = BULK_DIR / alloy / temperature_name
            bulk_path = bulk_dir / f"{config_name}.data"

            atoms = generate_bulk(
                composition=composition,
                size=BULK_SIZE,
                lattice_constant=lattice_constant,
                seed=seed,
            )

            bulk_report = composition_report(atoms)
            write_bulk(atoms, bulk_path)

            write_metadata(
                bulk_path.with_suffix(".json"),
                {
                    "structure_type": "bulk",
                    "structure_file": str(bulk_path),
                    "alloy": alloy,
                    "temperature_K": temperature,
                    "configuration": config_name,
                    "composition_target": composition,
                    "random_seed": seed,
                    "lattice_constant_angstrom": lattice_constant,
                    "lattice_parameter_source": str(SUMMARY_CSV),
                    "supercell_repeats": BULK_SIZE,
                    "report": asdict(bulk_report),
                },
            )

            n_bulk += 1

            print(f"\n{bulk_path}")
            print(f"  Atoms: {len(atoms)}")
            print(f"  Seed: {seed}")
            print_composition(bulk_report)

            # ---------------------------------------------------------
            # Generate all slab orientations from this exact bulk.
            # ---------------------------------------------------------

            for hkl in MILLER_INDICES:
                hkl_str = "".join(map(str, hkl))
                slab_name = f"{config_name}_{hkl_str}"

                slab_dir = SLAB_DIR / alloy / temperature_name
                slab_path = slab_dir / f"{slab_name}.data"

                slab = generate_slab(
                    bulk_atoms=atoms,
                    miller_indices=hkl,
                    layers=SLAB_LAYERS,
                    vacuum=VACUUM,
                )

                report = slab_report(slab)
                write_slab(slab, slab_path)

                write_metadata(
                    slab_path.with_suffix(".json"),
                    {
                        "structure_type": "slab",
                        "structure_file": str(slab_path),
                        "alloy": alloy,
                        "temperature_K": temperature,
                        "configuration": config_name,
                        "parent_bulk_file": str(bulk_path),
                        "parent_random_seed": seed,
                        "composition_target": composition,
                        "lattice_constant_angstrom": lattice_constant,
                        "lattice_parameter_source": str(SUMMARY_CSV),
                        "miller_indices": hkl,
                        "layers": SLAB_LAYERS,
                        "vacuum_angstrom": VACUUM,
                        "report": asdict(report),
                    },
                )

                n_slabs += 1

                print(f"\n{slab_path}")
                print(f"  Miller indices: {hkl}")
                print(f"  Atoms: {report.n_atoms}")
                print(
                    "  Composition: "
                    + ", ".join(
                        f"{element}={data['fraction']:.4f}"
                        for element, data in asdict(
                            report.composition
                        ).items()
                    )
                )
                print(
                    f"  Surface area: "
                    f"{report.surface_area_angstrom2:.3f} Å²"
                )
                print(
                    f"  Atomic thickness: "
                    f"{report.atomic_thickness_angstrom:.3f} Å"
                )
                print(
                    f"  Cell height: "
                    f"{report.cell_height_normal_angstrom:.3f} Å"
                )

    print("\nGeneration complete.")
    print(f"  Bulk structures: {n_bulk}")
    print(f"  Slab structures: {n_slabs}")
    print(f"  Total structures: {n_bulk + n_slabs}")


if __name__ == "__main__":
    main()
