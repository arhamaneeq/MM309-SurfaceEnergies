
from dataclasses import asdict
import json
from pathlib import Path

from src.structures.bulk import generate_bulk, write_bulk, composition_report
from src.structures.slabs import generate_slab, write_slab, slab_report


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

ROOT = Path("studies/main/structures")
BULK_DIR = ROOT / "bulk"
SLAB_DIR = ROOT / "slabs"

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
LATTICE_CONSTANT = 3.60

SLAB_LAYERS = 12
VACUUM = 10.0


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def write_metadata(path: Path, metadata: dict) -> None:
    """Write JSON metadata alongside a generated structure."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)


# ---------------------------------------------------------------------
# Generation
# ---------------------------------------------------------------------

def main() -> None:
    BULK_DIR.mkdir(parents=True, exist_ok=True)
    SLAB_DIR.mkdir(parents=True, exist_ok=True)

    n_bulk = 0
    n_slabs = 0

    for name, composition in COMPOSITIONS.items():
        for seed in SEEDS:
            config_name = f"{name}_config{SEEDS.index(seed) + 1}"
            bulk_path = BULK_DIR / f"{config_name}.data"

            # Generate the parent bulk exactly once.
            atoms = generate_bulk(
                composition=composition,
                size=BULK_SIZE,
                lattice_constant=LATTICE_CONSTANT,
                seed=seed,
            )

            bulk_report = composition_report(atoms)

            write_bulk(atoms, bulk_path)

            write_metadata(
                bulk_path.with_suffix(".json"),
                {
                    "structure_type": "bulk",
                    "structure_file": str(bulk_path),
                    "configuration": config_name,
                    "composition_target": composition,
                    "random_seed": seed,
                    "lattice_constant_angstrom": LATTICE_CONSTANT,
                    "supercell_repeats": BULK_SIZE,
                    "report": asdict(bulk_report),
                },
            )

            n_bulk += 1

            print(f"\n{bulk_path}")
            print(f"  Atoms: {len(atoms)}")
            print(f"  Seed: {seed}")

            for element, data in asdict(bulk_report).items():
                print(
                    f"  {element}: {data['count']:3d} atoms "
                    f"({100 * data['fraction']:7.3f} at.%)"
                )

            # Generate each slab from this exact in-memory bulk.
            for hkl in MILLER_INDICES:
                hkl_str = "".join(map(str, hkl))
                slab_name = f"{config_name}_{hkl_str}"
                slab_path = SLAB_DIR / f"{slab_name}.data"

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
                        "configuration": config_name,
                        "parent_bulk_file": str(bulk_path),
                        "parent_random_seed": seed,
                        "composition_target": composition,
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
                    f"{report.surface_area_angstrom2:.3f} A^2"
                )
                print(
                    f"  Atomic thickness: "
                    f"{report.atomic_thickness_angstrom:.3f} A"
                )
                print(
                    f"  Cell height: "
                    f"{report.cell_height_normal_angstrom:.3f} A"
                )

    print("\nGeneration complete.")
    print(f"  Bulk structures: {n_bulk}")
    print(f"  Slab structures: {n_slabs}")
    print(f"  Total structures: {n_bulk + n_slabs}")


if __name__ == "__main__":
    main()
