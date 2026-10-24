
#!/usr/bin/env python3
"""Generate LAMMPS input decks for finite-temperature slab calculations.

Run from the repository root:
    python studies/main/slab_lammps.py

For each generated slab, the input deck performs:
  1. Energy minimization with fixed cell dimensions.
  2. NVT equilibration at the target temperature.
  3. NVT production, averaging potential energy, temperature, and
     an excess-potential-energy-based surface-energy estimate.

The bulk reference energy per atom is read from the bulk lattice-parameter
summary for the matching alloy and temperature.

This script generates input decks only; it does not launch LAMMPS.
"""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
ROOT = PROJECT_ROOT / "studies/main"

STRUCTURE_DIR = ROOT / "structures" / "slabs"
RUNS_DIR = ROOT / "runs" / "slab_calculations"

SUMMARY_CSV = ROOT / "results" / "bulk_lattice_parameters_summary.csv"

LIBRARY_FILE = ROOT / "library.meam"
ALLOY_FILE = ROOT / "FeNiAl.meam"


# ---------------------------------------------------------------------
# Simulation configuration
# ---------------------------------------------------------------------

TIMESTEP_PS = 0.001  # LAMMPS metal units: 1 fs
EQUILIBRATION_STEPS = 5_000
PRODUCTION_STEPS = 10_000

THERMO_EVERY = 100
AVERAGE_EVERY = 10
AVERAGE_REPEATS = 100
AVERAGE_START = AVERAGE_EVERY * AVERAGE_REPEATS

TARGET_TEMPERATURES = (300, 500, 700)
RANDOM_SEED_BASE = 831_000

# Convert eV/angstrom^2 to J/m^2.
EV_PER_ANGSTROM2_TO_J_PER_M2 = 16.02176634

# ASE/LAMMPS atom type order used by this project:
# type 1 = Al, type 2 = Fe, type 3 = Ni.
MASS_BY_TYPE = {
    1: 26.9815,
    2: 55.8450,
    3: 58.6934,
}

PAIR_COEFF = (
    "pair_coeff * * library.meam "
    "Fe Ni Al FeNiAl.meam Al Fe Ni"
)


# ---------------------------------------------------------------------
# Bulk reference data
# ---------------------------------------------------------------------

def read_bulk_references(path: Path) -> dict[tuple[str, int], dict]:
    """Read bulk energy and lattice-parameter references from the summary."""
    if not path.is_file():
        raise FileNotFoundError(
            f"Bulk reference summary not found: {path}. "
            "Run bulk_extract.py first."
        )

    required = {
        "alloy",
        "temperature_K",
        "mean_a_fcc_A",
        "mean_e_atom_eV_per_atom",
    }

    references = {}

    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)

        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")

        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                "Bulk summary is missing columns: "
                + ", ".join(sorted(missing))
            )

        for line_number, row in enumerate(reader, start=2):
            alloy = row["alloy"].strip()

            try:
                temperature = int(row["temperature_K"])
                lattice_constant = float(row["mean_a_fcc_A"])
                bulk_energy = float(row["mean_e_atom_eV_per_atom"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid bulk reference on CSV line {line_number}"
                ) from exc

            if lattice_constant <= 0:
                raise ValueError(
                    f"Non-positive lattice parameter on line {line_number}"
                )

            if not all(
                map(
                    lambda x: x == x and abs(x) != float("inf"),
                    (lattice_constant, bulk_energy),
                )
            ):
                raise ValueError(
                    f"Non-finite bulk reference on CSV line {line_number}"
                )

            key = (alloy, temperature)

            if key in references:
                raise ValueError(
                    f"Duplicate bulk reference for {alloy} at "
                    f"{temperature} K"
                )

            references[key] = {
                "alloy": alloy,
                "temperature_K": temperature,
                "lattice_constant_A": lattice_constant,
                "bulk_energy_eV_per_atom": bulk_energy,
                "source_row": dict(row),
            }

    if not references:
        raise ValueError(f"No bulk references found in {path}")

    return references


# ---------------------------------------------------------------------
# Input generation
# ---------------------------------------------------------------------

def make_input(
    *,
    alloy: str,
    temperature: int,
    seed: int,
    miller_indices: tuple[int, int, int],
    surface_area: float,
    bulk_energy: float,
) -> str:
    """Return a LAMMPS input deck for one slab configuration."""

    hkl_code = int("".join(map(str, miller_indices)))
    alloy_code = sum(ord(character) for character in alloy)

    velocity_seed = (
        RANDOM_SEED_BASE
        + seed % 100_000
        + temperature * 10
        + hkl_code
        + alloy_code
    )

    return f"""# Generated by studies/main/slab_lammps.py
# Alloy: {alloy}
# Target temperature: {temperature} K
# Miller indices: {miller_indices}
# Parent configuration seed: {seed}

clear
units metal
dimension 3
boundary p p f
atom_style atomic
newton on

read_data structure.data

mass 1 {MASS_BY_TYPE[1]}
mass 2 {MASS_BY_TYPE[2]}
mass 3 {MASS_BY_TYPE[3]}

pair_style meam
{PAIR_COEFF}

neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes

thermo {THERMO_EVERY}
thermo_style custom step time temp pe etotal atoms lx ly lz
thermo_modify norm no flush yes

# Minimize atomic positions while retaining the slab cell and vacuum.
min_style cg
minimize 1.0e-8 1.0e-6 2000 20000

# The slab cell has vacuum along z. Keep the box fixed and use NVT.
timestep {TIMESTEP_PS}
velocity all create {temperature} {velocity_seed} mom yes rot no dist gaussian

fix EQ all nvt temp {temperature} {temperature} 0.1
run {EQUILIBRATION_STEPS}
unfix EQ

# Reference energy per atom is the mean bulk potential energy
# for this alloy and temperature.
variable e_bulk equal {bulk_energy:.16g}
variable area equal {surface_area:.16g}

variable e_atom equal pe/atoms
variable temp_now equal temp

# Two exposed surfaces are present in the slab.
# gamma = (E_slab - N * e_bulk) / (2 * A)
# The result is an excess-potential-energy estimate, not a full
# finite-temperature surface free energy.
variable gamma_eV_A2 equal (pe - atoms*v_e_bulk)/(2.0*v_area)
variable gamma_J_m2 equal v_gamma_eV_A2*{EV_PER_ANGSTROM2_TO_J_PER_M2:.12g}

fix PROD all nvt temp {temperature} {temperature} 0.1

fix AVG all ave/time {AVERAGE_EVERY} {AVERAGE_REPEATS} {AVERAGE_START} &
    v_e_atom v_temp_now v_gamma_eV_A2 v_gamma_J_m2 &
    file slab_averages.dat

run {PRODUCTION_STEPS}

write_data final.data
write_restart final.restart

unfix AVG
unfix PROD

print "Completed slab NVT equilibration/production at {temperature} K."
"""


# ---------------------------------------------------------------------
# Job generation
# ---------------------------------------------------------------------

def main() -> None:
    if not STRUCTURE_DIR.is_dir():
        raise FileNotFoundError(
            f"{STRUCTURE_DIR} not found. Run slab_struct.py first."
        )

    if not LIBRARY_FILE.is_file() or not ALLOY_FILE.is_file():
        raise FileNotFoundError(
            "Expected library.meam and FeNiAl.meam under "
            f"{ROOT}."
        )

    references = read_bulk_references(SUMMARY_CSV)

    structures = sorted(STRUCTURE_DIR.rglob("*.data"))
    if not structures:
        raise FileNotFoundError(
            f"No slab *.data files found under {STRUCTURE_DIR}"
        )

    # The summary's completeness flag is deliberately not used here:
    # this module generates jobs from available structures and references.
    # Workflow-level readiness checks belong in the orchestration script.

    RUNS_DIR.mkdir(parents=True, exist_ok=True)

    manifest: list[dict] = []
    missing_references: list[str] = []

    for structure_path in structures:
        metadata_path = structure_path.with_suffix(".json")

        if not metadata_path.is_file():
            raise FileNotFoundError(
                f"Missing metadata for slab: {structure_path}"
            )

        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

        alloy = metadata.get("alloy")
        temperature = int(metadata.get("temperature_K", -1))
        seed = int(metadata.get("parent_random_seed", -1))
        hkl = tuple(metadata.get("miller_indices", ()))

        if not alloy or temperature < 0 or seed < 0:
            raise ValueError(
                f"Incomplete slab metadata: {metadata_path}"
            )

        if len(hkl) != 3:
            raise ValueError(
                f"Invalid Miller indices in {metadata_path}: {hkl}"
            )

        hkl = tuple(int(value) for value in hkl)

        if temperature not in TARGET_TEMPERATURES:
            raise ValueError(
                f"Unexpected temperature {temperature} K in "
                f"{metadata_path}"
            )

        reference_key = (alloy, temperature)
        if reference_key not in references:
            missing_references.append(
                f"{alloy} at {temperature} K ({structure_path})"
            )
            continue

        reference = references[reference_key]

        report = metadata.get("report", {})
        surface_area = report.get("surface_area_angstrom2")

        if surface_area is None:
            raise ValueError(
                f"Missing report.surface_area_angstrom2 in "
                f"{metadata_path}"
            )

        surface_area = float(surface_area)
        if surface_area <= 0:
            raise ValueError(
                f"Invalid surface area in {metadata_path}: {surface_area}"
            )

        configuration = metadata.get("configuration", structure_path.stem)
        hkl_str = "".join(map(str, hkl))

        run_name = (
            f"{alloy}_T{temperature}_{configuration}_{hkl_str}"
        )
        run_dir = RUNS_DIR / run_name
        run_dir.mkdir(parents=True, exist_ok=True)

        shutil.copy2(structure_path, run_dir / "structure.data")
        shutil.copy2(LIBRARY_FILE, run_dir / LIBRARY_FILE.name)
        shutil.copy2(ALLOY_FILE, run_dir / ALLOY_FILE.name)
        shutil.copy2(metadata_path, run_dir / "structure.json")

        # Preserve the exact bulk reference used to calculate gamma.
        (run_dir / "bulk_reference.json").write_text(
            json.dumps(reference, indent=2) + "\n",
            encoding="utf-8",
        )

        input_path = run_dir / "in.lammps"
        input_path.write_text(
            make_input(
                alloy=alloy,
                temperature=temperature,
                seed=seed,
                miller_indices=hkl,
                surface_area=surface_area,
                bulk_energy=reference["bulk_energy_eV_per_atom"],
            ),
            encoding="utf-8",
        )

        manifest.append({
            "run": run_name,
            "alloy": alloy,
            "temperature_K": temperature,
            "configuration": configuration,
            "parent_seed": seed,
            "miller_indices": hkl,
            "surface_area_angstrom2": surface_area,
            "bulk_reference_energy_eV_per_atom": (
                reference["bulk_energy_eV_per_atom"]
            ),
            "directory": str(run_dir),
            "input": str(input_path),
        })

        print(f"Wrote {input_path}")

    manifest_path = RUNS_DIR / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    print("\nSlab job generation complete.")
    print(f"  Jobs generated: {len(manifest)}")
    print(f"  Missing bulk references: {len(missing_references)}")
    print(f"  Run directory: {RUNS_DIR}")
    print(f"  Manifest: {manifest_path}")

    if missing_references:
        print("\nSkipped slabs with missing bulk references:")
        for item in missing_references:
            print(f"  {item}")

    print("\nThese are input decks only; execute them with your LAMMPS runner.")


if __name__ == "__main__":
    main()
