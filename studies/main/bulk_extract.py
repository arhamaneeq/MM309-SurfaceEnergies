#!/usr/bin/env python3
"""Extract temperature-dependent equivalent FCC lattice parameters from LAMMPS.

Expected input:
  studies/main/runs/bulk_equilibration/<alloy>_seed<seed>_T<T>/bulk_averages.dat

The file is written by fix ave/time in bulk_lammps.py. With the current input
(100 100 10000), it contains one averaged row per completed production run.
This script also supports multiple rows if the averaging frequency is changed.

Outputs:
  studies/main/results/bulk_lattice_parameters_by_run.csv
  studies/main/results/bulk_lattice_parameters_summary.csv

The per-run CSV preserves all available averaged rows and records their mean
and standard deviation. The summary CSV aggregates run-level means across
independent seeds for each alloy and temperature.

Run from the repository root:
  python studies/main/extract_bulk_lattice_parameters.py
"""

from __future__ import annotations

import csv
import re
import statistics
from pathlib import Path
from typing import Iterable

ROOT = Path("studies/main")
RUNS_DIR = ROOT / "runs" / "bulk_equilibration"
OUTPUT_DIR = ROOT / "results"

AVERAGES_FILENAME = "bulk_averages.dat"
RUN_OUTPUT = OUTPUT_DIR / "bulk_lattice_parameters_by_run.csv"
SUMMARY_OUTPUT = OUTPUT_DIR / "bulk_lattice_parameters_summary.csv"

# Expected fix ave/time output columns, in the order supplied to fix AVG:
# v_a_fcc v_e_atom v_temp_now v_press_now v_volume_now
VALUE_COLUMNS = (
    "a_fcc_A",
    "e_atom_eV_per_atom",
    "measured_temperature_K",
    "pressure_bar",
    "volume_A3",
)


def read_averages(path: Path) -> list[dict[str, float]]:
    """Read numeric fix ave/time rows, ignoring blank/comment lines."""
    rows: list[dict[str, float]] = []

    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        fields = stripped.split()
        # fix ave/time output usually begins with the timestep followed by
        # one value per requested input variable.
        try:
            values = [float(field) for field in fields]
        except ValueError:
            # Ignore a possible non-numeric header, but don't silently ignore
            # malformed data rows.
            if not rows:
                continue
            raise ValueError(
                f"Non-numeric data in {path}:{line_number}: {stripped}"
            )

        if len(values) != 1 + len(VALUE_COLUMNS):
            raise ValueError(
                f"Expected timestep plus {len(VALUE_COLUMNS)} values "
                f"({1 + len(VALUE_COLUMNS)} columns), got {len(values)} "
                f"in {path}:{line_number}: {stripped}"
            )

        row = {"timestep": values[0]}
        row.update(dict(zip(
            (
                "a_fcc_A",
                "e_atom_eV_per_atom",
                "measured_temperature_K",
                "pressure_bar",
                "volume_A3",
            ),
            values[1:],
        )))
        rows.append(row)

    return rows


def parse_run_name(run_name: str) -> tuple[str, int, int] | None:
    """Parse <alloy>_seed<seed>_T<temperature>; tolerate underscores in alloy."""
    match = re.fullmatch(r"(?P<alloy>.+)_seed(?P<seed>\d+)_T(?P<temp>\d+)", run_name)
    if not match:
        return None
    return match.group("alloy"), int(match.group("seed")), int(match.group("temp"))


def mean_or_none(values: Iterable[float]) -> float | None:
    values = list(values)
    return statistics.fmean(values) if values else None


def stdev_or_none(values: Iterable[float]) -> float | None:
    values = list(values)
    return statistics.stdev(values) if len(values) >= 2 else None


def main() -> None:
    if not RUNS_DIR.is_dir():
        raise FileNotFoundError(
            f"Run directory not found: {RUNS_DIR}\n"
            "Generate the bulk jobs first with: python -m studies.main.bulk_lammps"
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    per_run_records: list[dict[str, object]] = []
    missing: list[str] = []
    invalid: list[str] = []

    run_dirs = sorted(path for path in RUNS_DIR.iterdir() if path.is_dir())
    for run_dir in run_dirs:
        parsed = parse_run_name(run_dir.name)
        if parsed is None:
            continue

        alloy, seed, temperature = parsed
        averages_path = run_dir / AVERAGES_FILENAME
        if not averages_path.is_file():
            missing.append(run_dir.name)
            continue

        try:
            rows = read_averages(averages_path)
        except (OSError, ValueError) as exc:
            invalid.append(f"{run_dir.name}: {exc}")
            continue

        if not rows:
            missing.append(f"{run_dir.name} (file has no data rows)")
            continue

        lattice_values = [row["a_fcc_A"] for row in rows]
        run_mean = statistics.fmean(lattice_values)
        run_sd = statistics.stdev(lattice_values) if len(lattice_values) >= 2 else None

        for row in rows:
            per_run_records.append({
                "alloy": alloy,
                "seed": seed,
                "temperature_K": temperature,
                "run": run_dir.name,
                "timestep": int(row["timestep"]),
                **{column: row[column] for column in VALUE_COLUMNS},
                "run_mean_a_fcc_A": run_mean,
                "run_sd_a_fcc_A": run_sd,
                "n_average_rows_in_run": len(rows),
                "complete_marker_present": (run_dir / "COMPLETE").exists(),
            })

    per_run_fields = [
        "alloy", "seed", "temperature_K", "run", "timestep",
        *VALUE_COLUMNS, "run_mean_a_fcc_A", "run_sd_a_fcc_A",
        "n_average_rows_in_run", "complete_marker_present",
    ]
    with RUN_OUTPUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=per_run_fields)
        writer.writeheader()
        writer.writerows(per_run_records)

    # Aggregate run-level means, so runs with more averaging rows do not get
    # unintentionally weighted more heavily than runs with fewer rows.
    groups: dict[tuple[str, int], list[dict[str, object]]] = {}
    for record in per_run_records:
        key = (str(record["alloy"]), int(record["temperature_K"]))
        groups.setdefault(key, [])
        # One record per run is sufficient for the group; keep only the first
        # timestep row from each run.
        if not any(existing["run"] == record["run"] for existing in groups[key]):
            groups[key].append(record)

    summary_fields = [
        "alloy", "temperature_K", "n_runs", "n_seeds",
        "mean_a_fcc_A", "sd_between_run_means_A",
        "min_run_mean_a_fcc_A", "max_run_mean_a_fcc_A",
        "mean_e_atom_eV_per_atom", "mean_temperature_K",
        "mean_pressure_bar", "all_runs_marked_complete",
    ]
    summary_records: list[dict[str, object]] = []

    for (alloy, temperature), records in sorted(groups.items()):
        lattice_means = [float(r["run_mean_a_fcc_A"]) for r in records]
        energy_means = [float(r["e_atom_eV_per_atom"]) for r in records]
        temp_means = [float(r["measured_temperature_K"]) for r in records]
        pressure_means = [float(r["pressure_bar"]) for r in records]
        seeds = {int(r["seed"]) for r in records}

        summary_records.append({
            "alloy": alloy,
            "temperature_K": temperature,
            "n_runs": len(records),
            "n_seeds": len(seeds),
            "mean_a_fcc_A": statistics.fmean(lattice_means),
            "sd_between_run_means_A": (
                statistics.stdev(lattice_means) if len(lattice_means) >= 2 else ""
            ),
            "min_run_mean_a_fcc_A": min(lattice_means),
            "max_run_mean_a_fcc_A": max(lattice_means),
            "mean_e_atom_eV_per_atom": statistics.fmean(energy_means),
            "mean_temperature_K": statistics.fmean(temp_means),
            "mean_pressure_bar": statistics.fmean(pressure_means),
            "all_runs_marked_complete": all(
                bool(r["complete_marker_present"]) for r in records
            ),
        })

    with SUMMARY_OUTPUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_records)

    print(f"Wrote per-run data: {RUN_OUTPUT} ({len(per_run_records)} averaged rows)")
    print(f"Wrote summary data: {SUMMARY_OUTPUT} ({len(summary_records)} alloy-temperature groups)")

    if missing:
        print(f"\nPending/missing averages: {len(missing)}")
        for name in missing:
            print(f"  PENDING  {name}")

    if invalid:
        print(f"\nInvalid averages: {len(invalid)}")
        for message in invalid:
            print(f"  INVALID  {message}")

    if not per_run_records:
        print("\nNo usable averages found yet. Run this script again as jobs finish.")


if __name__ == "__main__":
    main()
