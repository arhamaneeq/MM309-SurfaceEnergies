#!/usr/bin/env python3
"""Collect raw bulk and slab LAMMPS averages into tidy CSV files.

Run from anywhere:
    python studies/main/collect_results.py

Writes studies/main/results/bulk_samples.csv, slab_samples.csv, and
run_inventory.csv. This script only collects and labels existing data; it
does not calculate statistics or interpret results.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any, Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[2]
MAIN_DIR = PROJECT_ROOT / "studies" / "main"
BULK_RUNS_DIR = MAIN_DIR / "runs" / "bulk_equilibration"
SLAB_RUNS_DIR = MAIN_DIR / "runs" / "slab_calculations"
RESULTS_DIR = MAIN_DIR / "results"

BULK_COLUMNS = [
    "timestep", "a_fcc_A", "e_atom_eV_per_atom", "temperature_K",
    "pressure_bar", "volume_A3",
]
SLAB_COLUMNS = [
    "timestep", "e_atom_eV_per_atom", "temperature_K",
    "gamma_eV_per_A2", "gamma_J_per_m2",
]


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def parse_run_name(name: str, kind: str) -> dict[str, Any]:
    result: dict[str, Any] = {"run": name}
    if kind == "bulk":
        match = re.match(r"(?P<alloy>.+)_seed(?P<seed>\d+)_T(?P<temp>\d+)$", name)
        if match:
            result.update(alloy=match["alloy"], seed=int(match["seed"]),
                          temperature_K=int(match["temp"]))
    else:
        match = re.match(
            r"(?P<alloy>.+)_T(?P<temp>\d+)_(?P<config>.+)_(?P<hkl>\d{3})$",
            name,
        )
        if match:
            result.update(alloy=match["alloy"], temperature_K=int(match["temp"]),
                          configuration=match["config"], miller_indices=match["hkl"])
    return result


def parse_numeric_rows(path: Path, expected_values: int) -> list[list[float]]:
    """Read numeric LAMMPS fix ave/time rows, ignoring comments and headers."""
    rows: list[list[float]] = []
    if not path.is_file():
        return rows
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            try:
                values = [float(field) for field in fields]
            except ValueError:
                continue
            if len(values) == expected_values:
                rows.append(values)
    return rows


def write_csv(path: Path, fieldnames: list[str],
              rows: Iterable[dict[str, Any]]) -> int:
    count = 0
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def collect_bulk() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    samples: list[dict[str, Any]] = []
    inventory: list[dict[str, Any]] = []
    if not BULK_RUNS_DIR.is_dir():
        return samples, inventory
    for run_dir in sorted(p for p in BULK_RUNS_DIR.iterdir() if p.is_dir()):
        meta = read_json(run_dir / "structure.json")
        ids = parse_run_name(run_dir.name, "bulk")
        alloy = meta.get("alloy", ids.get("alloy", ""))
        seed = meta.get("seed", ids.get("seed", ""))
        temp = ids.get("temperature_K", "")
        averages = run_dir / "bulk_averages.dat"
        complete = (run_dir / "COMPLETE").is_file()
        data = parse_numeric_rows(averages, len(BULK_COLUMNS))
        inventory.append({
            "kind": "bulk", "run": run_dir.name, "alloy": alloy,
            "temperature_K": temp, "seed": seed, "configuration": "",
            "miller_indices": "", "complete_marker_present": complete,
            "averages_file_present": averages.is_file(),
            "numeric_sample_rows": len(data),
            "status": ("complete_with_data" if complete and data else
                       "complete_marker_but_no_data" if complete else
                       "data_without_complete_marker" if data else "no_data"),
            "run_directory": str(run_dir.relative_to(PROJECT_ROOT)),
        })
        for index, values in enumerate(data, start=1):
            row = {"run": run_dir.name, "alloy": alloy,
                   "temperature_K": temp, "seed": seed,
                   "sample_index": index, "complete_marker_present": complete}
            row.update(dict(zip(BULK_COLUMNS, values)))
            samples.append(row)
    return samples, inventory


def collect_slabs() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    samples: list[dict[str, Any]] = []
    inventory: list[dict[str, Any]] = []
    if not SLAB_RUNS_DIR.is_dir():
        return samples, inventory
    for run_dir in sorted(p for p in SLAB_RUNS_DIR.iterdir() if p.is_dir()):
        meta = read_json(run_dir / "structure.json")
        ids = parse_run_name(run_dir.name, "slab")
        report = meta.get("report", {})
        if not isinstance(report, dict):
            report = {}
        alloy = meta.get("alloy", ids.get("alloy", ""))
        temp = meta.get("temperature_K", ids.get("temperature_K", ""))
        config = meta.get("configuration", ids.get("configuration", ""))
        hkl = meta.get("miller_indices", ids.get("miller_indices", ""))
        if isinstance(hkl, (list, tuple)):
            hkl = "".join(str(x) for x in hkl)
        area = report.get("surface_area_angstrom2", "")
        seed = meta.get("parent_random_seed", "")
        averages = run_dir / "slab_averages.dat"
        complete = (run_dir / "COMPLETE").is_file()
        data = parse_numeric_rows(averages, len(SLAB_COLUMNS))
        inventory.append({
            "kind": "slab", "run": run_dir.name, "alloy": alloy,
            "temperature_K": temp, "seed": seed, "configuration": config,
            "miller_indices": hkl, "complete_marker_present": complete,
            "averages_file_present": averages.is_file(),
            "numeric_sample_rows": len(data),
            "status": ("complete_with_data" if complete and data else
                       "complete_marker_but_no_data" if complete else
                       "data_without_complete_marker" if data else "no_data"),
            "run_directory": str(run_dir.relative_to(PROJECT_ROOT)),
        })
        for index, values in enumerate(data, start=1):
            row = {"run": run_dir.name, "alloy": alloy,
                   "temperature_K": temp, "configuration": config,
                   "parent_seed": seed, "miller_indices": hkl,
                   "surface_area_A2": area, "sample_index": index,
                   "complete_marker_present": complete}
            row.update(dict(zip(SLAB_COLUMNS, values)))
            samples.append(row)
    return samples, inventory


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    bulk_samples, bulk_inventory = collect_bulk()
    slab_samples, slab_inventory = collect_slabs()
    bulk_fields = [
        "run", "alloy", "temperature_K", "seed", "sample_index",
        "complete_marker_present", *BULK_COLUMNS,
    ]
    slab_fields = [
        "run", "alloy", "temperature_K", "configuration", "parent_seed",
        "miller_indices", "surface_area_A2", "sample_index",
        "complete_marker_present", *SLAB_COLUMNS,
    ]
    inventory_fields = [
        "kind", "run", "alloy", "temperature_K", "seed", "configuration",
        "miller_indices", "complete_marker_present", "averages_file_present",
        "numeric_sample_rows", "status", "run_directory",
    ]
    n_bulk = write_csv(RESULTS_DIR / "bulk_samples.csv", bulk_fields, bulk_samples)
    n_slab = write_csv(RESULTS_DIR / "slab_samples.csv", slab_fields, slab_samples)
    n_runs = write_csv(RESULTS_DIR / "run_inventory.csv", inventory_fields,
                       [*bulk_inventory, *slab_inventory])
    print("Collection complete; no analysis performed.")
    print(f"  Bulk samples: {n_bulk} rows -> studies/main/results/bulk_samples.csv")
    print(f"  Slab samples: {n_slab} rows -> studies/main/results/slab_samples.csv")
    print(f"  Run inventory: {n_runs} runs -> studies/main/results/run_inventory.csv")
    print("Check run_inventory.csv for missing or empty averages files.")


if __name__ == "__main__":
    main()

