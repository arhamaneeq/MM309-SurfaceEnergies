#!/usr/bin/env python3
"""Numerically process collected bulk/slab LAMMPS samples; create CSVs, no plots.

Expected inputs (created by collect_results.py):
  studies/main/results/bulk_samples.csv
  studies/main/results/slab_samples.csv
  studies/main/results/run_inventory.csv

Run from anywhere:
  python studies/main/process_results.py

Outputs are written to studies/main/results/processed/. This script does not
modify raw inputs and does not claim that COMPLETE markers prove validity.
Sample-level statistics describe recorded time-series rows; run-level statistics
treat simulation runs/configurations as the replicate units where appropriate.
"""
from __future__ import annotations

import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = PROJECT_ROOT / "studies" / "main" / "results"
PROCESSED_DIR = RESULTS_DIR / "processed"
BULK_FILE = RESULTS_DIR / "bulk_samples.csv"
SLAB_FILE = RESULTS_DIR / "slab_samples.csv"
INVENTORY_FILE = RESULTS_DIR / "run_inventory.csv"

BULK_METRICS = [
    "a_fcc_A", "e_atom_eV_per_atom", "temperature_K", "pressure_bar", "volume_A3"
]
SLAB_METRICS = [
    "e_atom_eV_per_atom", "temperature_K", "gamma_eV_per_A2", "gamma_J_per_m2"
]
MIN_TREND_POINTS = 3


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(
            f"Required input not found: {path}\n"
            "Run collect_results.py first, and check PROJECT_ROOT / file placement."
        )
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def number(value: Any) -> float | None:
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def as_int(value: Any) -> int | None:
    val = number(value)
    return int(val) if val is not None and val.is_integer() else None


def mean_or_blank(values: list[float]) -> float | str:
    return statistics.fmean(values) if values else ""


def std_or_blank(values: list[float]) -> float | str:
    return statistics.stdev(values) if len(values) >= 2 else ""


def sem_or_blank(values: list[float]) -> float | str:
    return statistics.stdev(values) / math.sqrt(len(values)) if len(values) >= 2 else ""


def ci95_halfwidth(values: list[float]) -> float | str:
    # Normal-approximation interval; clearly labelled as such in output.
    if len(values) < 2:
        return ""
    return 1.96 * statistics.stdev(values) / math.sqrt(len(values))


def cv_percent(values: list[float]) -> float | str:
    if len(values) < 2:
        return ""
    m = statistics.fmean(values)
    if abs(m) < 1e-15:
        return ""
    return 100.0 * statistics.stdev(values) / abs(m)


def safe_range(values: list[float]) -> float | str:
    return max(values) - min(values) if values else ""


def lag1_autocorrelation(values: list[float]) -> float | str:
    if len(values) < 3:
        return ""
    m = statistics.fmean(values)
    denom = sum((x - m) ** 2 for x in values)
    if denom <= 0:
        return ""
    return sum((values[i] - m) * (values[i - 1] - m)
               for i in range(1, len(values))) / denom


def effective_sample_size(values: list[float]) -> float | str:
    """Crude AR(1)-based ESS diagnostic, not a general MCMC estimator."""
    rho = lag1_autocorrelation(values)
    if rho == "":
        return ""
    # Clamp to avoid nonsensical ESS from noisy estimates of negative/large rho.
    r = max(-0.99, min(0.99, float(rho)))
    return len(values) * (1.0 - r) / (1.0 + r)


def half_difference(values: list[float]) -> float | str:
    """Second-half mean minus first-half mean; diagnostic, not a convergence proof."""
    if len(values) < 4:
        return ""
    cut = len(values) // 2
    return statistics.fmean(values[cut:]) - statistics.fmean(values[:cut])


def linear_trend(xs: list[float], ys: list[float]) -> dict[str, float | str]:
    """Ordinary least-squares slope/intercept/R²; blank if fewer than 3 x points."""
    pairs = [(x, y) for x, y in zip(xs, ys)
             if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < MIN_TREND_POINTS or len({x for x, _ in pairs}) < MIN_TREND_POINTS:
        return {"slope_per_K": "", "intercept": "", "r_squared": ""}
    xvals = [x for x, _ in pairs]
    yvals = [y for _, y in pairs]
    xm, ym = statistics.fmean(xvals), statistics.fmean(yvals)
    sxx = sum((x - xm) ** 2 for x in xvals)
    if sxx <= 0:
        return {"slope_per_K": "", "intercept": "", "r_squared": ""}
    slope = sum((x - xm) * (y - ym) for x, y in pairs) / sxx
    intercept = ym - slope * xm
    sst = sum((y - ym) ** 2 for y in yvals)
    sse = sum((y - (intercept + slope * x)) ** 2 for x, y in pairs)
    r2 = 1.0 - sse / sst if sst > 0 else ""
    return {"slope_per_K": slope, "intercept": intercept, "r_squared": r2}


def grouped(rows: list[dict[str, Any]], keys: tuple[str, ...]) -> dict[tuple[str, ...], list[dict[str, Any]]]:
    out: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        out[tuple(str(row.get(k, "") or "") for k in keys)].append(row)
    return out


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> int:
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row.keys()))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def load_structure_metadata(run_directory: str) -> dict[str, Any]:
    path = PROJECT_ROOT / run_directory / "structure.json"
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def sample_values(rows: list[dict[str, str]], metric: str) -> list[float]:
    return [v for row in rows if (v := number(row.get(metric))) is not None]


def summarize_metric(values: list[float], prefix: str) -> dict[str, Any]:
    return {
        f"{prefix}_n": len(values),
        f"{prefix}_mean": mean_or_blank(values),
        f"{prefix}_std": std_or_blank(values),
        f"{prefix}_sem": sem_or_blank(values),
        f"{prefix}_ci95_normal_halfwidth": ci95_halfwidth(values),
        f"{prefix}_min": min(values) if values else "",
        f"{prefix}_max": max(values) if values else "",
        f"{prefix}_range": safe_range(values),
        f"{prefix}_cv_percent": cv_percent(values),
        f"{prefix}_lag1_autocorrelation": lag1_autocorrelation(values),
        f"{prefix}_ess_ar1_approx": effective_sample_size(values),
        f"{prefix}_second_minus_first_half_mean": half_difference(values),
    }


def process_kind(
    samples: list[dict[str, str]],
    inventory: list[dict[str, str]],
    kind: str,
    metrics: list[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Return per-run summaries, per-condition summaries, and data-quality checks."""
    sample_by_run: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in samples:
        sample_by_run[row.get("run", "")].append(row)
    inv_by_run = {row.get("run", ""): row for row in inventory if row.get("kind") == kind}

    run_summaries: list[dict[str, Any]] = []
    quality: list[dict[str, Any]] = []
    for run in sorted(set(sample_by_run) | {r for r, v in inv_by_run.items()}):
        rows = sample_by_run.get(run, [])
        inv = inv_by_run.get(run, {})
        meta = load_structure_metadata(inv.get("run_directory", ""))
        base: dict[str, Any] = {
            "kind": kind, "run": run,
            "alloy": (rows[0].get("alloy", "") if rows else inv.get("alloy", "")),
            "temperature_K": (rows[0].get("temperature_K", "") if rows else inv.get("temperature_K", "")),
            "seed": (rows[0].get("seed", rows[0].get("parent_seed", "")) if rows else inv.get("seed", "")),
            "configuration": (rows[0].get("configuration", "") if rows else inv.get("configuration", "")),
            "miller_indices": (rows[0].get("miller_indices", "") if rows else inv.get("miller_indices", "")),
            "complete_marker_present": inv.get("complete_marker_present", ""),
            "inventory_status": inv.get("status", "samples_not_in_inventory"),
            "sample_rows": len(rows),
        }
        report = meta.get("report", {}) if isinstance(meta.get("report", {}), dict) else {}
        base["surface_area_A2"] = report.get("surface_area_angstrom2", "")
        base["n_atoms_from_structure_metadata"] = report.get("n_atoms", "")
        for metric in metrics:
            values = sample_values(rows, metric)
            base.update(summarize_metric(values, metric))
        run_summaries.append(base)

        timesteps = [number(r.get("timestep")) for r in rows]
        valid_steps = [x for x in timesteps if x is not None]
        diffs = [valid_steps[i] - valid_steps[i - 1] for i in range(1, len(valid_steps))]
        temp_values = sample_values(rows, "temperature_K")
        nonfinite_by_metric = {
            m: sum(1 for r in rows if r.get(m, "").strip() and number(r.get(m)) is None)
            for m in metrics
        }
        flags: list[str] = []
        if not rows:
            flags.append("no_numeric_samples")
        if inv.get("complete_marker_present", "").lower() not in ("true", "1", "yes"):
            flags.append("missing_COMPLETE_marker")
        if inv.get("averages_file_present", "").lower() not in ("true", "1", "yes"):
            flags.append("missing_averages_file")
        if len(valid_steps) != len(rows):
            flags.append("invalid_or_missing_timestep")
        if any(d <= 0 for d in diffs):
            flags.append("non_increasing_or_duplicate_timestep")
        if temp_values and number(base["temperature_K"]) is not None:
            target = number(base["temperature_K"])
            if target is not None and abs(statistics.fmean(temp_values) - target) > max(10.0, 0.05 * abs(target)):
                flags.append("mean_temperature_far_from_target")
        if any(v for v in nonfinite_by_metric.values()):
            flags.append("invalid_numeric_values")
        quality.append({
            "kind": kind, "run": run, "inventory_status": inv.get("status", ""),
            "complete_marker_present": inv.get("complete_marker_present", ""),
            "averages_file_present": inv.get("averages_file_present", ""),
            "sample_rows": len(rows), "valid_timestep_rows": len(valid_steps),
            "timestep_min": min(valid_steps) if valid_steps else "",
            "timestep_max": max(valid_steps) if valid_steps else "",
            "timestep_nonincreasing_count": sum(1 for d in diffs if d <= 0),
            "mean_sample_temperature_K": mean_or_blank(temp_values),
            "target_temperature_K": base["temperature_K"],
            "nonfinite_or_unparseable_values_by_metric": json.dumps(nonfinite_by_metric, sort_keys=True),
            "quality_flags": ";".join(flags) if flags else "no_basic_flags",
            "note": "Automated checks only; inspect input decks/logs and physical convergence separately.",
        })

    # Condition summaries retain both pooled time-series rows and run-level replicate means.
    condition_keys = ("alloy", "temperature_K") if kind == "bulk" else (
        "alloy", "temperature_K", "miller_indices"
    )
    condition_groups = grouped(run_summaries, condition_keys)
    condition_summaries: list[dict[str, Any]] = []
    for key, runs in sorted(condition_groups.items()):
        base = dict(zip(condition_keys, key))
        base["n_runs"] = len(runs)
        base["runs_with_COMPLETE_marker"] = sum(
            str(r.get("complete_marker_present", "")).lower() in ("true", "1", "yes") for r in runs
        )
        base["runs_with_samples"] = sum((as_int(r.get("sample_rows")) or 0) > 0 for r in runs)
        for metric in metrics:
            # Run-level means: the replicate-level summary, avoiding false precision from
            # treating correlated fix ave/time rows as independent experiments.
            run_means = [v for r in runs if (v := number(r.get(f"{metric}_mean"))) is not None]
            base.update({
                f"{metric}_between_run_mean": mean_or_blank(run_means),
                f"{metric}_between_run_std": std_or_blank(run_means),
                f"{metric}_between_run_sem": sem_or_blank(run_means),
                f"{metric}_between_run_ci95_normal_halfwidth": ci95_halfwidth(run_means),
                f"{metric}_between_run_min": min(run_means) if run_means else "",
                f"{metric}_between_run_max": max(run_means) if run_means else "",
            })
        # Add pooled sample-row metrics for comparison, explicitly named as pooled.
        matching_rows = [
            row for row in samples
            if all(str(row.get(k, "") or "") == key[i] for i, k in enumerate(condition_keys))
        ]
        for metric in metrics:
            vals = sample_values(matching_rows, metric)
            base.update({
                f"{metric}_pooled_sample_n": len(vals),
                f"{metric}_pooled_sample_mean": mean_or_blank(vals),
                f"{metric}_pooled_sample_std": std_or_blank(vals),
                f"{metric}_pooled_sample_min": min(vals) if vals else "",
                f"{metric}_pooled_sample_max": max(vals) if vals else "",
            })
        condition_summaries.append(base)
    return run_summaries, condition_summaries, quality


def add_temperature_trends(condition_rows: list[dict[str, Any]], kind: str,
                           metrics: list[str]) -> list[dict[str, Any]]:
    """Fit descriptive linear trends to condition-level between-run means."""
    keys = ("alloy",) if kind == "bulk" else ("alloy", "miller_indices")
    buckets = grouped(condition_rows, keys)
    out: list[dict[str, Any]] = []
    for key, rows in sorted(buckets.items()):
        for metric in metrics:
            pairs = [(number(r.get("temperature_K")), number(r.get(f"{metric}_between_run_mean")))
                     for r in rows]
            pairs = [(x, y) for x, y in pairs if x is not None and y is not None]
            trend = linear_trend([x for x, _ in pairs], [y for _, y in pairs])
            out.append({
                **dict(zip(keys, key)), "metric": metric,
                "n_temperature_conditions": len(pairs),
                "temperature_min_K": min((x for x, _ in pairs), default=""),
                "temperature_max_K": max((x for x, _ in pairs), default=""),
                **trend,
                "interpretation": "descriptive linear fit only; not a phase-transition model",
            })
    return out


def make_orientation_contrasts(slab_conditions: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Differences in condition means relative to (100), when available."""
    by_alloy_temp: dict[tuple[str, str], dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in slab_conditions:
        by_alloy_temp[(str(row.get("alloy", "")), str(row.get("temperature_K", "")))][
            str(row.get("miller_indices", ""))
        ] = row
    out: list[dict[str, Any]] = []
    for (alloy, temp), planes in sorted(by_alloy_temp.items()):
        ref = planes.get("100")
        if not ref:
            continue
        for hkl, row in sorted(planes.items()):
            if hkl == "100":
                continue
            a = number(row.get("gamma_J_per_m2_between_run_mean"))
            b = number(ref.get("gamma_J_per_m2_between_run_mean"))
            out.append({
                "alloy": alloy, "temperature_K": temp,
                "comparison": f"({hkl}) minus (100)", "reference_plane": "100",
                "gamma_mean_hkl_J_per_m2": "" if a is None else a,
                "gamma_mean_100_J_per_m2": "" if b is None else b,
                "delta_gamma_J_per_m2": "" if a is None or b is None else a - b,
                "note": "Difference of condition-level run means; not a paired test.",
            })
    return out


def main() -> None:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    bulk_samples = read_csv(BULK_FILE)
    slab_samples = read_csv(SLAB_FILE)
    inventory = read_csv(INVENTORY_FILE)

    bulk_inv = [r for r in inventory if r.get("kind") == "bulk"]
    slab_inv = [r for r in inventory if r.get("kind") == "slab"]
    bulk_runs, bulk_conditions, bulk_quality = process_kind(
        bulk_samples, bulk_inv, "bulk", BULK_METRICS
    )
    slab_runs, slab_conditions, slab_quality = process_kind(
        slab_samples, slab_inv, "slab", SLAB_METRICS
    )

    bulk_trends = add_temperature_trends(bulk_conditions, "bulk", BULK_METRICS)
    slab_trends = add_temperature_trends(slab_conditions, "slab", SLAB_METRICS)
    orientation = make_orientation_contrasts(slab_conditions)

    outputs = [
        ("bulk_run_summary.csv", bulk_runs),
        ("bulk_condition_summary.csv", bulk_conditions),
        ("slab_run_summary.csv", slab_runs),
        ("slab_condition_summary.csv", slab_conditions),
        ("bulk_temperature_trends.csv", bulk_trends),
        ("slab_temperature_trends.csv", slab_trends),
        ("slab_orientation_contrasts.csv", orientation),
        ("data_quality_checks.csv", [*bulk_quality, *slab_quality]),
    ]
    print("Numerical processing complete; no plots created.")
    for filename, rows in outputs:
        count = write_csv(PROCESSED_DIR / filename, rows)
        print(f"  {filename}: {count} rows")
    print(f"Output directory: {PROCESSED_DIR.relative_to(PROJECT_ROOT)}")
    print()
    print("Interpretation cautions:")
    print("- CI half-widths use a normal approximation and are descriptive, especially for small n.")
    print("- Pooled sample rows are time-series observations, not independent replicates.")
    print("- Between-run summaries use run means as replicate units; verify run independence.")
    print("- AR(1) ESS is a rough diagnostic, not a rigorous uncertainty correction.")
    print("- Linear temperature fits are descriptive and do not establish phase stability.")
    print("- Gamma values are processed as recorded in slab_averages.dat; this script does not")
    print("  recompute gamma because the collected CSV lacks a guaranteed atom count/reference mapping.")
    print("- Basic quality flags cannot establish equilibration, FCC stability, or physical validity.")


if __name__ == "__main__":
    main()
