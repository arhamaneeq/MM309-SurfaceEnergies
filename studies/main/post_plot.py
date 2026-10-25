#!/usr/bin/env python3
"""Generate numerical-results plots using nominal thermostat setpoints.

Run from anywhere:
    python studies/main/plot_results.py

Inputs:
  studies/main/results/bulk_samples.csv
  studies/main/results/slab_samples.csv
  studies/main/results/processed/bulk_run_summary.csv
  studies/main/results/processed/slab_run_summary.csv
  studies/main/results/processed/slab_condition_summary.csv

Outputs:
  studies/main/results/plots/*.png and *.pdf

Temperature-dependent figures group results by nominal thermostat setpoint
(300, 500, 700 K by default), not by the instantaneous/averaged temperature.
Actual recorded temperatures remain in the raw CSVs for diagnostics.
Ternary plots show measured alloy-composition points only; no interpolation.
"""
from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = PROJECT_ROOT / "studies" / "main" / "results"
PROCESSED_DIR = RESULTS_DIR / "processed"
PLOTS_DIR = RESULTS_DIR / "plots"

# Edit this list if the project uses different thermostat setpoints.
NOMINAL_TEMPERATURES_K = (300.0, 500.0, 700.0)
# A measured temperature farther than this from every setpoint is not binned.
TEMPERATURE_BIN_TOLERANCE_K = 50.0
ELEMENTS = ("Fe", "Ni", "Al")


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(
            f"Required input not found: {path}\n"
            "Run collect_results.py and process_results.py before plotting."
        )
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def num(value: Any) -> float | None:
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def safe_stem(value: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in str(value))


def save(fig: Any, stem: str) -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(PLOTS_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(PLOTS_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  {stem}.png / {stem}.pdf")


def nominal_temperature(value: Any) -> float | None:
    """Map a measured temperature to its nearest configured thermostat setpoint."""
    measured = num(value)
    if measured is None or not NOMINAL_TEMPERATURES_K:
        return None
    target = min(NOMINAL_TEMPERATURES_K, key=lambda t: abs(t - measured))
    return target if abs(target - measured) <= TEMPERATURE_BIN_TOLERANCE_K else None


def with_nominal_temperature(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    out = []
    for row in rows:
        target = nominal_temperature(row.get("temperature_K"))
        out.append({**row, "_nominal_temperature_K": target})
    return out


def grouped(rows: list[dict[str, Any]], keys: tuple[str, ...]):
    groups: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[tuple(str(row.get(k, "") or "") for k in keys)].append(row)
    return groups


def plot_run_spread(
    rows: list[dict[str, Any]], metric: str, ylabel: str, stem: str,
    title: str, group_keys: tuple[str, ...],
) -> None:
    """Plot per-run means against nominal temperature with between-run SD."""
    rows = with_nominal_temperature(rows)
    groups = grouped(rows, group_keys)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    made = False
    for key, group in sorted(groups.items()):
        by_temp: dict[float, list[float]] = defaultdict(list)
        for row in group:
            x, y = row.get("_nominal_temperature_K"), num(row.get(metric))
            if x is not None and y is not None:
                by_temp[float(x)].append(y)
        xs = sorted(by_temp)
        if not xs:
            continue
        ys, errors = [], []
        for x in xs:
            vals = by_temp[x]
            mu = sum(vals) / len(vals)
            ys.append(mu)
            errors.append(
                math.sqrt(sum((v - mu) ** 2 for v in vals) / (len(vals) - 1))
                if len(vals) >= 2 else 0.0
            )
        label = " | ".join(f"{k}={v}" for k, v in zip(group_keys, key) if v)
        ax.errorbar(xs, ys, yerr=errors, marker="o", capsize=3,
                    linewidth=1.5, label=label)
        made = True
    if not made:
        plt.close(fig)
        print(f"  SKIP {stem}: no usable run summaries")
        return
    ax.set_xticks(NOMINAL_TEMPERATURES_K)
    ax.set_xlabel("Nominal thermostat setpoint (K)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8, frameon=False)
    save(fig, stem)


def plot_surface_energy_temperature(slab_runs: list[dict[str, str]]) -> None:
    """One figure per alloy; orientation curves grouped by nominal setpoint."""
    rows = with_nominal_temperature(slab_runs)
    by_alloy = grouped(rows, ("alloy",))
    for (alloy,), alloy_rows in sorted(by_alloy.items()):
        by_plane = grouped(alloy_rows, ("miller_indices",))
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        made = False
        for (hkl,), plane_rows in sorted(by_plane.items()):
            by_temp: dict[float, list[float]] = defaultdict(list)
            for row in plane_rows:
                t = row.get("_nominal_temperature_K")
                g = num(row.get("gamma_J_per_m2_mean"))
                if t is not None and g is not None:
                    by_temp[float(t)].append(g)
            xs = sorted(by_temp)
            if not xs:
                continue
            ys, errors = [], []
            for t in xs:
                vals = by_temp[t]
                mu = sum(vals) / len(vals)
                ys.append(mu)
                errors.append(
                    math.sqrt(sum((v - mu) ** 2 for v in vals) / (len(vals) - 1))
                    if len(vals) >= 2 else 0.0
                )
            ax.errorbar(xs, ys, yerr=errors, marker="o", capsize=3,
                        linewidth=1.7, label=f"({hkl})")
            made = True
        if not made:
            plt.close(fig)
            continue
        ax.set_xticks(NOMINAL_TEMPERATURES_K)
        ax.set_xlabel("Nominal thermostat setpoint (K)")
        ax.set_ylabel("Surface energy (J m$^{-2}$)")
        ax.set_title(f"Surface energy vs temperature — {alloy}")
        ax.grid(True, alpha=0.25)
        ax.legend(title="Orientation", frameon=False)
        save(fig, f"slab_surface_energy_vs_temperature_{safe_stem(alloy)}")


def parse_composition(alloy: str) -> dict[str, float] | None:
    """Parse alloy labels; FeNiAl denotes the equiatomic ternary composition."""
    label = str(alloy).strip()
    # The project's equiatomic structures are labelled simply "FeNiAl",
    # rather than "Fe33.3Ni33.3Al33.3".
    if re.fullmatch(r"FeNiAl", label, flags=re.IGNORECASE):
        return {el: 100.0 / 3.0 for el in ELEMENTS}

    matches = re.findall(r"([A-Z][a-z]?)(\d+(?:\.\d+)?)", label)
    comp = {element: float(amount) for element, amount in matches}
    if set(comp) != set(ELEMENTS):
        return None
    total = sum(comp.values())
    if total <= 0:
        return None
    return {el: 100.0 * comp[el] / total for el in ELEMENTS}


def composition_coordinate_label(composition: dict[str, float]) -> str:
    """Return normalized (Fe, Ni, Al) atomic fractions for a LaTeX point label."""
    fractions = [composition[el] / 100.0 for el in ELEMENTS]
    return r"$(" + r",\;".join(f"{x:.3f}" for x in fractions) + r")$"


def ternary_xy(composition: dict[str, float]) -> tuple[float, float]:
    fe, ni, al = (composition[e] / 100.0 for e in ELEMENTS)
    return ni + 0.5 * al, (math.sqrt(3) / 2.0) * al


def draw_ternary_axes(ax: Any) -> None:
    height = math.sqrt(3) / 2.0
    vertices = [(0, 0), (1, 0), (0.5, height), (0, 0)]
    ax.plot([p[0] for p in vertices], [p[1] for p in vertices], color="black", lw=1.2)
    for pct in (20, 40, 60, 80):
        # Constant Fe, Ni, and Al composition gridlines.
        segments = [
            ({"Fe": pct, "Ni": 100-pct, "Al": 0},
             {"Fe": pct, "Ni": 0, "Al": 100-pct}),
            ({"Fe": 100-pct, "Ni": pct, "Al": 0},
             {"Fe": 0, "Ni": pct, "Al": 100-pct}),
            ({"Fe": 100-pct, "Ni": 0, "Al": pct},
             {"Fe": 0, "Ni": 100-pct, "Al": pct}),
        ]
        for a, b in segments:
            x1, y1 = ternary_xy(a)
            x2, y2 = ternary_xy(b)
            ax.plot([x1, x2], [y1, y2], color="0.82", lw=0.6, zorder=0)
    ax.text(-0.035, -0.035, "Fe", ha="right", va="top", fontsize=11)
    ax.text(1.035, -0.035, "Ni", ha="left", va="top", fontsize=11)
    ax.text(0.5, height + 0.045, "Al", ha="center", va="bottom", fontsize=11)
    ax.set_xlim(-0.10, 1.10)
    ax.set_ylim(-0.08, height + 0.12)
    ax.set_aspect("equal")
    ax.axis("off")


def plot_ternary_surface_energy(slab_conditions: list[dict[str, str]]) -> None:
    """One ternary scatter per nominal temperature and orientation; no interpolation."""
    valid = []
    for row in slab_conditions:
        comp = parse_composition(row.get("alloy", ""))
        gamma = num(row.get("gamma_J_per_m2_between_run_mean"))
        target = nominal_temperature(row.get("temperature_K"))
        if comp is None or gamma is None or target is None:
            continue
        valid.append({**row, "_composition": comp, "_gamma": gamma,
                      "_nominal_temperature_K": target})
    if not valid:
        print("  SKIP ternary plots: no parseable compositions/energies in temperature bins")
        return

    print(f"  Ternary points retained: {len(valid)} / {len(slab_conditions)} condition rows")
    print("  Ternary point counts by nominal temperature and orientation:")
    for key, group in sorted(grouped(valid, ("_nominal_temperature_K", "miller_indices")).items()):
        print(f"    T={float(key[0]):g} K, ({key[1]}): {len(group)} point(s)")

    for (temp, hkl), rows in sorted(grouped(
        valid, ("_nominal_temperature_K", "miller_indices")
    ).items()):
        fig, ax = plt.subplots(figsize=(6.8, 5.8))
        draw_ternary_axes(ax)
        xs, ys, values = [], [], []
        for row in rows:
            x, y = ternary_xy(row["_composition"])
            xs.append(x)
            ys.append(y)
            values.append(row["_gamma"])
        scatter = ax.scatter(xs, ys, c=values, cmap="viridis", s=145,
                             edgecolors="black", linewidths=0.7, zorder=3)
        for row, x, y, value in zip(rows, xs, ys, values):
            coord_label = composition_coordinate_label(row["_composition"])
            ax.annotate(
                coord_label,
                (x, y), xytext=(5, 5), textcoords="offset points",
                fontsize=8, ha="left", va="bottom",
            )
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.045, pad=0.03)
        cbar.set_label(r"Surface free energy, $\gamma$ (J m$^{-2}$)")
        ax.set_title(
            rf"$\gamma_{{{hkl}}}$ at $T = {float(temp):g}\,\mathrm{{K}}$"
        )
        save(fig, f"ternary_gamma_T{int(float(temp))}_hkl{safe_stem(hkl)}")


def plot_time_series(
    samples: list[dict[str, str]], metric: str, ylabel: str, prefix: str,
    max_runs: int = 100,
) -> None:
    by_run = grouped(samples, ("run",))
    count = 0
    for (run,), rows in sorted(by_run.items()):
        points = sorted(
            (x, y) for row in rows
            if (x := num(row.get("timestep"))) is not None
            and (y := num(row.get(metric))) is not None
        )
        if len(points) < 2:
            continue
        fig, ax = plt.subplots(figsize=(7.2, 3.8))
        ax.plot([p[0] for p in points], [p[1] for p in points], linewidth=1.0)
        ax.set_xlabel("Timestep")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{run}: {metric}")
        ax.grid(True, alpha=0.25)
        save(fig, f"{prefix}_{safe_stem(run)}")
        count += 1
        if count >= max_runs:
            print(f"  Reached time-series figure limit ({max_runs}) for {prefix}")
            break


def main() -> None:
    bulk_samples = read_csv(RESULTS_DIR / "bulk_samples.csv")
    slab_samples = read_csv(RESULTS_DIR / "slab_samples.csv")
    bulk_runs = read_csv(PROCESSED_DIR / "bulk_run_summary.csv")
    slab_runs = read_csv(PROCESSED_DIR / "slab_run_summary.csv")
    slab_conditions = read_csv(PROCESSED_DIR / "slab_condition_summary.csv")

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    print("Generating plots grouped by nominal thermostat setpoint.")
    print(f"Setpoints: {NOMINAL_TEMPERATURES_K}; bin tolerance: "
          f"+/-{TEMPERATURE_BIN_TOLERANCE_K:g} K")

    plot_run_spread(
        bulk_runs, "a_fcc_A_mean", "Lattice parameter (Å)",
        "bulk_lattice_parameter_vs_temperature",
        "Bulk FCC lattice parameter vs temperature", ("alloy",),
    )
    plot_run_spread(
        bulk_runs, "e_atom_eV_per_atom_mean", "Potential energy (eV atom$^{-1}$)",
        "bulk_energy_vs_temperature",
        "Bulk potential energy per atom vs temperature", ("alloy",),
    )
    plot_run_spread(
        bulk_runs, "pressure_bar_mean", "Pressure (bar)",
        "bulk_pressure_vs_temperature", "Bulk pressure vs temperature", ("alloy",),
    )
    plot_run_spread(
        bulk_runs, "volume_A3_mean", "Volume (Å$^3$)",
        "bulk_volume_vs_temperature", "Bulk volume vs temperature", ("alloy",),
    )
    plot_surface_energy_temperature(slab_runs)
    plot_ternary_surface_energy(slab_conditions)

    plot_time_series(bulk_samples, "a_fcc_A", "Lattice parameter (Å)",
                     "bulk_lattice_timeseries")
    plot_time_series(bulk_samples, "e_atom_eV_per_atom", "Energy (eV atom$^{-1}$)",
                     "bulk_energy_timeseries")
    plot_time_series(slab_samples, "gamma_J_per_m2", "Surface energy (J m$^{-2}$)",
                     "slab_gamma_timeseries")
    plot_time_series(slab_samples, "temperature_K", "Temperature (K)",
                     "slab_temperature_timeseries")

    print(f"Plots written to: {PLOTS_DIR.relative_to(PROJECT_ROOT)}")
    print("Each figure is saved as PNG (300 dpi) and PDF.")
    print("Ternary plots use nominal-temperature bins and measured composition points only.")
    print("No raw or processed CSVs are modified.")


if __name__ == "__main__":
    main()
