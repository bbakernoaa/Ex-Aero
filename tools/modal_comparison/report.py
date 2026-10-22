#!/usr/bin/env python3
"""Program: modal_comparison_report. Author: EX-aero contributors.

Abstract: Run paired dry boxes at dt=1 and 0.5 s; publish CSV and a report.
History: 2026-10 initial MAM4xx comparison harness.
Usage: report.py --executable BUILD/tests/exaero_modal_mam4_comparison
Input: executable, committed reaction YAML, optional generated MKPP manifest.
Output: $DATA/modal.comparison.{csv,half.csv,summary.csv,md,metadata.json}.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)


def read_results(path: Path) -> pd.DataFrame:
    """Read and validate a complete comparison table.

    Parameters
    ----------
    path : Path
        CSV produced by the comparison executable.

    Returns
    -------
    pandas.DataFrame
        Unique experiment/case/time/field records with finite values.

    Examples
    --------
    >>> frame = read_results(Path(os.environ["DATA"]) / "modal.comparison.csv")
    """
    frame = pd.read_csv(path)
    keys = ["experiment", "case", "time_s", "field"]
    if (
        frame.duplicated(keys).any()
        or not np.isfinite(frame[["generic", "mam4xx"]]).all().all()
    ):
        raise ValueError("duplicate or nonfinite comparison results")
    if (frame[["generic", "mam4xx"]] < 0).any().any():
        raise ValueError("negative comparison state")
    if len(frame) != 9 * 3 * 5 * 37 or set(frame.time_s) != {0, 1, 60, 300, 600}:
        raise ValueError("incomplete experiment/case/time/field coverage")
    initial = frame.loc[frame.time_s == 0]
    np.testing.assert_array_equal(initial.generic, initial.mam4xx)
    controls = frame.experiment.isin(
        ["coagulation_control", "nucleation_linear", "nucleation_quadratic"]
    )
    np.testing.assert_allclose(
        frame.loc[controls, "generic"],
        frame.loc[controls, "mam4xx"],
        rtol=1e-10,
        atol=1e-25,
    )
    return frame


def summarize(coarse: pd.DataFrame, fine: pd.DataFrame) -> pd.DataFrame:
    """Compare final states and their sensitivity to a halved time step.

    Parameters
    ----------
    coarse, fine : pandas.DataFrame
        Results at dt=1 s and dt=0.5 s, respectively.

    Returns
    -------
    pandas.DataFrame
        One row per experiment, case and field. Symmetric relative differences
        divide by max(abs(generic), abs(mam4xx), 1e-30), without hiding zeros.

    Examples
    --------
    >>> summary = summarize(coarse_results, half_step_results)
    """
    keys = ["experiment", "case", "field", "unit"]
    coarse_end = coarse.loc[coarse.time_s == 600, keys + ["generic", "mam4xx"]]
    fine_end = fine.loc[fine.time_s == 600, keys + ["generic", "mam4xx"]]
    joined = coarse_end.merge(
        fine_end, on=keys, suffixes=("_1s", "_half"), validate="one_to_one"
    )
    joined["absolute_difference"] = abs(joined.generic_half - joined.mam4xx_half)
    joined["symmetric_relative_difference"] = joined.absolute_difference / np.maximum(
        np.maximum(abs(joined.generic_half), abs(joined.mam4xx_half)), 1e-30
    )
    for model in ("generic", "mam4xx"):
        joined[f"{model}_step_sensitivity"] = abs(
            joined[f"{model}_1s"] - joined[f"{model}_half"]
        ) / np.maximum(
            np.maximum(abs(joined[f"{model}_1s"]), abs(joined[f"{model}_half"])), 1e-30
        )
    return joined


def main() -> None:
    """Run the bounded benchmark and write provenance and reviewable results.

    Parameters
    ----------
    None
        Arguments and DATA come from the driver environment.

    Returns
    -------
    None
        Writes artifacts beneath DATA; raises on a failed executable or check.

    Examples
    --------
    >>> # DATA=/scratch/run python report.py --executable build/tests/...
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--mkpp-manifest", type=Path)
    parser.add_argument(
        "--mkpp-source", type=Path, help="Compiler checkout for revision provenance"
    )
    args = parser.parse_args()
    if "DATA" not in os.environ:
        parser.error("FATAL ERROR: DATA must identify the output directory")
    destination = Path(os.environ["DATA"]).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    executable = args.executable.resolve(strict=True)
    repo = Path(__file__).resolve().parents[2]
    env = {**os.environ, "OMP_NUM_THREADS": "2", "OMP_PROC_BIND": "false"}
    for suffix, flags in (("csv", []), ("half.csv", ["--half-step"])):
        with (destination / f"modal.comparison.{suffix}").open("w") as output:
            subprocess.run(
                [str(executable), *flags],
                stdout=output,
                check=True,
                env=env,
                timeout=600,
            )
    coarse = read_results(destination / "modal.comparison.csv")
    fine = read_results(destination / "modal.comparison.half.csv")
    summary = summarize(coarse, fine)
    summary.to_csv(destination / "modal.comparison.summary.csv", index=False)
    keys = ["experiment", "case", "field", "unit"]
    first = coarse.loc[coarse.time_s == 1].merge(
        coarse.loc[coarse.time_s == 0, keys + ["generic", "mam4xx"]],
        on=keys,
        suffixes=("_1s", "_initial"),
        validate="one_to_one",
    )
    first["generic_tendency_per_s"] = first.generic_1s - first.generic_initial
    first["mam4xx_tendency_per_s"] = first.mam4xx_1s - first.mam4xx_initial
    first["unit"] = first.unit + "/s"
    first[keys + ["generic_tendency_per_s", "mam4xx_tendency_per_s"]].to_csv(
        destination / "modal.comparison.first_step.csv", index=False
    )
    metadata = {
        "exaero_commit": subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
        ).strip(),
        "exaero_dirty": bool(
            subprocess.check_output(
                ["git", "-C", str(repo), "status", "--porcelain"], text=True
            ).strip()
        ),
        "mam4xx_commit": subprocess.check_output(
            ["git", "-C", str(repo / "extern/mam4xx"), "rev-parse", "HEAD"], text=True
        ).strip(),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "reaction_yaml_sha256": hashlib.sha256(
            (repo / "tools/modal_comparison/reactions.yaml").read_bytes()
        ).hexdigest(),
        "chemistry_rhs": sorted(coarse.chemistry_rhs.unique().tolist()),
        "dt_s": [1, 0.5],
        "duration_s": 600,
        "cases": [
            {
                "temperature_K": 250,
                "pressure_Pa": 60000,
                "acid_molecules_cm3": 1e5,
                "size_scale": 0.7,
            },
            {
                "temperature_K": 280,
                "pressure_Pa": 90000,
                "acid_molecules_cm3": 1e6,
                "size_scale": 1,
            },
            {
                "temperature_K": 300,
                "pressure_Pa": 101325,
                "acid_molecules_cm3": 1e7,
                "size_scale": 1.5,
            },
        ],
        "rh": 0,
        "source_sha256": {
            name: hashlib.sha256((repo / name).read_bytes()).hexdigest()
            for name in (
                "tests/modal_mam4_comparison.cpp",
                "tests/modal_mam4_reference.cpp",
                "tests/modal_mam4_reference.hpp",
                "tools/modal_comparison/reactions.yaml",
                "tools/modal_comparison/report.py",
            )
        },
    }
    if args.mkpp_manifest:
        metadata["mkpp_manifest"] = json.loads(args.mkpp_manifest.read_text())
        generated = args.mkpp_manifest.parent / "modal_mam4_comparison"
        metadata["generated_rates_sha256"] = hashlib.sha256(
            (generated / "rates.hpp").read_bytes()
        ).hexdigest()
    if args.mkpp_source:
        metadata["mkpp_commit"] = subprocess.check_output(
            ["git", "-C", str(args.mkpp_source), "rev-parse", "HEAD"], text=True
        ).strip()
        metadata["mkpp_dirty"] = bool(
            subprocess.check_output(
                ["git", "-C", str(args.mkpp_source), "status", "--porcelain"], text=True
            ).strip()
        )
    if "mkpp_generated" in metadata["chemistry_rhs"] and not args.mkpp_manifest:
        parser.error(
            "FATAL ERROR: pass --mkpp-manifest to record generated-kernel provenance"
        )
    (destination / "modal.comparison.metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    selected = summary.loc[
        summary.field.isin(
            [
                "H2SO4_SO4_equivalent",
                "num_aitken_i",
                "num_primary_carbon_i",
                "BC_primary_carbon_i",
            ]
        )
    ]
    columns = [
        "experiment",
        "case",
        "field",
        "symmetric_relative_difference",
        "generic_step_sensitivity",
        "mam4xx_step_sensitivity",
    ]
    lines = [
        "# Generic modal / MAM4xx dry-box comparison",
        "",
        f"Chemistry RHS: {', '.join(metadata['chemistry_rhs'])}.",
        "All cases passed initial-state equality, finite/nonnegative state and per-species material budgets.",
        "Matched Wang linear/quadratic formation and native-coagulation bridge controls passed parity checks.",
        "",
        "Final states at 600 s use dt=0.5 s. Relative differences use max(abs(generic), abs(native), 1e-30).",
        "Step sensitivity compares dt=1 s and 0.5 s within each path; it is a diagnostic, not a convergence proof.",
        "Differences are measured results, not acceptance tolerances or a MAM4-equivalence claim.",
        "",
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    lines.extend(
        "| "
        + " | ".join(
            f"{value:.6g}" if isinstance(value, (float, np.floating)) else str(value)
            for value in row
        )
        + " |"
        for row in selected[columns].itertuples(index=False, name=None)
    )
    lines.extend(
        [
            "",
            "## Scope",
            "",
            "Native comparison calls pinned MAM4xx coagulation, dry H2SO4 uptake, Wang PBL cluster formation, and finite-step primary-carbon aging kernels.",
            "Both paths use a common explicit split driver and a controlled 100-molecule seed; this does not execute native amicphys or native cluster survival/growth to Aitken.",
            "Aging-only has a budgeted 1e-15 kg/m3/s sulfate coating source; native and generic aging formulations differ.",
            "MKPP, when enabled, executes the frozen stage RHS and checks it against independent host rates; this driver does not exercise MKPP's implicit solver or fully symbolic geometry feedback.",
            "Wet condensation, thermodynamic organic partitioning, growth renaming, cloud processes, deposition, optics and GPU parity remain outside these dry-box cases.",
            "",
            "See modal.comparison.summary.csv for every mass, number, gas and dry-radius field.",
        ]
    )
    (destination / "modal.comparison.md").write_text("\n".join(lines) + "\n")
    print(f"Comparison report: {destination / 'modal.comparison.md'}")


if __name__ == "__main__":
    main()
