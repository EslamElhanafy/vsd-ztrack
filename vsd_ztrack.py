#!/usr/bin/env python3
"""
VSD-ZTrack: residue-level z-position analysis for voltage-sensor domains.

This script calculates the position of selected residues along a Cartesian axis
(default: z) relative to a hydrophobic constriction site (HCS) reference residue.
It was designed for voltage-gated ion channels, but it can be used for any
protein system where a reference residue and target residues are known.

Default metric
--------------
For each selected target residue and each analyzed trajectory frame:

    relative_position_A = axis(COM of target residue) - axis(CA of HCS residue)

By default:
    * HCS reference = C-alpha atom of the HCS residue.
    * Target residue position = side-chain heavy-atom center of mass.
    * Axis = z.

Two use modes
-------------
1. Residue-map mode
   The script reads a CSV table containing family/gene/accession/VSD/HCS/R1-R6
   residue annotations. This is the recommended mode for Nav, Cav, Kv, HCN, CNG,
   Hv, and TRPP channels included in the distributed residue map.

2. Manual-target mode
   The user supplies the HCS residue and arbitrary target residues from the
   command line. This is useful for mutants, noncanonical gating charges, control
   residues, drug-binding residues, or any residue set not encoded in the map.

Examples
--------
Analyze all mapped Nav1.7/SCN9A gating charges:

    python vsd_ztrack.py --psf step5_assembly.psf --traj traj.dcd \
        --gene SCN9A --vsd all --time-per-frame 0.24 --y-limits -20 20

Analyze only user-selected residues relative to a manual HCS:

    python vsd_ztrack.py --psf step5_assembly.psf --traj traj.dcd \
        --hcs-residue 163 --target-residues 214 217 220 223 \
        --target-labels R1 R2 R3 R4 --time-per-frame 0.24

Analyze selected residues while still using the mapped HCS for SCN9A VSD1:

    python vsd_ztrack.py --psf step5_assembly.psf --traj traj.dcd \
        --gene SCN9A --vsd VSD1 --target-residues 214 217 220 223

Author: Eslam Elhanafy and contributors
Suggested license: MIT or BSD-3-Clause
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import matplotlib

# Headless backend for HPC/server use.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from MDAnalysis import Universe
from MDAnalysis.analysis import align

R_COLUMNS = ("r1", "r2", "r3", "r4", "r5", "r6")
REQUIRED_COLUMNS = (
    "gene",
    "accession",
    "family",
    "vsd_label",
    "hcs_residue",
    *R_COLUMNS,
)
AXIS_TO_INDEX = {"x": 0, "y": 1, "z": 2}
DEFAULT_TRACE_COLORS = ["#ba508d", "#0de3ff", "#f72dd7", "#f2a533", "#00cc99", "#ff9900"]
DEFAULT_RESIDUE_MAP = Path(__file__).resolve().parent / "data" / "vsd_residue_map_with_accessions.csv"


@dataclass(frozen=True)
class MeasurementTarget:
    """One target-residue position measured relative to one HCS reference."""

    family: str
    gene: str
    accession: str
    vsd_label: str
    target_label: str
    target_order: int
    target_resid_original: int
    target_resid_selected: int
    hcs_resid_original: int
    hcs_resid_selected: int
    target_atoms: object
    hcs_atoms: object
    target_reference: str
    target_atom_name: str
    hcs_reference: str
    hcs_atom_name: str


def parse_residue(value: object) -> Optional[int]:
    """Parse one residue value from the residue map."""
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text == "" or text == "-" or text.lower() in {"nan", "na", "none", "null"}:
        return None
    try:
        return int(float(text))
    except ValueError as exc:
        raise ValueError(f"Could not parse residue value {value!r}") from exc


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize column names to lower-case labels."""
    return df.rename(columns={col: str(col).strip().lower() for col in df.columns})


def load_residue_map(csv_path: Path) -> pd.DataFrame:
    """Load and validate the VSD residue map."""
    if not csv_path.exists():
        raise FileNotFoundError(f"Residue-map file not found: {csv_path}")

    df = normalize_columns(pd.read_csv(csv_path))
    missing = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError("Residue-map CSV is missing required columns: " + ", ".join(missing))

    for col in ("gene", "accession", "family", "vsd_label"):
        df[col] = df[col].astype(str).str.strip()

    df["hcs_residue"] = df["hcs_residue"].apply(parse_residue).astype("Int64")
    for col in R_COLUMNS:
        df[col] = df[col].apply(parse_residue).astype("Int64")

    if df["hcs_residue"].isna().any():
        bad = df[df["hcs_residue"].isna()][["gene", "accession", "vsd_label"]]
        raise ValueError(f"HCS residue is missing for these entries:\n{bad}")

    return df


def available_entries_table(df: pd.DataFrame) -> pd.DataFrame:
    """Return a compact table of available residue-map entries."""
    cols = ["family", "gene", "accession", "vsd_label", "hcs_residue", *R_COLUMNS]
    return df[cols].sort_values(["family", "gene", "vsd_label"]).reset_index(drop=True)


def _case_match(series: pd.Series, value: Optional[str]) -> pd.Series:
    """Case-insensitive exact match helper."""
    if value is None:
        return pd.Series(True, index=series.index)
    return series.astype(str).str.upper() == value.upper()


def filter_residue_map(
    df: pd.DataFrame,
    family: Optional[str],
    gene: Optional[str],
    accession: Optional[str],
    vsd: str,
) -> pd.DataFrame:
    """Filter residue-map rows by family, gene, accession, and VSD label."""
    mask = (
        _case_match(df["family"], family)
        & _case_match(df["gene"], gene)
        & _case_match(df["accession"], accession)
    )
    if vsd.lower() != "all":
        mask &= _case_match(df["vsd_label"], vsd)

    selected = df[mask].copy()
    if selected.empty:
        raise ValueError(
            "No residue-map entries matched your filters. Use --list-map to inspect "
            "valid family/gene/accession/VSD combinations."
        )
    return selected.sort_values(["family", "gene", "vsd_label"]).reset_index(drop=True)


def manual_entry_from_args(args: argparse.Namespace) -> pd.DataFrame:
    """Create a one-row pseudo residue-map entry for fully manual analyses."""
    if args.hcs_residue is None:
        raise ValueError("Manual mode requires --hcs-residue.")
    if not args.target_residues:
        raise ValueError("Manual mode requires --target-residues.")

    vsd_label = args.vsd if args.vsd.lower() != "all" else args.custom_vsd_label
    row = {
        "family": args.family or "Custom",
        "gene": args.gene or "CUSTOM",
        "accession": args.accession or "NA",
        "vsd_label": vsd_label,
        "hcs_residue": int(args.hcs_residue),
    }
    for col in R_COLUMNS:
        row[col] = pd.NA
    return pd.DataFrame([row])


def parse_target_labels(args: argparse.Namespace, n_targets: int) -> list[str]:
    """Return target labels that match --target-residues."""
    if args.target_labels:
        if len(args.target_labels) != n_targets:
            raise ValueError(
                f"--target-labels contains {len(args.target_labels)} labels, but "
                f"--target-residues contains {n_targets} residues."
            )
        return [str(label) for label in args.target_labels]
    return [f"Residue_{resid}" for resid in args.target_residues]


def mapped_target_specs(row: pd.Series, args: argparse.Namespace) -> list[tuple[str, int]]:
    """Return target labels/residues from mapped R1-R6 columns."""
    requested = None
    if args.target_r_labels:
        requested = {label.lower().replace("-", "") for label in args.target_r_labels}

    specs: list[tuple[str, int]] = []
    for col in R_COLUMNS:
        if requested is not None and col.lower() not in requested:
            continue
        value = row[col]
        if pd.isna(value):
            continue
        specs.append((col.upper(), int(value)))
    return specs


def custom_target_specs(args: argparse.Namespace) -> list[tuple[str, int]]:
    """Return target labels/residues from --target-residues."""
    labels = parse_target_labels(args, len(args.target_residues))
    return [(label, int(resid)) for label, resid in zip(labels, args.target_residues)]


def build_residue_selection(
    resid_original: int,
    residue_selector: str,
    resid_offset: int,
    reference_mode: str,
    atom_name: str,
    no_hydrogens: bool,
    exclude_segid: Iterable[str],
    extra_selection: Optional[str],
) -> tuple[str, int]:
    """
    Build an MDAnalysis atom-selection string for a residue position.

    reference_mode controls which atoms represent the residue:
        ca            C-alpha atom only
        sidechain_com side-chain atoms; hydrogens optionally excluded
        residue_com   full residue; hydrogens optionally excluded
    """
    resid_selected = int(resid_original) + int(resid_offset)
    base = f"protein and {residue_selector} {resid_selected}"

    if reference_mode == "ca":
        selection = f"{base} and name {atom_name}"
    elif reference_mode == "sidechain_com":
        selection = f"{base} and not backbone"
        if no_hydrogens:
            selection += " and not name H*"
    elif reference_mode == "residue_com":
        selection = base
        if no_hydrogens:
            selection += " and not name H*"
    else:
        raise ValueError(f"Unsupported residue reference mode: {reference_mode!r}")

    for segid in exclude_segid:
        if segid:
            selection += f" and not segid {segid}"

    if extra_selection:
        selection += f" and ({extra_selection})"

    return selection, resid_selected


def residue_instances(atomgroup) -> set[tuple[str, int, str]]:
    """Return unique residue instances represented by an AtomGroup."""
    instances = set()
    for residue in atomgroup.residues:
        instances.add((str(residue.segid), int(residue.resid), str(residue.resname)))
    return instances


def select_or_record_missing(universe: Universe, selection: str, label: str, missing_messages: list[str]):
    """Select atoms and record missing selections for later reporting."""
    atoms = universe.select_atoms(selection)
    if len(atoms) == 0:
        missing_messages.append(f"{label}: empty selection -> {selection}")
    return atoms


def prepare_measurement_targets(
    universe: Universe,
    entries: pd.DataFrame,
    args: argparse.Namespace,
) -> list[MeasurementTarget]:
    """Convert map/manual rows into MDAnalysis AtomGroup measurement targets."""
    targets: list[MeasurementTarget] = []
    missing_messages: list[str] = []

    for _, row in entries.iterrows():
        hcs_original = int(args.hcs_residue) if args.hcs_residue is not None else int(row["hcs_residue"])
        hcs_selection, hcs_selected = build_residue_selection(
            resid_original=hcs_original,
            residue_selector=args.residue_selector,
            resid_offset=args.resid_offset,
            reference_mode=args.hcs_reference,
            atom_name=args.hcs_atom_name,
            no_hydrogens=args.no_hydrogens,
            exclude_segid=args.exclude_segid,
            extra_selection=args.extra_selection,
        )
        hcs_atoms = select_or_record_missing(
            universe,
            hcs_selection,
            f"{row['gene']} {row['vsd_label']} HCS {hcs_original}",
            missing_messages,
        )

        if len(hcs_atoms) > 0:
            hcs_instances = residue_instances(hcs_atoms)
            if len(hcs_instances) > 1:
                logging.warning(
                    "HCS selection contains atoms from %d residue instances for %s %s HCS %s. "
                    "Use --extra-selection or --exclude-segid if this is not intended. Selection: %s",
                    len(hcs_instances),
                    row["gene"],
                    row["vsd_label"],
                    hcs_original,
                    hcs_selection,
                )
            if args.hcs_reference == "ca" and len(hcs_atoms) > 1:
                logging.warning(
                    "HCS CA selection contains %d atoms for %s %s HCS %s. "
                    "Use --extra-selection or --exclude-segid to restrict the protein. Selection: %s",
                    len(hcs_atoms),
                    row["gene"],
                    row["vsd_label"],
                    hcs_original,
                    hcs_selection,
                )

        if args.target_residues:
            specs = custom_target_specs(args)
        else:
            specs = mapped_target_specs(row, args)

        if not specs:
            raise ValueError(
                f"No target residues were available for {row['gene']} {row['vsd_label']}. "
                "Provide --target-residues or check the residue map."
            )

        for order, (target_label, target_original) in enumerate(specs, start=1):
            target_selection, target_selected = build_residue_selection(
                resid_original=target_original,
                residue_selector=args.residue_selector,
                resid_offset=args.resid_offset,
                reference_mode=args.target_reference,
                atom_name=args.target_atom_name,
                no_hydrogens=args.no_hydrogens,
                exclude_segid=args.exclude_segid,
                extra_selection=args.extra_selection,
            )
            target_atoms = select_or_record_missing(
                universe,
                target_selection,
                f"{row['gene']} {row['vsd_label']} target {target_label} residue {target_original}",
                missing_messages,
            )

            if len(target_atoms) > 0:
                target_instances = residue_instances(target_atoms)
                if len(target_instances) > 1:
                    logging.warning(
                        "Target selection contains atoms from %d residue instances for %s %s %s residue %s. "
                        "Use --extra-selection or --exclude-segid if this is not intended. Selection: %s",
                        len(target_instances),
                        row["gene"],
                        row["vsd_label"],
                        target_label,
                        target_original,
                        target_selection,
                    )

            targets.append(
                MeasurementTarget(
                    family=str(row["family"]),
                    gene=str(row["gene"]),
                    accession=str(row["accession"]),
                    vsd_label=str(row["vsd_label"]),
                    target_label=str(target_label),
                    target_order=order,
                    target_resid_original=int(target_original),
                    target_resid_selected=int(target_selected),
                    hcs_resid_original=int(hcs_original),
                    hcs_resid_selected=int(hcs_selected),
                    target_atoms=target_atoms,
                    hcs_atoms=hcs_atoms,
                    target_reference=args.target_reference,
                    target_atom_name=args.target_atom_name,
                    hcs_reference=args.hcs_reference,
                    hcs_atom_name=args.hcs_atom_name,
                )
            )

    if missing_messages:
        message = "\n".join(missing_messages)
        if args.skip_missing:
            logging.warning("Some selections were empty and will be skipped:\n%s", message)
            targets = [
                target
                for target in targets
                if len(target.target_atoms) > 0 and len(target.hcs_atoms) > 0
            ]
        else:
            raise ValueError(
                "At least one atom selection was empty. This usually indicates a residue-numbering "
                "or segment-selection mismatch. Use --resid-offset, --residue-selector, "
                "--extra-selection, --exclude-segid, --no-sidechain-only, or --skip-missing.\n\n"
                + message
            )

    if not targets:
        raise ValueError("No valid targets remain after selection filtering.")

    return targets


def trajectory_time(ts, args: argparse.Namespace) -> float:
    """Return output time for a trajectory frame."""
    if args.time_per_frame is not None:
        return float(ts.frame) * float(args.time_per_frame)
    raw_time = getattr(ts, "time", None)
    if raw_time is None or not np.isfinite(raw_time):
        return float(ts.frame)
    return float(raw_time) * float(args.time_scale)


def analyze_trajectory(
    universe: Universe,
    targets: list[MeasurementTarget],
    args: argparse.Namespace,
    reference_universe: Optional[Universe] = None,
) -> pd.DataFrame:
    """Iterate over the trajectory and calculate target-HCS relative positions."""
    axis_index = AXIS_TO_INDEX[args.axis]
    rows = []
    trajectory_slice = universe.trajectory[args.start : args.stop : args.step]

    if args.align_protein and reference_universe is None:
        raise ValueError("--align-protein requires a prepared reference universe.")
    align_weights = None if args.align_weights == "none" else args.align_weights

    for ts in trajectory_slice:
        time_value = trajectory_time(ts, args)
        alignment_rmsd_before = np.nan
        alignment_rmsd_after = np.nan

        if args.align_protein:
            alignment_rmsd_before, alignment_rmsd_after = align.alignto(
                mobile=universe,
                reference=reference_universe,
                select=args.align_selection,
                weights=align_weights,
                match_atoms=True,
            )

        for target in targets:
            hcs_axis = target.hcs_atoms.center_of_mass()[axis_index]
            target_axis = target.target_atoms.center_of_mass()[axis_index]
            relative_position = target_axis - hcs_axis

            rows.append(
                {
                    "frame": int(ts.frame),
                    "time": time_value,
                    "time_unit": args.time_unit,
                    "axis": args.axis,
                    "family": target.family,
                    "gene": target.gene,
                    "accession": target.accession,
                    "vsd_label": target.vsd_label,
                    "hcs_resid_original": target.hcs_resid_original,
                    "hcs_resid_selected": target.hcs_resid_selected,
                    "target_label": target.target_label,
                    "target_order": target.target_order,
                    "target_resid_original": target.target_resid_original,
                    "target_resid_selected": target.target_resid_selected,
                    "target_reference": target.target_reference,
                    "target_atom_name": target.target_atom_name if target.target_reference == "ca" else "not_applicable",
                    "hcs_reference": target.hcs_reference,
                    "hcs_atom_name": target.hcs_atom_name if target.hcs_reference == "ca" else "not_applicable",
                    "hcs_reference_position_A": hcs_axis,
                    "target_position_A": target_axis,
                    "relative_position_A": relative_position,
                    # Backward-compatible aliases for earlier script versions.
                    "charge_label": target.target_label,
                    "charge_resid_original": target.target_resid_original,
                    "charge_resid_selected": target.target_resid_selected,
                    "charge_com_position_A": target_axis,
                    "charge_com_mode": target.target_reference,
                    "aligned_to_reference": bool(args.align_protein),
                    "align_selection": args.align_selection if args.align_protein else "none",
                    "align_reference_frame": args.align_reference_frame if args.align_protein else np.nan,
                    "alignment_rmsd_before_A": alignment_rmsd_before,
                    "alignment_rmsd_after_A": alignment_rmsd_after,
                }
            )

    if not rows:
        raise ValueError("No trajectory frames were analyzed. Check --start/--stop/--step.")
    return pd.DataFrame(rows)


def summarize_displacement(long_df: pd.DataFrame, avg_window: int) -> pd.DataFrame:
    """Summarize final-window minus initial-window relative displacement."""
    summary_rows = []
    grouping_cols = [
        "family",
        "gene",
        "accession",
        "vsd_label",
        "hcs_resid_original",
        "target_label",
        "target_resid_original",
    ]
    for keys, group in long_df.groupby(grouping_cols, sort=True):
        group = group.sort_values("frame")
        n = len(group)
        n_window = min(int(avg_window), n)
        first_mean = group["relative_position_A"].iloc[:n_window].mean()
        last_mean = group["relative_position_A"].iloc[-n_window:].mean()
        row = dict(zip(grouping_cols, keys))
        row.update(
            {
                "n_frames": n,
                "n_average_window": n_window,
                "initial_mean_A": first_mean,
                "final_mean_A": last_mean,
                "delta_final_minus_initial_A": last_mean - first_mean,
            }
        )
        summary_rows.append(row)
    return pd.DataFrame(summary_rows)


def parse_pair(values: Optional[list[float]]) -> Optional[tuple[float, float]]:
    """Parse a two-value range argument."""
    if values is None:
        return None
    if len(values) != 2:
        raise ValueError("Range arguments must contain exactly two values.")
    return float(values[0]), float(values[1])


def plot_vsd_timeseries(long_df: pd.DataFrame, args: argparse.Namespace) -> None:
    """Create one time-series plot per family/gene/accession/VSD group."""
    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    y_limits = parse_pair(args.y_limits)
    x_limits = parse_pair(args.x_limits)

    if args.font_family:
        plt.rcParams["font.family"] = args.font_family

    group_cols = ["family", "gene", "accession", "vsd_label"]
    trace_colors = args.colors if args.colors else DEFAULT_TRACE_COLORS
    prefix = args.prefix if args.prefix is not None else Path(args.traj).stem

    for keys, group in long_df.groupby(group_cols, sort=True):
        family, gene, accession, vsd_label = keys
        fig, ax = plt.subplots(figsize=(args.figure_width, args.figure_height), dpi=args.dpi)

        for idx, (target_label, target_group) in enumerate(
            sorted(
                group.groupby("target_label"),
                key=lambda item: int(item[1]["target_order"].iloc[0]),
            )
        ):
            target_group = target_group.sort_values("frame")
            resid = int(target_group["target_resid_original"].iloc[0])
            if args.legend_labels == "residue":
                label = str(resid)
            elif args.legend_labels == "target":
                label = str(target_label)
            else:
                label = f"{target_label} ({resid})"

            ax.plot(
                target_group["time"],
                target_group["relative_position_A"],
                linewidth=args.linewidth,
                label=label,
                color=trace_colors[idx % len(trace_colors)],
            )

        ax.axhline(0.0, linestyle="--", linewidth=args.zero_linewidth, color="black")

        if args.legacy_axis_labels:
            ylabel = f"{args.axis.upper()}_Position (Å)"
        else:
            if args.hcs_reference == "ca":
                ylabel = f"{args.axis.upper()}(target) - {args.axis.upper()}(HCS Cα) (Å)"
            else:
                ylabel = f"{args.axis.upper()}(target) - {args.axis.upper()}(HCS reference) (Å)"

        ax.set_xlabel(f"Time ({args.time_unit})", fontsize=args.axis_label_fontsize)
        ax.set_ylabel(ylabel, fontsize=args.axis_label_fontsize)

        if args.title:
            title = args.title.format(
                family=family,
                gene=gene,
                accession=accession,
                vsd=vsd_label,
                prefix=prefix,
                axis=args.axis.upper(),
            )
        else:
            title = f"{vsd_label}_{prefix}"
        ax.set_title(title, fontsize=args.title_fontsize, pad=args.title_pad)

        if y_limits:
            ax.set_ylim(*y_limits)
        if x_limits:
            ax.set_xlim(*x_limits)
        else:
            ax.set_xlim(float(group["time"].min()), float(group["time"].max()))

        ax.tick_params(axis="both", labelsize=args.tick_fontsize, width=args.spine_linewidth, length=6)
        for spine in ax.spines.values():
            spine.set_linewidth(args.spine_linewidth)

        if args.upper_left_label:
            ax.text(
                args.panel_label_x,
                args.panel_label_y,
                args.upper_left_label,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=args.panel_label_fontsize,
            )

        leg = ax.legend(
            frameon=False,
            bbox_to_anchor=(args.legend_x, args.legend_y),
            loc=args.legend_loc,
            fontsize=args.legend_fontsize,
            handlelength=args.legend_handlelength,
            labelspacing=args.legend_labelspacing,
            borderaxespad=0.0,
        )
        legend_handles = getattr(leg, "legend_handles", None)
        if legend_handles is None:
            legend_handles = getattr(leg, "legendHandles", [])
        for legobj in legend_handles:
            legobj.set_linewidth(args.legend_linewidth)

        fig.tight_layout(rect=(0.0, 0.0, args.tight_right, 1.0))

        hcs_suffix = "HCS_CA" if args.hcs_reference == "ca" else f"HCS_{args.hcs_reference}"
        descriptive_stem = f"{gene}_{accession}_{vsd_label}_{args.axis}_relative_to_{hcs_suffix}"
        legacy_stem = f"{vsd_label}_{args.axis.upper()}_position"

        for fmt in args.plot_formats:
            fig.savefig(output_dir / f"{descriptive_stem}.{fmt}", dpi=args.dpi)
            if args.save_legacy_plot_names:
                fig.savefig(output_dir / f"{legacy_stem}.{fmt}", dpi=args.dpi)
        plt.close(fig)


def write_run_metadata(args: argparse.Namespace, entries: pd.DataFrame, output_dir: Path) -> None:
    """Write run settings for reproducibility."""
    metadata_path = output_dir / "run_metadata.txt"
    with metadata_path.open("w", encoding="utf-8") as handle:
        handle.write("VSD-ZTrack analysis metadata\n")
        handle.write("============================\n\n")
        for key, value in sorted(vars(args).items()):
            handle.write(f"{key}: {value}\n")
        handle.write("\nSelected map/manual entries\n")
        handle.write("---------------------------\n")
        entries.to_csv(handle, index=False)


def build_arg_parser() -> argparse.ArgumentParser:
    """Create the command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Calculate target-residue z-position or x/y/z displacement relative to "
            "the HCS C-alpha or another HCS reference. Supports mapped VSD gating "
            "charges and user-selected residues."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--psf", required=False, type=Path, help="Topology file, e.g., PSF/PDB/GRO.")
    parser.add_argument("--traj", required=False, type=Path, help="Trajectory file, e.g., DCD/XTC/TRR.")
    parser.add_argument(
        "--residue-map",
        type=Path,
        default=DEFAULT_RESIDUE_MAP,
        help="CSV containing family/gene/accession/VSD/HCS/R1-R6 residue numbers.",
    )
    parser.add_argument("--family", default=None, help="Channel family filter, e.g., Nav, Cav, Kv.")
    parser.add_argument("--gene", default=None, help="Gene filter, e.g., SCN9A, CACNA1C, KCNQ1.")
    parser.add_argument("--accession", default=None, help="UniProt accession filter, e.g., Q15858.")
    parser.add_argument("--vsd", default="all", help="VSD label, e.g., VSD1, VSD2, VSD3, VSD4, or all.")
    parser.add_argument("--list-map", action="store_true", help="Print available residue-map entries and exit.")

    parser.add_argument(
        "--target-residues",
        nargs="+",
        type=int,
        default=None,
        help=(
            "User-selected residue numbers to analyze instead of mapped R1-R6 residues. "
            "Example: --target-residues 214 217 220 223."
        ),
    )
    parser.add_argument(
        "--target-labels",
        nargs="+",
        default=None,
        help=(
            "Optional labels for --target-residues. Must have the same number of entries. "
            "Example: --target-labels R1 R2 R3 R4."
        ),
    )
    parser.add_argument(
        "--target-r-labels",
        nargs="+",
        default=None,
        help="Subset mapped R labels without typing residue numbers, e.g., --target-r-labels R1 R2 R3.",
    )
    parser.add_argument(
        "--hcs-residue",
        type=int,
        default=None,
        help=(
            "Manual HCS residue number. Overrides the residue-map HCS. Required if no "
            "residue-map selection is used."
        ),
    )
    parser.add_argument(
        "--custom-vsd-label",
        default="Custom",
        help="VSD/domain label used for fully manual analyses when --vsd all.",
    )

    parser.add_argument(
        "--residue-selector",
        choices=("resid", "resnum"),
        default="resid",
        help="MDAnalysis residue-number selection keyword.",
    )
    parser.add_argument(
        "--resid-offset",
        type=int,
        default=0,
        help="Integer added to all input residue numbers before atom selection.",
    )
    parser.add_argument(
        "--hcs-reference",
        choices=("ca", "sidechain_com", "residue_com"),
        default="ca",
        help="Atoms used to represent the HCS reference residue.",
    )
    parser.add_argument("--hcs-atom-name", default="CA", help="Atom name used when --hcs-reference ca.")
    parser.add_argument(
        "--target-reference",
        choices=("sidechain_com", "residue_com", "ca"),
        default="sidechain_com",
        help="Atoms used to represent each target residue.",
    )
    parser.add_argument("--target-atom-name", default="CA", help="Atom name used when --target-reference ca.")
    parser.add_argument(
        "--sidechain-only",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Backward-compatible alias for target representation. If explicitly set to "
            "--no-sidechain-only and --target-reference is unchanged, target-reference "
            "is converted to residue_com."
        ),
    )
    parser.add_argument(
        "--no-hydrogens",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Exclude hydrogen atoms from COM calculations.",
    )
    parser.add_argument("--exclude-segid", nargs="*", default=[], help="One or more segment IDs to exclude.")
    parser.add_argument("--extra-selection", default=None, help='Extra MDAnalysis selection text, e.g., "segid PROA".')
    parser.add_argument("--skip-missing", action="store_true", help="Skip empty selections instead of stopping.")

    parser.add_argument(
        "--align-protein",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Align each frame to a reference frame before analysis.",
    )
    parser.add_argument("--align-selection", default="protein and name CA", help="Selection used for alignment.")
    parser.add_argument("--align-reference-frame", type=int, default=0, help="Reference frame for alignment.")
    parser.add_argument("--align-weights", choices=("mass", "none"), default="mass", help="Alignment weights.")

    parser.add_argument("--axis", choices=("x", "y", "z"), default="z", help="Coordinate axis used for displacement.")
    parser.add_argument("--time-per-frame", type=float, default=None, help="Output time increment per frame.")
    parser.add_argument("--time-scale", type=float, default=0.001, help="Scale applied to MDAnalysis ts.time.")
    parser.add_argument("--time-unit", default="ns", help="Output time-unit label.")
    parser.add_argument("--start", type=int, default=None, help="First trajectory frame to analyze.")
    parser.add_argument("--stop", type=int, default=None, help="Stop before this trajectory frame.")
    parser.add_argument("--step", type=int, default=1, help="Trajectory stride.")
    parser.add_argument("--avg-window", type=int, default=100, help="First/last analyzed frames used for summary.")

    parser.add_argument("--outdir", default="vsd_ztrack_results", help="Output directory.")
    parser.add_argument("--prefix", default=None, help="Output filename prefix. Default uses trajectory stem.")
    parser.add_argument("--plot-formats", nargs="+", default=["png", "svg"], help="Plot formats.")
    parser.add_argument("--dpi", type=int, default=100, help="Raster image resolution.")
    parser.add_argument("--figure-width", type=float, default=20.48, help="Figure width in inches.")
    parser.add_argument("--figure-height", type=float, default=6.14, help="Figure height in inches.")
    parser.add_argument("--linewidth", type=float, default=1.6, help="Trace linewidth.")
    parser.add_argument("--zero-linewidth", type=float, default=1.2, help="Zero reference-line width.")
    parser.add_argument("--colors", nargs="+", default=DEFAULT_TRACE_COLORS, help="Trace colors.")
    parser.add_argument("--font-family", default="DejaVu Sans", help="Matplotlib font family.")
    parser.add_argument("--title-fontsize", type=float, default=24, help="Plot title font size.")
    parser.add_argument("--title-pad", type=float, default=12, help="Title padding.")
    parser.add_argument("--axis-label-fontsize", type=float, default=18, help="Axis-label font size.")
    parser.add_argument("--tick-fontsize", type=float, default=16, help="Tick-label font size.")
    parser.add_argument("--panel-label-fontsize", type=float, default=18, help="Upper-left panel-label font size.")
    parser.add_argument("--legend-fontsize", type=float, default=18, help="Legend font size.")
    parser.add_argument("--legend-linewidth", type=float, default=2.5, help="Legend trace linewidth.")
    parser.add_argument("--legend-handlelength", type=float, default=2.2, help="Legend handle length.")
    parser.add_argument("--legend-labelspacing", type=float, default=1.1, help="Vertical legend-label spacing.")
    parser.add_argument("--legend-x", type=float, default=1.015, help="Legend x-position in axes coordinates.")
    parser.add_argument("--legend-y", type=float, default=0.72, help="Legend y-position in axes coordinates.")
    parser.add_argument("--legend-loc", default="center left", help="Matplotlib legend location.")
    parser.add_argument("--tight-right", type=float, default=0.91, help="Right boundary for tight_layout rectangle.")
    parser.add_argument("--spine-linewidth", type=float, default=1.2, help="Axis spine and tick linewidth.")
    parser.add_argument("--panel-label-x", type=float, default=0.0, help="Upper-left label x-position.")
    parser.add_argument("--panel-label-y", type=float, default=1.0, help="Upper-left label y-position.")
    parser.add_argument(
        "--legend-labels",
        choices=("residue", "target", "both"),
        default="residue",
        help="Legend text style.",
    )
    parser.add_argument(
        "--legacy-axis-labels",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use labels such as Z_Position (Å), matching the original plots.",
    )
    parser.add_argument(
        "--save-legacy-plot-names",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Also save plots as VSDx_Z_position.<format>.",
    )
    parser.add_argument("--x-limits", nargs=2, type=float, default=None, metavar=("XMIN", "XMAX"))
    parser.add_argument("--y-limits", nargs=2, type=float, default=None, metavar=("YMIN", "YMAX"))
    parser.add_argument("--upper-left-label", default=None, help="Optional panel label.")
    parser.add_argument(
        "--title",
        default=None,
        help="Optional title template with {family}, {gene}, {accession}, {vsd}, {prefix}, {axis}.",
    )
    parser.add_argument("--no-plots", action="store_true", help="Write CSV outputs only; do not make figures.")
    parser.add_argument("--verbose", action="store_true", help="Print additional progress information.")
    return parser


def resolve_entries(args: argparse.Namespace) -> pd.DataFrame:
    """Resolve residue-map or fully manual entries from CLI options."""
    # Fully manual mode: no map filters supplied, but HCS and targets are supplied.
    using_map_filters = any([args.family, args.gene, args.accession, args.vsd.lower() != "all"])
    fully_manual = args.hcs_residue is not None and args.target_residues and not using_map_filters

    if fully_manual:
        return manual_entry_from_args(args)

    if not args.residue_map.exists():
        if args.hcs_residue is not None and args.target_residues:
            return manual_entry_from_args(args)
        raise FileNotFoundError(
            f"Residue map not found: {args.residue_map}. Provide --residue-map, or use "
            "manual mode with --hcs-residue and --target-residues."
        )

    residue_map = load_residue_map(args.residue_map)
    if args.list_map:
        table = available_entries_table(residue_map)
        with pd.option_context("display.max_rows", None, "display.max_columns", None):
            print(table.to_string(index=False))
        raise SystemExit(0)

    return filter_residue_map(
        residue_map,
        family=args.family,
        gene=args.gene,
        accession=args.accession,
        vsd=args.vsd,
    )


def main(argv: Optional[list[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )

    if args.psf is None or args.traj is None:
        if args.list_map:
            if not args.residue_map.exists():
                parser.error(f"Residue map not found: {args.residue_map}")
            residue_map = load_residue_map(args.residue_map)
            table = available_entries_table(residue_map)
            with pd.option_context("display.max_rows", None, "display.max_columns", None):
                print(table.to_string(index=False))
            return 0
        parser.error("--psf and --traj are required unless --list-map is used.")

    if args.target_labels and not args.target_residues:
        parser.error("--target-labels can only be used with --target-residues.")

    if args.target_r_labels and args.target_residues:
        parser.error("Use either --target-r-labels for mapped R1-R6 or --target-residues, not both.")

    # Backward compatibility with older flag. argparse does not tell whether the
    # user explicitly set it, so convert only for the common requested case.
    if args.sidechain_only is False and args.target_reference == "sidechain_com":
        args.target_reference = "residue_com"

    if not args.psf.exists():
        raise FileNotFoundError(f"Topology file not found: {args.psf}")
    if not args.traj.exists():
        raise FileNotFoundError(f"Trajectory file not found: {args.traj}")

    entries = resolve_entries(args)

    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Loading Universe: topology=%s trajectory=%s", args.psf, args.traj)
    universe = Universe(str(args.psf), str(args.traj))

    reference_universe = None
    if args.align_protein:
        logging.info("Preparing alignment reference frame %s", args.align_reference_frame)
        reference_universe = Universe(str(args.psf), str(args.traj))
        n_ref_frames = len(reference_universe.trajectory)
        if args.align_reference_frame < 0 or args.align_reference_frame >= n_ref_frames:
            raise ValueError(
                f"--align-reference-frame must be between 0 and {n_ref_frames - 1}; "
                f"received {args.align_reference_frame}."
            )
        reference_universe.trajectory[args.align_reference_frame]
        n_mobile = len(universe.select_atoms(args.align_selection))
        n_reference = len(reference_universe.select_atoms(args.align_selection))
        if n_mobile == 0 or n_reference == 0:
            raise ValueError(f"Alignment selection is empty: {args.align_selection!r}")
        if n_mobile != n_reference:
            raise ValueError(
                "Alignment selection produced different atom counts in mobile and reference universes: "
                f"mobile={n_mobile}, reference={n_reference}."
            )

    targets = prepare_measurement_targets(universe, entries, args)
    logging.info("Prepared %d target-HCS measurements.", len(targets))

    long_df = analyze_trajectory(universe, targets, args, reference_universe=reference_universe)
    summary_df = summarize_displacement(long_df, args.avg_window)

    prefix = args.prefix if args.prefix is not None else Path(args.traj).stem
    long_path = output_dir / f"{prefix}_vsd_ztrack_positions_long.csv"
    summary_path = output_dir / f"{prefix}_vsd_ztrack_summary.csv"
    selected_path = output_dir / f"{prefix}_selected_entries.csv"

    long_df.to_csv(long_path, index=False)
    summary_df.to_csv(summary_path, index=False)
    entries.to_csv(selected_path, index=False)
    write_run_metadata(args, entries, output_dir)

    if not args.no_plots:
        plot_vsd_timeseries(long_df, args)

    print(f"Done. Wrote: {long_path}")
    print(f"Done. Wrote: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
