#!/usr/bin/env bash
# TRPV1 homotetramer — custom residues, all four subunits in one run.
#
# TRP channels are homotetramers: each subunit (PROA–PROD) carries one
# voltage-sensing-like domain.  --segments processes each subunit
# independently so HCS and target selections are always scoped to the
# correct chain.
#
# Because TRPV1 is not in the built-in Kv/Nav/Cav residue map, the HCS
# and target residues are supplied manually via --hcs-residue and
# --target-residues.  Replace the residue numbers below with those
# matching your sequence/numbering (human TRPV1, UniProt Q8NER1,
# cryo-EM PDB 3J5P as reference):
#
#   HCS:  I679  (hydrophobic gate, S4–S5 linker)
#   S4 gating-charge equivalents: R557  R560  R563  K571
#
# Outputs (in vsd_ztrack_TRPV1_homotetramer/):
#   *_vsd_ztrack_positions_long.csv
#   *_vsd_ztrack_summary.csv
#   CUSTOM_NA_Custom_PROA_z_relative_to_HCS_CA.{png,svg}  (x4 segments)
#   run_metadata.txt
set -euo pipefail

cd /path/to/simulation/directory

python /path/to/vsd-ztrack/vsd_ztrack.py \
  --psf step5_assembly.psf \
  --traj production.dcd \
  --hcs-residue 679 \
  --target-residues 557 560 563 571 \
  --target-labels R1 R2 R3 R4 \
  --family TRP \
  --gene TRPV1 \
  --vsd VSD1 \
  --segments PROA PROB PROC PROD \
  --time-per-frame 0.24 \
  --time-unit ns \
  --y-limits -20 20 \
  --outdir vsd_ztrack_TRPV1_homotetramer
