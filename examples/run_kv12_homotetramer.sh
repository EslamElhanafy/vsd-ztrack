#!/usr/bin/env bash
# Kv1.2 (KCNA2) homotetramer — analyze all four subunits in one run.
#
# Kv channels are homotetramers: each subunit (PROA–PROD) carries one VSD
# with the same gating-charge residues.  --segments processes each subunit
# independently, scoping both the HCS and target atom selections to that
# segment so the four subunits never interfere with each other.
#
# Residues are taken from the built-in residue map (family Kv, gene KCNA2):
#   HCS: I233   R1: R294   R2: R297   R3: R300   R4: R303
#
# Outputs (in vsd_ztrack_KCNA2_homotetramer/):
#   *_vsd_ztrack_positions_long.csv  — per-frame z-positions, one row per
#                                       (frame, segment, gating charge)
#   *_vsd_ztrack_summary.csv         — initial/final mean and displacement
#   KCNA2_P16389_VSD1_PROA_z_relative_to_HCS_CA.{png,svg}  (x4 segments)
#   run_metadata.txt
set -euo pipefail

cd /path/to/simulation/directory

python /path/to/vsd-ztrack/vsd_ztrack.py \
  --psf step5_assembly.psf \
  --traj production.dcd \
  --gene KCNA2 \
  --vsd VSD1 \
  --segments PROA PROB PROC PROD \
  --time-per-frame 0.24 \
  --time-unit ns \
  --y-limits -20 20 \
  --outdir vsd_ztrack_KCNA2_homotetramer
