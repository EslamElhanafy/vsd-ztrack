#!/usr/bin/env bash
set -euo pipefail

cd /media/eslam/18T/cationpi/20260501_Nav7_down/y163l_f787l_9mut/04_sim/prod_anton3

python /path/to/vsd-ztrack/vsd_ztrack.py \
  --psf step5_assembly.psf \
  --traj 20260525_nav7_y163l_f787l_9mut_-200mV_3us.dcd \
  --hcs-residue 163 \
  --target-residues 214 217 220 223 \
  --target-labels R1 R2 R3 R4 \
  --family Nav \
  --gene SCN9A \
  --vsd VSD1 \
  --time-per-frame 0.24 \
  --time-unit ns \
  --y-limits -20 20 \
  --prefix 20260525_nav7_y163l_f787l_9mut_-200mV_3us \
  --upper-left-label nav7_y163l_f787l_9mut_-200mV_3us \
  --outdir vsd_ztrack_SCN9A_VSD1_custom_residues
