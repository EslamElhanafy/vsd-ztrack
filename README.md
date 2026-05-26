# VSD-ZTrack

**VSD-ZTrack** is a command-line MDAnalysis tool for measuring residue-level displacement along the membrane-normal axis during molecular dynamics simulations.

The default analysis is designed for voltage-sensor domains (VSDs): it measures the position of each S4 gating-charge residue relative to the hydrophobic constriction site (HCS). The script can also analyze any user-defined residue set, so it is not restricted to canonical R1-R6 gating charges.

## Recommended script/repository name

Suggested GitHub repository name:

```text
vsd-ztrack
```

Suggested script name:

```text
vsd_ztrack.py
```

Rationale: the name is short, searchable, and describes the core analysis: tracking residue displacement along the z-axis in voltage-sensor domains.

## Core metric

For each trajectory frame and each target residue:

```text
relative_position_A = z(target residue) - z(HCS reference)
```

By default:

```text
HCS reference = C-alpha atom of the HCS residue
Target residue = side-chain heavy-atom center of mass
Axis = z
```

This gives:

```text
Z-distance = Z(COM of selected target residue side chain) - Z(CA of HCS)
```

The metric is useful for comparing the axial position of S4 gating charges or other residues relative to the HCS during equilibrium MD, SMD, or voltage-dependent simulations.

## Main features

- Supports mapped VSD analyses using `data/vsd_residue_map_with_accessions.csv`.
- Supports Nav, Cav, Kv, CNG, HCN, Hv, and TRPP entries included in the residue map.
- Supports arbitrary user-selected residues through `--target-residues`.
- Supports mapped R-label subsetting through `--target-r-labels R1 R2 R3`.
- Measures target residues relative to HCS Cα by default.
- Optional protein alignment before analysis using `--align-protein`.
- Optional segment filtering through `--extra-selection` or `--exclude-segid`.
- Outputs long-form CSV, summary CSV, selected-entry CSV, metadata, and publication-style plots.
- Works in non-interactive environments such as HPC clusters.

## Installation

### Option 1: Python virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

### Option 2: conda or mamba

```bash
mamba create -n vsd-ztrack python=3.11 mdanalysis numpy pandas matplotlib -c conda-forge
mamba activate vsd-ztrack
```

or:

```bash
conda create -n vsd-ztrack python=3.11 mdanalysis numpy pandas matplotlib -c conda-forge
conda activate vsd-ztrack
```

## Repository structure

```text
vsd-ztrack/
├── README.md
├── LICENSE
├── requirements.txt
├── .gitignore
├── vsd_ztrack.py
├── data/
│   └── vsd_residue_map_with_accessions.csv
└── examples/
    ├── run_nav17_all_gating_charges.sh
    ├── run_nav17_custom_residues.sh
    └── run_nav17_custom_residues_aligned.sh
```

## Quick start

### 1. List available channels and VSDs

```bash
python vsd_ztrack.py --list-map
```

### 2. Run mapped Nav1.7 / SCN9A analysis for all VSDs

```bash
cd /media/eslam/18T/cationpi/20260501_Nav7_down/y163l_f787l_9mut/04_sim/prod_anton3

python ~/Downloads/vsd-ztrack/vsd_ztrack.py \
  --psf step5_assembly.psf \
  --traj 20260525_nav7_y163l_f787l_9mut_-200mV_3us.dcd \
  --gene SCN9A \
  --vsd all \
  --time-per-frame 0.24 \
  --time-unit ns \
  --y-limits -20 20 \
  --prefix 20260525_nav7_y163l_f787l_9mut_-200mV_3us \
  --upper-left-label nav7_y163l_f787l_9mut_-200mV_3us \
  --outdir vsd_ztrack_SCN9A_all
```

This uses the mapped HCS and mapped R1-R6 residues for each SCN9A VSD.

## New user-selected residue option

Use `--target-residues` when you want to choose the residues manually instead of using the mapped R1-R6 residues.

### Example: measure selected Nav1.7 VSD1 residues relative to HCS 163

```bash
cd /media/eslam/18T/cationpi/20260501_Nav7_down/y163l_f787l_9mut/04_sim/prod_anton3

python ~/Downloads/vsd-ztrack/vsd_ztrack.py \
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
```

This measures:

```text
Z(COM of residue 214 side chain) - Z(CA of HCS 163)
Z(COM of residue 217 side chain) - Z(CA of HCS 163)
Z(COM of residue 220 side chain) - Z(CA of HCS 163)
Z(COM of residue 223 side chain) - Z(CA of HCS 163)
```

The labels `R1 R2 R3 R4` are optional but useful for legends and CSV outputs.

### Example: use mapped HCS, but manually choose target residues

```bash
python ~/Downloads/vsd-ztrack/vsd_ztrack.py \
  --psf step5_assembly.psf \
  --traj 20260525_nav7_y163l_f787l_9mut_-200mV_3us.dcd \
  --gene SCN9A \
  --vsd VSD1 \
  --target-residues 214 217 220 223 \
  --target-labels R1 R2 R3 R4 \
  --time-per-frame 0.24 \
  --time-unit ns \
  --outdir vsd_ztrack_SCN9A_VSD1_custom_targets
```

Here, the HCS is taken from the residue map for SCN9A VSD1, while the target residues are supplied by the user.

### Example: analyze only selected mapped R-labels

```bash
python ~/Downloads/vsd-ztrack/vsd_ztrack.py \
  --psf step5_assembly.psf \
  --traj 20260525_nav7_y163l_f787l_9mut_-200mV_3us.dcd \
  --gene SCN9A \
  --vsd VSD1 \
  --target-r-labels R1 R2 R3 \
  --time-per-frame 0.24 \
  --time-unit ns \
  --outdir vsd_ztrack_SCN9A_VSD1_R1_R2_R3
```

This uses the residue map but measures only R1-R3.

## Alignment option

Use alignment when the protein drifts or rotates in the simulation box.

```bash
python ~/Downloads/vsd-ztrack/vsd_ztrack.py \
  --psf step5_assembly.psf \
  --traj 20260525_nav7_y163l_f787l_9mut_-200mV_3us.dcd \
  --hcs-residue 163 \
  --target-residues 214 217 220 223 \
  --target-labels R1 R2 R3 R4 \
  --family Nav \
  --gene SCN9A \
  --vsd VSD1 \
  --align-protein \
  --align-selection "protein and name CA" \
  --align-reference-frame 0 \
  --time-per-frame 0.24 \
  --time-unit ns \
  --outdir vsd_ztrack_SCN9A_VSD1_custom_residues_aligned
```

If the system has duplicate protein segments, restrict the analysis first. For example:

```bash
--extra-selection "segid PROA"
```

or exclude a segment:

```bash
--exclude-segid PROB
```

`--extra-selection` is usually safer because it explicitly selects the protein segment to analyze.

## Residue representation options

### Default target representation

```bash
--target-reference sidechain_com
```

This uses the side-chain heavy-atom center of mass.

### Whole-residue target COM

```bash
--target-reference residue_com
```

or the older compatibility flag:

```bash
--no-sidechain-only
```

### Target Cα only

```bash
--target-reference ca --target-atom-name CA
```

### HCS side-chain COM instead of HCS Cα

```bash
--hcs-reference sidechain_com
```

### HCS whole-residue COM

```bash
--hcs-reference residue_com
```

## Important time-factor check

For a 3 µs trajectory, the time-per-frame value depends on the number of frames in the DCD.

```bash
python - <<'PY'
import MDAnalysis as mda
u = mda.Universe("step5_assembly.psf", "20260525_nav7_y163l_f787l_9mut_-200mV_3us.dcd")
print("Number of frames:", len(u.trajectory))
print("time_per_frame for 3000 ns:", 3000 / (len(u.trajectory) - 1), "ns/frame")
PY
```

Use the printed value with:

```bash
--time-per-frame VALUE
```

## Outputs

The output directory contains:

```text
<prefix>_vsd_ztrack_positions_long.csv
<prefix>_vsd_ztrack_summary.csv
<prefix>_selected_entries.csv
run_metadata.txt
<gene>_<accession>_<VSD>_z_relative_to_HCS_CA.png
<gene>_<accession>_<VSD>_z_relative_to_HCS_CA.svg
VSD1_Z_position.png
VSD1_Z_position.svg
```

### Long-form CSV columns

Important columns include:

```text
frame
time
time_unit
axis
family
gene
accession
vsd_label
hcs_resid_original
hcs_resid_selected
target_label
target_resid_original
target_resid_selected
target_reference
hcs_reference
hcs_reference_position_A
target_position_A
relative_position_A
aligned_to_reference
alignment_rmsd_before_A
alignment_rmsd_after_A
```

### Summary CSV columns

The summary file reports the first-window and last-window means:

```text
initial_mean_A
final_mean_A
delta_final_minus_initial_A
```

The averaging window is controlled by:

```bash
--avg-window 100
```

## Troubleshooting

### Empty residue selection

If you see an error such as:

```text
empty selection -> protein and resid 163 and name CA
```

check whether the residue numbering in the topology matches the residue map. Possible fixes:

```bash
--resid-offset -2
```

or:

```bash
--residue-selector resnum
```

### Duplicate HCS Cα atom warning

If you see:

```text
HCS CA selection contains 2 atoms
```

your topology likely has more than one protein segment with the same residue number. Identify the segment:

```bash
python - <<'PY'
import MDAnalysis as mda
u = mda.Universe("step5_assembly.psf", "20260525_nav7_y163l_f787l_9mut_-200mV_3us.dcd")
sel = u.select_atoms("protein and resid 163 and name CA")
for a in sel:
    print(a.index, a.segid, a.resid, a.resname, a.name)
PY
```

Then rerun with:

```bash
--extra-selection "segid YOUR_SEGID"
```

## How to upload this package to GitHub

### 1. Create a new GitHub repository

Create an empty repository named:

```text
vsd-ztrack
```

Do not initialize it with a README if you are uploading this local folder, because this package already includes a README.

### 2. Prepare the local folder

```bash
cd ~/Downloads/vsd-ztrack
```

Check the files:

```bash
ls -la
```

### 3. Initialize Git

```bash
git init
git add README.md LICENSE requirements.txt .gitignore vsd_ztrack.py data examples
git commit -m "Initial release of VSD-ZTrack"
git branch -M main
```

### 4. Connect to your GitHub repository

Replace `YOUR_USERNAME` with your GitHub username:

```bash
git remote add origin https://github.com/YOUR_USERNAME/vsd-ztrack.git
```

### 5. Push to GitHub

```bash
git push -u origin main
```

### 6. Future updates

After editing the script or README:

```bash
git status
git add vsd_ztrack.py README.md
git commit -m "Add custom target residue analysis option"
git push
```

## Recommended repository description

```text
Residue-map-driven MDAnalysis tool for tracking VSD gating-charge and custom residue displacement relative to the hydrophobic constriction site in ion-channel simulations.
```

## Citation suggestion

If this tool supports a publication, cite the GitHub repository and specify the version or commit hash used for the analysis.
