# VSD-ZTrack

VSD-ZTrack is a command-line analysis tool for tracking residue displacement along a selected Cartesian axis in molecular-dynamics trajectories. It was designed for voltage-sensor domains, where S4 gating-charge residues are often compared with the hydrophobic constriction site, but it can also analyze any user-defined residue set.

<img width="1536" height="864" alt="9bd0e7b5-58f1-4bdf-9ec4-68afa79bba16-0" src="https://github.com/user-attachments/assets/d2ca6b9b-6565-4a78-b6a6-4bc1433d4e52" />

The default analysis reports the position of each target residue relative to a reference residue:

```text
relative_position = coordinate(target residue) - coordinate(reference residue)
```

For voltage-sensor applications, the reference residue is usually the hydrophobic constriction site, and the targets are usually S4 gating-charge residues. By default, the reference is represented by its C-alpha atom, and each target residue is represented by the side-chain heavy-atom center of mass.

## Main use cases

- Track S4 gating-charge displacement relative to the hydrophobic constriction site.
- Compare voltage-sensor conformational changes across channel families.
- Analyze Nav, Cav, Kv, HCN, CNG, Hv, TRPP, or other VSD-containing proteins.
- Measure custom residues defined by the user.
- Run reproducible non-interactive analyses on workstations or HPC clusters.

## Features

- Residue-map-driven analysis for mapped VSD systems.
- Fully manual mode for any reference residue and target-residue list.
- Optional analysis of all VSDs for one selected channel.
- Optional subsetting of mapped gating-charge labels.
- Optional protein alignment before the distance calculation.
- Optional segment filtering for multichain or duplicated-protein systems.
- Configurable coordinate axis, time scaling, frame stride, averaging window, plot size, colors, and output formats.
- Long-form and summary data outputs for downstream plotting or statistical analysis.

## Installation

Create a clean Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Alternatively, use conda or mamba:

```bash
mamba create -n vsd-ztrack python=3.11 mdanalysis numpy pandas matplotlib -c conda-forge
mamba activate vsd-ztrack
```

## Repository layout

```text
vsd-ztrack/
├── README
├── LICENSE
├── requirements
├── analysis script
├── data/
│   └── residue map
├── example_data/
│   ├── example residue-map subset
│   ├── example target definitions
│   ├── example long-form output
│   └── example summary output
└── examples/
    ├── mapped VSD run
    ├── custom residue run
    └── aligned custom residue run
```

The repository includes example tabular data so users can inspect the expected residue-map format and output format without needing to download a molecular-dynamics trajectory.

## Quick start

### List available mapped systems

```bash
python vsd_ztrack.py --list-map
```

### Run a mapped VSD analysis

Use this when the residue map already contains the channel, VSD label, hydrophobic constriction site, and gating-charge residues.

```bash
python vsd_ztrack.py \
  --psf <topology_file> \
  --traj <trajectory_file> \
  --gene <gene_name> \
  --vsd all \
  --time-per-frame <time_between_saved_frames> \
  --time-unit ns \
  --y-limits -20 20 \
  --outdir <output_directory>
```

Example with a mapped sodium-channel entry:

```bash
python vsd_ztrack.py \
  --psf <topology_file> \
  --traj <trajectory_file> \
  --gene SCN9A \
  --vsd all \
  --time-per-frame 0.24 \
  --time-unit ns \
  --y-limits -20 20 \
  --outdir results_mapped_vsd
```

This command uses the mapped hydrophobic constriction site and mapped gating-charge residues for each selected VSD.

## Custom residue mode

Use custom residue mode when the residues are not available in the map, when the topology uses modified numbering, or when you want to measure residues other than canonical gating charges.

```bash
python vsd_ztrack.py \
  --psf <topology_file> \
  --traj <trajectory_file> \
  --hcs-residue <reference_residue_number> \
  --target-residues <residue_1> <residue_2> <residue_3> \
  --target-labels <label_1> <label_2> <label_3> \
  --custom-vsd-label <domain_label> \
  --time-per-frame <time_between_saved_frames> \
  --time-unit ns \
  --outdir results_custom_residues
```

Example:

```bash
python vsd_ztrack.py \
  --psf <topology_file> \
  --traj <trajectory_file> \
  --hcs-residue 163 \
  --target-residues 214 217 220 223 \
  --target-labels R1 R2 R3 R4 \
  --custom-vsd-label VSD1 \
  --time-per-frame 0.24 \
  --time-unit ns \
  --outdir results_custom_residues
```

This calculates:

```text
coordinate(target residue) - coordinate(reference residue)
```

for each selected target residue at each trajectory frame.

## Use mapped HCS with custom target residues

This mode uses the residue map to define the hydrophobic constriction site while allowing the user to define the target residues manually.

```bash
python vsd_ztrack.py \
  --psf <topology_file> \
  --traj <trajectory_file> \
  --gene SCN9A \
  --vsd VSD1 \
  --target-residues 214 217 220 223 \
  --target-labels R1 R2 R3 R4 \
  --time-per-frame 0.24 \
  --time-unit ns \
  --outdir results_mapped_reference_custom_targets
```

## Analyze selected mapped gating charges

Use mapped labels when you want only part of the S4 gating-charge set.

```bash
python vsd_ztrack.py \
  --psf <topology_file> \
  --traj <trajectory_file> \
  --gene SCN9A \
  --vsd VSD1 \
  --target-r-labels R1 R2 R3 \
  --time-per-frame 0.24 \
  --time-unit ns \
  --outdir results_selected_mapped_targets
```

## Optional protein alignment

Use alignment when the protein drifts, translates, or rotates in the simulation box. The default alignment selection uses protein C-alpha atoms.

```bash
python vsd_ztrack.py \
  --psf <topology_file> \
  --traj <trajectory_file> \
  --gene SCN9A \
  --vsd VSD1 \
  --align-protein \
  --align-selection "protein and name CA" \
  --align-reference-frame 0 \
  --time-per-frame 0.24 \
  --time-unit ns \
  --outdir results_aligned
```

If the system contains more than one protein segment with overlapping residue numbers, restrict the analysis to the intended segment:

```bash
--extra-selection "segid <segment_id>"
```

or exclude an unwanted segment:

```bash
--exclude-segid <segment_id>
```

Explicit inclusion with `--extra-selection` is usually safer than exclusion because it selects the intended segment directly.

## Reference and target representation

### Default behavior

```text
Reference residue = C-alpha atom
Target residue = side-chain heavy-atom center of mass
Axis = z
```

### Target residue options

Use side-chain heavy-atom center of mass:

```bash
--target-reference sidechain_com
```

Use whole-residue center of mass:

```bash
--target-reference residue_com
```

Use target C-alpha atom:

```bash
--target-reference ca --target-atom-name CA
```

### Reference residue options

Use reference C-alpha atom:

```bash
--hcs-reference ca --hcs-atom-name CA
```

Use reference side-chain heavy-atom center of mass:

```bash
--hcs-reference sidechain_com
```

Use reference whole-residue center of mass:

```bash
--hcs-reference residue_com
```

## Coordinate axis

The default coordinate axis is `z`, which is usually the membrane-normal axis in channel simulations.

```bash
--axis z
```

Other axes can be selected if needed:

```bash
--axis x
--axis y
```

## Time scaling

If the trajectory stores physical time correctly, the script can use the trajectory time and apply `--time-scale`.

If the trajectory does not store the desired time directly, provide the time between saved frames:

```bash
--time-per-frame <time_between_saved_frames>
```

For example, if a trajectory represents 3000 ns and contains 12501 saved frames:

```text
time_between_saved_frames = 3000 / (12501 - 1)
```

Use the resulting value as:

```bash
--time-per-frame 0.24
```

## Frame subsetting

Analyze only part of the trajectory:

```bash
--start <first_frame> --stop <stop_before_this_frame> --step <stride>
```

Example:

```bash
--start 0 --stop 5000 --step 10
```

## Output files

The output directory contains:

```text
long-form trajectory table
summary table
selected-entry table
run metadata
plots for each analyzed VSD or custom domain
```

The long-form table contains one row per analyzed frame and target residue. Important columns include:

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

The summary table reports the mean value at the beginning and end of the trajectory:

```text
initial_mean_A
final_mean_A
delta_final_minus_initial_A
```

The averaging window is controlled with:

```bash
--avg-window <number_of_frames>
```

## Example data included in this repository

The `example_data` directory contains small tabular examples:

- a subset of the residue-map format,
- a custom target-residue table,
- an example long-form output table,
- an example summary output table.

These files are for documentation and format checking only. They are not MD trajectories.

## Troubleshooting

### Empty residue selection

If the script reports an empty selection, check whether the residue numbering in the topology matches the residue map. Common fixes are:

```bash
--resid-offset <integer_offset>
```

or:

```bash
--residue-selector resnum
```

### Duplicate reference atom selection

If the script reports more than one matching reference atom, the topology probably contains multiple protein segments with overlapping residue numbers. Identify the matching atoms and then rerun with a segment filter:

```bash
--extra-selection "segid <segment_id>"
```

### Output appears shifted by whole-protein motion

Rerun with alignment:

```bash
--align-protein --align-selection "protein and name CA"
```

## Minimal citation text

If this tool supports a publication, cite this publication
Elhanafy, E., Akbari Ahangar, A., Roth, R., Gamal El-Din, T. M., Bankston, J. R., & Li, J. (2025). The differential impacts of equivalent gating-charge mutations in voltage-gated sodium channels. Journal of General Physiology, 157(2), e202413669. https://doi.org/10.1085/jgp.202413669

## Repository description

Residue-map-driven MDAnalysis tool for tracking VSD gating-charge and custom residue displacement relative to a reference residue in ion-channel molecular-dynamics simulations.
