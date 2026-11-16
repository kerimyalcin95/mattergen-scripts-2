⚠️ The README.md was initially generated with AI assistance and has not yet been manually reviewed. Its contents should therefore be considered a draft and may contain inaccuracies, omissions, or outdated information that will be corrected during review.

# MatterGen-Scripts

A collection of Python scripts for benchmarking MatterGen-generated crystal structures against reference crystal-structure databases.

The project focuses on structural properties that can be derived directly from crystal structure files. The analyses cover basic structure statistics, crystallographic symmetry, local coordination, bond geometry, radial distribution functions (RDF), and crystal fingerprints.

The scripts are designed as independent command-line tools with shared utilities for:

- structure loading
- file discovery
- parallel processing
- result/error handling
- statistics
- plotting
- command-line argument parsing

---

## Project Overview

The benchmark pipeline is organized as follows:

| Script | Purpose |
| --- | --- |
| `00_download_database.py` | Download structures from supported databases |
| `01_basic_statistics.py` | Basic structural statistics |
| `02_symmetry_analysis.py` | Crystallographic symmetry |
| `03_coordination_analysis.py` | CrystalNN/VoronoiNN coordination |
| `04_bond_analysis.py` | Bond lengths, bond angles, nearest neighbors |
| `05_rdf_analysis.py` | Radial distribution functions |
| `06_fingerprint_analysis.py` | Crystal fingerprints and structural similarity |

Each analysis can be executed independently.

---


# Table of Contents

- [MatterGen-Scripts](#mattergen-scripts)
- [Quickstart: `launcher.py`](#quickstart-launcherpy)
- [Requirements](#requirements)
- [Installation](#installation)
- [Project Structure](#project-structure)
- [Database Downloading](#database-downloading)
- [Dataset Organization](#dataset-organization)
- [Common Analysis CLI](#common-analysis-cli)
- [General Analysis Workflow](#general-analysis-workflow)
- [Output Organization](#output-organization)
- [01_basic_statistics.py](#01_basic_statisticspy)
- [02_symmetry_analysis.py](#02_symmetry_analysispy)
- [03_coordination_analysis.py](#03_coordination_analysispy)
- [04_bond_analysis.py](#04_bond_analysispy)
- [05_rdf_analysis.py](#05_rdf_analysispy)
- [06_fingerprint_analysis.py](#06_fingerprint_analysispy)
- [Structural Benchmark Methodology](#structural-benchmark-methodology)
- [RDF vs. Fingerprints](#rdf-vs-fingerprints)
- [Limitations](#limitations)
- [Troubleshooting](#troubleshooting)
- [License](#license)

# Quickstart: `launcher.py`

The recommended way to run the project is through the interactive launcher.

Start it from the project directory:

```bash
python launcher.py
```

The launcher guides you through the complete workflow:

1. Select an action:
   - Download datasets
   - Analyze datasets
   - Download + analyze
   - Exit
2. Select the reference databases when downloading.
3. Enter one or more chemical systems, for example `Al-O`, `Fe-O`, or `Y-Al-O`.
4. Select the datasets to analyze.
5. Select one or more analyses.
6. Configure the worker count.
7. Configure analysis-specific parameters when required.
8. Select the data and results directories.
9. Review the complete configuration.
10. Confirm the workflow.

The launcher calls each numbered script separately with the appropriate arguments. It does **not** search structure directories itself and does **not** recursively discover files. File discovery and empty-directory handling remain the responsibility of the individual analysis scripts.

### Launcher workflow

```text
launcher.py
    │
    ├── Download datasets
    │       └── 00_download_database.py
    │
    └── Analyze datasets
            ├── 01_basic_statistics.py
            ├── 02_symmetry_analysis.py
            ├── 03_coordination_analysis.py
            ├── 04_bond_analysis.py
            ├── 05_rdf_analysis.py
            └── 06_fingerprint_analysis.py
```

### Analysis parameters

| Analysis | Parameters |
| --- | --- |
| Basic statistics | `--input`, `--output`, `--workers` |
| Symmetry | `--input`, `--output`, `--workers`, `--symprec` |
| Coordination | `--input`, `--output`, `--workers` |
| Bond | `--input`, `--output`, `--workers` |
| RDF | `--input`, `--output`, `--workers`, `--rmax`, `--bin-size` |
| Fingerprint | `--input`, `--output`, `--workers` |

The launcher intentionally does not pass `--recursive`. Each analysis receives the exact dataset directory selected by the launcher.

### Example selection

```text
Action:
    Download + analyze

Chemical systems:
    Al-O
    Fe-O

Download databases:
    Materials Project
    OQMD
    JARVIS

Analysis datasets:
    MatterGen
    Materials Project
    OQMD
    JARVIS

Analyses:
    Basic statistics
    Symmetry
    Coordination
    Bond
    RDF
    Fingerprint
```

After confirming the displayed configuration, the launcher executes the selected commands.


## Features

### Database Downloading

`00_download_database.py` provides a common entry point for supported database downloaders.

Currently registered databases are:

- Materials Project (`mp`)
- OQMD (`oqmd`)
- JARVIS (`jarvis`)
- Alexandria (`alexandria`)
- Crystallography Open Database (`cod`)

The downloader creates the database/chemical-system output directory through the selected downloader and then delegates the actual download to that implementation.

### Crystal Structure Analysis

The analysis scripts support:

- CIF
- XYZ
- EXTXYZ

Depending on the analysis, periodic information may be required. CIF and EXTXYZ are therefore preferred for periodic crystal-structure analyses.


### Output

Analyses write machine-readable CSV data and figures in PDF/PNG format where implemented.

Invalid structures are represented in the result tables instead of stopping the complete dataset analysis.

### Parallel Processing

Structure-level analyses use the shared `common.parallel.parallel_map` helper and expose a common `--workers` argument.

---

# Requirements

## Python

Python 3.11 or newer is recommended.

The project uses the following packages:

| Package | Main use |
| --- | --- |
| `numpy` | Numerical operations |
| `scipy` | Scientific/statistical operations |
| `pandas` | Tabular data |
| `matplotlib` | Plotting |
| `seaborn` | Statistical plotting |
| `pymatgen` | Crystal structures and symmetry |
| `pymatviz` | RDF calculation |
| `matminer` | Crystal fingerprints |
| `ase` | Structure-file support |
| `spglib` | Crystallographic symmetry backend |
| `scikit-learn` | Cosine similarity and PCA |
| `joblib` | Parallel processing backend |
| `tqdm` | Progress reporting |
| `python-dotenv` | `.env` configuration for database downloads |
| `requests` | HTTP requests |
| `mp-api` | Materials Project API |
| `questionary` | Interactive launcher prompts |
| `rich` | Launcher console output and formatting |

Install the project dependencies with:

```bash
pip install -r requirements.txt
```

If `requirements.txt` is not available, install all project dependencies with:

```bash
pip install numpy scipy pandas matplotlib seaborn pymatgen pymatviz matminer ase spglib scikit-learn joblib tqdm python-dotenv requests mp-api questionary rich
```

This includes the dependencies required by the launcher, database downloaders, and analysis scripts.


# Installation

Clone the repository:

```bash
git clone <repository-url>
cd MatterGen-Scripts
```

Install dependencies:

```bash
pip install -r requirements.txt
```

For the Materials Project downloader, an API key can be supplied directly or through a `.env` file using `MP_API_KEY`.

Example `.env` entry:

```text
MP_API_KEY=your_api_key
```

Do not commit API keys or other credentials to version control.

---

# Project Structure

A typical project layout is:

```text
MatterGen-Scripts/
├── 00_download_database.py
├── 01_basic_statistics.py
├── 02_symmetry_analysis.py
├── 03_coordination_analysis.py
├── 04_bond_analysis.py
├── 05_rdf_analysis.py
├── 06_fingerprint_analysis.py
├── common/
│   ├── cli.py
│   ├── io.py
│   ├── parallel.py
│   ├── plotting.py
│   ├── results.py
│   └── statistics.py
├── downloaders/
│   ├── materials_project.py
│   ├── oqmd.py
│   ├── jarvis.py
│   ├── alexandria.py
│   └── cod.py
├── requirements.txt
├── README.md
├── data/
└── results/
```

The exact contents of `common/` and `downloaders/` may change as the project evolves.

---

# Database Downloading

## `00_download_database.py`

The downloader is the common entry point for database acquisition.

### Command

```bash
python 00_download_database.py \
    --database mp \
    --chemsys Al-O Fe-O Y-Al-O \
    --output data \
    --api-key YOUR_API_KEY
```

### Arguments

| Argument | Required | Description |
| --- | --- | --- |
| `--database` | Yes | Database key registered in `DOWNLOADERS` |
| `--chemsys` | Yes | One or more chemical systems |
| `--output` | Yes | Root output directory |
| `--api-key` | No | API key; defaults to `MP_API_KEY` |
| `--provider` | No | OPTIMADE provider, when supported by the downloader |

### Database Keys

| Key | Database |
| --- | --- |
| `mp` | Materials Project |
| `oqmd` | OQMD |
| `jarvis` | JARVIS |
| `alexandria` | Alexandria |
| `cod` | Crystallography Open Database |

The downloader delegates database-specific behavior to the corresponding class in `downloaders/`.

---

# Dataset Organization

A convenient organization is:

```text
data/
├── MatterGen/
│   ├── Al-O/
│   │   └── generated_crystals_cif/
│   ├── Fe-O/
│   │   └── generated_crystals_cif/
│   └── Y-Al-O/
│       └── generated_crystals_cif/
├── MaterialsProject/
│   ├── Al-O/
│   ├── Fe-O/
│   └── Y-Al-O/
├── OQMD/
├── JARVIS/
├── Alexandria/
└── COD/
```

The analysis scripts operate on the directory passed through `--input`.

Each chemical system can therefore be analyzed separately and compared afterward.

---

# Common Analysis CLI

The analysis scripts use the shared parser from `common.cli`.

The common interface is:

```bash
python <script>.py \
    --input <input_directory> \
    --output <output_directory> \
    --workers -1
```

The scripts also expose `--recursive` through the common parser.

### Common Arguments

| Argument | Description |
| --- | --- |
| `--input` | Input directory containing supported structure files |
| `--output` | Root directory for analysis results |
| `--workers` | Number of parallel workers |
| `--recursive` | Enable recursive input-file discovery |

For analyses that have additional parameters, those parameters are documented in the corresponding section below.

---

# General Analysis Workflow

The analysis scripts follow the common pattern:

```text
Input directory
       │
       ▼
Find supported structure files
       │
       ▼
Load structures
       │
       ▼
Analyze structures in parallel
       │
       ▼
Collect per-structure results
       │
       ▼
Print summary
       │
       ├──► CSV results
       │
       └──► PDF/PNG figures
```

The shared utilities keep file loading, parallel execution, statistics, plotting, and result/error handling consistent across analyses.

---

# Output Organization

Analyses that use `create_analysis_output_directory()` create a dedicated analysis directory below the requested output root.

Typical output:

```text
results/
└── <database>/
    └── <chemical-system>/
        ├── symmetry/
        ├── coordination/
        ├── bond/
        ├── rdf/
        └── fingerprint/
```

`01_basic_statistics.py` currently writes directly to the output directory passed through `--output`.

The exact files depend on the analysis and on whether valid structures were available.

---

# 01_basic_statistics.py

## Purpose

Computes basic structural statistics for each crystal structure.

The script uses:

- `pymatgen.Structure`
- `SpacegroupAnalyzer`
- shared structure loading
- shared parallel processing
- shared result/statistics helpers
- shared plotting utilities

## Computed Properties

For each structure:

- chemical formula
- reduced formula
- number of atoms
- number of unique elements
- lattice parameters `a`, `b`, `c`
- lattice angles `alpha`, `beta`, `gamma`
- cell volume
- density
- space-group symbol

Space-group detection uses `symprec=0.1`.

## CSV Output

The main analysis file is:

```text
basic_statistics.csv
```

The CSV contains the structure-level results and the selected numeric analysis columns.

## Figures

The script generates:

```text
num_atoms_histogram.pdf
num_atoms_histogram.png
volume_histogram.pdf
volume_histogram.png
density_histogram.pdf
density_histogram.png
lattice_parameters.pdf
lattice_parameters.png
spacegroup_distribution.pdf
spacegroup_distribution.png
composition_distribution.pdf
composition_distribution.png
```

## Example

```bash
python 01_basic_statistics.py \
    --input data/MaterialsProject/Al-O \
    --output results/MaterialsProject/Al-O \
    --workers -1
```

---

# 02_symmetry_analysis.py

## Purpose

Analyzes crystallographic symmetry for every structure.

## Additional Argument

`--symprec` controls the symmetry tolerance.

Default:

```text
--symprec 0.1
```

Example:

```bash
python 02_symmetry_analysis.py \
    --input data/MaterialsProject/Al-O \
    --output results/MaterialsProject/Al-O \
    --workers -1 \
    --symprec 0.1
```

## Computed Properties

- space-group number
- space-group symbol
- crystal system
- Bravais lattice
- point group
- Hall symbol
- number of symmetry operations
- symmetry tolerance used
- validity status

## Output

The analysis CSV is:

```text
symmetry/symmetry_analysis.csv
```

The script also produces summary statistics and figures for:

- crystal systems
- space groups
- space-group numbers
- point groups
- symmetry-operation counts

---

# 03_coordination_analysis.py

## Purpose

Compares local coordination environments using two independent neighbor-finding methods.

## Methods

### CrystalNN

`CrystalNN` is used with:

- `distance_cutoffs=None`
- `x_diff_weight=0`
- `porous_adjustment=False`

### VoronoiNN

`VoronoiNN` is used as a geometry-based comparison method.

## Computed Properties

For both methods:

- mean coordination
- median coordination
- standard deviation
- minimum coordination
- maximum coordination

Additional comparison metrics:

- mean absolute difference
- maximum absolute difference
- agreement fraction

The agreement fraction is the fraction of sites for which the two methods return identical coordination numbers, using the valid site-level results.

## Output

Main CSV:

```text
coordination/coordination_analysis.csv
```

Figures include:

```text
crystalnn_coordination_histogram
voronoi_coordination_histogram
agreement_fraction_histogram
coordination_difference_histogram
coordination_scatter
coordination_boxplot
composition_distribution
```

Each figure is written in the formats supported by the shared plotting utilities.

## Example

```bash
python 03_coordination_analysis.py \
    --input data/MaterialsProject/Al-O \
    --output results/MaterialsProject/Al-O \
    --workers -1
```

---

# 04_bond_analysis.py

## Purpose

Analyzes local bond geometry using CrystalNN-derived neighbors.

## Computed Properties

For each structure:

### Bond lengths

- mean
- median
- standard deviation
- minimum
- maximum
- number of unique bonds

### Nearest-neighbor distances

- mean
- median
- standard deviation
- minimum
- maximum

### Coordination

- mean
- median
- standard deviation
- minimum
- maximum

### Bond angles

- mean
- median
- standard deviation
- minimum
- maximum
- number of angles

Bond angles are calculated for pairs of neighbors around each central atom.

## Example

```bash
python 04_bond_analysis.py \
    --input data/MaterialsProject/Al-O \
    --output results/MaterialsProject/Al-O \
    --workers -1
```

## Output

Main CSV:

```text
bond/bond_analysis.csv
```

Figures include:

```text
bond_length_histogram
bond_angle_histogram
nearest_neighbor_histogram
coordination_histogram
bond_statistics_boxplot
```

The exact PDF/PNG filenames are handled by `common.plotting.save_figure()`.

---

# 05_rdf_analysis.py

## Purpose

Computes radial distribution functions using `pymatviz`.

RDF analysis describes atomic-distance distributions over a specified radial range and complements the local bond analysis.

## Parameters

### `--rmax`

Maximum RDF distance in Å.

Default:

```text
10.0
```

### `--bin-size`

RDF histogram bin width in Å.

Default:

```text
0.02
```

## Example

```bash
python 05_rdf_analysis.py \
    --input data/MaterialsProject/Al-O \
    --output results/MaterialsProject/Al-O \
    --workers -1 \
    --rmax 10.0 \
    --bin-size 0.02
```

## Computed Properties

The script extracts:

- first RDF peak position
- first RDF peak height
- coordination-shell radius
- RDF integral
- maximum RDF value

The full RDF curves are also retained.

The first detected peak is obtained with `scipy.signal.find_peaks()` using `prominence=0.05`.

## Output

Main CSV:

```text
rdf/rdf_analysis.csv
```

Full RDF curves:

```text
rdf/rdf_curves.csv
```

Figures include:

- average RDF
- RDF heatmap
- first peak position distribution
- first peak height distribution
- first coordination-shell distribution

---

# 06_fingerprint_analysis.py

## Purpose

Computes several matminer-based crystal fingerprints and uses them for structural similarity analysis.

## Fingerprints

### CrystalNN fingerprint

Created from:

```text
SiteStatsFingerprint(CrystalNNFingerprint.from_preset("ops"))
```

### OPSite fingerprint

Created from:

```text
SiteStatsFingerprint(OPSiteFingerprint())
```

### Voronoi fingerprint

Created from:

```text
SiteStatsFingerprint(VoronoiFingerprint())
```

Each structure therefore receives three structure-level fingerprint vectors.

## Fingerprint Norms

The script computes the Euclidean norm of each fingerprint:

- `crystalnn_norm`
- `opsite_norm`
- `voronoi_norm`

## Cosine Similarity

Pairwise cosine similarity is calculated for the CrystalNN fingerprint matrix.

The script additionally computes:

- nearest-neighbor similarity
- mean similarity

The self-similarity value is excluded when determining these metrics.

## PCA

The CrystalNN fingerprint matrix is projected to two dimensions with:

```text
PCA(n_components=2, random_state=42)
```

The resulting file contains:

- `file`
- `pc1`
- `pc2`

and is merged with:

- `formula`
- `reduced_formula`

## Output

Main CSV:

```text
fingerprint/fingerprint_analysis.csv
```

Fingerprint matrices:

```text
crystalnn_fingerprints.csv
opsite_fingerprints.csv
voronoi_fingerprints.csv
```

Similarity matrix:

```text
similarity_matrix.csv
```

PCA coordinates:

```text
pca_coordinates.csv
```

Figures:

```text
fingerprint_similarity_histogram
fingerprint_heatmap
fingerprint_pca
fingerprint_boxplot
```

## Example

```bash
python 06_fingerprint_analysis.py \
    --input data/MaterialsProject/Al-O \
    --output results/MaterialsProject/Al-O \
    --workers -1
```

---

# Structural Benchmark Methodology

The analyses should be interpreted together because they describe different levels of structural information.

## 1. Basic Structural Statistics

Basic statistics describe properties such as:

- cell dimensions
- volume
- density
- number of atoms
- space-group distribution

These are useful first-level checks for obvious structural differences between datasets.

## 2. Crystallographic Symmetry

Symmetry describes the global periodic arrangement.

Important descriptors include:

- crystal system
- space group
- point group
- Hall symbol
- Bravais lattice
- symmetry-operation count

## 3. Local Coordination

Coordination describes which atoms are considered local neighbors and how many neighbors each site has.

Using both CrystalNN and VoronoiNN provides a comparison between chemically informed and geometric neighbor definitions.

## 4. Bond Geometry

Bond analysis describes:

- local distances
- nearest-neighbor distances
- coordination
- angles between neighboring atoms

## 5. RDF

RDF analysis describes the distribution of pair distances over the selected radial range.

It is a statistical descriptor and should not be interpreted as a direct list of chemically bonded pairs.

## 6. Crystal Fingerprints

Fingerprints represent local environments as numerical vectors and enable structure-level statistical comparisons.

---

# RDF vs. Fingerprints

| RDF | Fingerprints |
| --- | --- |
| Pair-distance distribution | Local-environment descriptor |
| Explicit radial coordinate | High-dimensional vector |
| Captures distance statistics | Captures local geometric/chemical information |
| Useful for dataset-level comparison | Useful for structure-level similarity |
| Sensitive to selected radial range/binning | Sensitive to descriptor definition |

Using both methods gives complementary information rather than duplicate measurements.

---


# Limitations

The analyses operate on structural data. They do not directly determine electronic or thermodynamic properties.

The following require additional calculations or database metadata:

- formation energy
- energy above hull
- band gap
- density of states
- elastic constants
- magnetic properties
- phonon properties
- electronic charge density

The results also depend on the algorithms and numerical parameters used by the underlying libraries.

In particular:

- symmetry results depend on `symprec`
- RDF results depend on `rmax` and `bin-size`
- coordination depends on the neighbor-finding algorithm
- fingerprints depend on the selected matminer descriptors
- similarity values depend on the chosen fingerprint representation

---

# Troubleshooting

## CrystalNN Radius Warning

You may encounter warnings such as:

```text
CrystalNN: cannot locate an appropriate radius...
```

This indicates that an appropriate radius was not found for one or more elements.

The underlying CrystalNN implementation may fall back to another available radius source. The warning should therefore be interpreted as a limitation of the neighbor determination for the affected species, not automatically as a failed structure.

## XYZ Files

Plain XYZ files generally do not contain lattice vectors.

Consequently, analyses that depend on periodic crystal information may require CIF or EXTXYZ input instead.

## Invalid Structures

Structure-level exceptions are caught by the analysis workers and represented as failed results where the corresponding script supports the shared result helpers.

Check the generated analysis CSV and failure information when diagnosing problematic structures.

## Memory Usage

If a large dataset causes excessive memory consumption, reduce the worker count:

```bash
python 06_fingerprint_analysis.py \
    --input data/MatterGen/Al-O/generated_crystals_cif \
    --output results/MatterGen/Al-O \
    --workers 4
```

---


# License

This project is intended for academic research and benchmarking of crystal-structure generation models.

Refer to the repository license for the applicable usage and redistribution terms.
