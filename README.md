# MDcentro - MD Trajectory Clustering & PCA Pipeline

This repository provides a Python workflow for analyzing **short molecular dynamics (MD) simulations**.  
It performs PCA dimensionality reduction, generates high-quality time-colored PCA plots, clusters conformations using MiniBatchKMeans, and extracts representative centroid structures.

## Features

- Processes short MD trajectories from `MD*` subdirectories  
- Removes water and aligns each frame to a reference backbone  
- Computes PCA (2D for visualization, 20D for clustering)  
- Generates clean, borderless PCA scatter plots  
- Performs fast k-means clustering (k=1 and k=2)  
- Saves the k=1 centroid structure as a PDB file  
- Calculates RMSD to assess clustering quality  
- Produces a summary Excel sheet for all systems  

## Conda Environment Setup

To ensure a stable and reproducible environment for running the MD PCA & clustering pipeline, create and activate the recommended Conda environment:

### Create environment with a stable Python version and required scientific libraries

```bash
conda create -y -n mdpca python=3.10 numpy scipy scikit-learn pandas matplotlib mdtraj h5py pip -c conda-forge
```
### Activate the environment

```bash
conda activate mdpca
```

## Usage

```bash
python3 process_md.py   --base-folder /path/to/MD_screening/   --ref-pdb /path/to/reference.pdb   --out-folder /path/to/output   --max-rmsd 2.0
```

## Input Arguments

- **`--base-folder`** — Directory containing `MD*` folders with `.h5` trajectories  
- **`--ref-pdb`** — Reference structure used for backbone alignment  
- **`--out-folder`** — Output directory for plots, centroid structures, and summary files  
- **`--max-rmsd`** — RMSD threshold (Å) for accepting k=1 clustering  

## Output Files

- High-resolution PCA plots  
- Representative centroid PDB structures  
- `clustering_summary.xlsx` with inertia, RMSD, and k=1 acceptance results  
- Additional text README explaining RMSD acceptance criteria  

## Example

MD trajectories (20 replicates) and the expected output files are provided for the enzyme Kemp HG3.R5. The variant has been described in our Nature Chemical Biology Paper (https://www.nature.com/articles/s41589-024-01712-3) and the MD simulations were conducted based on the high-throughput protocol described by Wang et al 2023 (https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00002).




---
