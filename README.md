# MDcentro -- MD Trajectory Clustering & PCA Pipeline

MDcentro provides a streamlined workflow for analyzing **short molecular
dynamics (MD) simulations**.\
It performs PCA dimensionality reduction, time‑colored PCA
visualization, clustering using MiniBatchKMeans, and extraction of
representative centroid structures.

## Features

-   Automatically discovers and processes `MD*` subdirectories\
-   Removes water and aligns trajectories to a reference backbone\
-   Computes PCA:
    -   **2D PCA** for visualization\
    -   **20D PCA** for clustering\
-   Produces high‑resolution PCA plots with KDE density overlays\
-   Performs fast MiniBatchKMeans clustering (`k=1` and `k=2`)\
-   Extracts the **k=1 centroid frame** and saves it as a PDB\
-   Calculates RMSD to evaluate centroid representativeness\
-   Generates:
    -   Per‑system PCA plots\
    -   Per‑system centroid PDBs\
    -   A global `clustering_summary.xlsx` report\
    -   A text README explaining RMSD acceptance rules

## Conda Environment Setup

### Create the environment

``` bash
conda env create -f mdcentro.yml
```

### Activate the environment

``` bash
conda activate mdcentro
```

## Usage

``` bash
python3 MDcentro.py     --base-folder /path/to/MD_screening/     --ref-pdb /path/to/reference.pdb     --out-folder /path/to/output     --max-rmsd 2.0
```

## Input Arguments

  -----------------------------------------------------------------------
  Argument                      Description
  ----------------------------- -----------------------------------------
  `--base-folder`               Folder containing `MD*` subdirectories
                                with `.h5` trajectories

  `--ref-pdb`                   Reference PDB used for backbone alignment

  `--out-folder`                Directory where plots, PDBs, and summary
                                files are written

  `--max-rmsd`                  RMSD threshold (Å) determining whether
                                k=1 clustering is acceptable

  `--samples-per-system`        Number of frames to sample for global PCA
                                (default: 5000)
  -----------------------------------------------------------------------

## Output Files

-   **PCA Plots:**
    -   Combined scatter + KDE density maps\
    -   Scaled-up font sizes for improved readability\
-   **Centroid PDB files** for each system\
-   **clustering_summary.xlsx** containing:
    -   Max RMSD to centroid\
    -   k=1 acceptability\
-   **clustering_summary_readme.txt** documenting clustering rules

## Example Data

Example MD trajectories and output files are available for the Shuffle library (N=199) between enzyme
**Kemp HG3.R5** and **Kemp HG3.17**, adescribed in:

-   **Nature Chemical Biology (2024)** --
    https://www.nature.com/articles/s41589-024-01712-3\
-   MD workflow based on **Wang et al., 2023** --
    https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00002
    
``` bash
python ./MDcentro.py --base-folder MD_trajectories --ref-pdb HG3_H2O.pdb --out-folder output_centroids
```
