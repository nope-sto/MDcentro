#!/usr/bin/env python3
import os
import glob
import argparse
import numpy as np
import mdtraj as md
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans


def parse_args():
    parser = argparse.ArgumentParser(
        description="Process MD trajectories, generate PCA plots, and compute k-means centroids."
    )

    parser.add_argument(
        "--base-folder", "-b",
        required=True,
        help="Base folder containing MD* subdirectories."
    )

    parser.add_argument(
        "--ref-pdb", "-r",
        required=True,
        help="Reference PDB file."
    )

    parser.add_argument(
        "--out-folder", "-o",
        required=True,
        help="Output folder for plots, centroids, and summary."
    )

    parser.add_argument(
        "--max-rmsd", "-m",
        type=float,
        default=2.0,
        help="Max RMSD threshold (Å) for accepting k=1 clustering."
    )

    return parser.parse_args()


def main():
    args = parse_args()

    base_folder = args.base_folder
    ref_pdb = args.ref_pdb
    out_folder = args.out_folder
    MAX_RMSD_FOR_K1 = args.max_rmsd

    os.makedirs(out_folder, exist_ok=True)

    ref = md.load(ref_pdb)
    bb_idx = ref.topology.select("backbone")

    folders = sorted(
        f for f in glob.glob(os.path.join(base_folder, "MD*"))
        if os.path.isdir(f)
    )

    results = []

    for folder in folders:
        name = os.path.basename(folder)
        print(f"Processing → {name}")

        traj_files = sorted(glob.glob(os.path.join(folder, "*.h5")))
        if not traj_files:
            continue

        plot_dir = os.path.join(out_folder, f"{name}_plots")
        os.makedirs(plot_dir, exist_ok=True)

        # Load trajectories
        trajs = []
        time_labels = []

        for f in traj_files:
            t = md.load(f)
            t = t.atom_slice(t.top.select("not water"))
            t.superpose(ref, atom_indices=bb_idx)
            trajs.append(t)
            time_labels.extend(range(t.n_frames))

        traj = trajs[0].join(trajs[1:])
        traj_bb = traj.atom_slice(bb_idx)
        time_array = np.array(time_labels)

        coords = traj_bb.xyz.reshape(traj_bb.n_frames, -1)
        X2 = PCA(n_components=2).fit_transform(coords)
        X = PCA(n_components=20).fit_transform(coords)

        # High-quality borderless plot
        fig, ax = plt.subplots(figsize=(7, 6), dpi=300)
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')

        sc = ax.scatter(
            X2[:, 0], X2[:, 1],
            c=time_array,
            s=12,
            cmap='viridis',
            alpha=0.5,
            linewidth=0,
            rasterized=True
        )

        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        cbar = plt.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label("Frame in original replicate", fontsize=12)
        cbar.ax.tick_params(labelsize=10)

        ax.set_title(name, fontsize=16, pad=20)
        plt.tight_layout()

        plt.savefig(
            os.path.join(plot_dir, f"{name}_time_colored.png"),
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.1,
            transparent=False
        )
        plt.close()

        # Fast clustering
        km1 = MiniBatchKMeans(n_clusters=1, batch_size=5000, max_iter=100, random_state=42).fit(X)
        km2 = MiniBatchKMeans(n_clusters=2, batch_size=5000, max_iter=100, random_state=42).fit(X)

        center_frame = np.argmin(np.linalg.norm(X - km1.cluster_centers_[0], axis=1))
        traj[center_frame].save(os.path.join(out_folder, f"{name}.pdb"))

        reference = traj_bb[center_frame]
        max_rmsd = md.rmsd(traj_bb, reference, 0).max()
        k1_ok = "Yes" if max_rmsd <= MAX_RMSD_FOR_K1 else "No"

        results.append({
            "System": name,
            "Inertia_k1": round(km1.inertia_, 3),
            "Inertia_k2": round(km2.inertia_, 3),
            "Max_RMSD_to_k1_center_Å": round(max_rmsd, 3),
            "k=1_acceptable": k1_ok
        })

    # Save summary
    df = pd.DataFrame(results).sort_values("System")
    out_excel = os.path.join(out_folder, "clustering_summary.xlsx")
    df.to_excel(out_excel, index=False)

    readme = (
        f"k=1_acceptable = 'Yes' if max backbone RMSD to k=1 center ≤ "
        f"{MAX_RMSD_FOR_K1} Å, else 'No'"
    )

    with open(os.path.join(out_folder, "clustering_summary_readme.txt"), "w") as f:
        f.write(readme + "\n")

    print(f"\nAll done. Summary → {out_excel}")


if __name__ == "__main__":
    main()
