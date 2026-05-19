import os
import glob
import argparse
import numpy as np
import mdtraj as md
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA, IncrementalPCA
from sklearn.cluster import MiniBatchKMeans
from scipy.stats import gaussian_kde
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

scale = 1.2
for key in [
    "font.size",
    "axes.labelsize",
    "axes.titlesize",
    "xtick.labelsize",
    "ytick.labelsize",
    "legend.fontsize",
    "figure.titlesize",
]:
    val = plt.rcParams.get(key)
    if isinstance(val, (int, float)):
        plt.rcParams[key] = val * scale

def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--base-folder","-b",required=True)
    parser.add_argument("--ref-pdb","-r",required=True)
    parser.add_argument("--out-folder","-o",required=True)
    parser.add_argument("--max-rmsd","-m",type=float,default=1.0)
    parser.add_argument("--samples-per-system","-s",type=int,default=5000)
    return parser.parse_args()

def prepare_trajectory(traj_files, ref, bb_idx):
    trajs = []
    times = []
    for f in traj_files:
        t = md.load(f)
        t = t.atom_slice(t.top.select("not water"))
        t.superpose(ref, atom_indices=bb_idx)
        trajs.append(t)
        times.extend(range(t.n_frames))
    traj = trajs[0].join(trajs[1:])
    traj_bb = traj.atom_slice(bb_idx)
    coords = traj_bb.xyz.reshape(traj_bb.n_frames, -1)
    return traj, traj_bb, np.array(times), coords

def plot_pca(X2, time_array, name, plot_dir, centroid=None, padding=0.05):
    x = X2[:, 0]
    y = X2[:, 1]

    xrange = x.max() - x.min()
    yrange = y.max() - y.min()
    pad_x = padding * xrange if xrange > 0 else 1e-6
    pad_y = padding * yrange if yrange > 0 else 1e-6

    xmin = x.min() - pad_x
    xmax = x.max() + pad_x
    ymin = y.min() - pad_y
    ymax = y.max() + pad_y

    nx = ny = 300
    xgrid = np.linspace(xmin, xmax, nx)
    ygrid = np.linspace(ymin, ymax, ny)
    xx, yy = np.meshgrid(xgrid, ygrid)

    kde = gaussian_kde(np.vstack([x, y]))
    density = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    point_density = kde(np.vstack([x, y]))
    levels = np.linspace(point_density.min(), density.max(), 15)

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(16, 6), dpi=300, constrained_layout=True,
        sharex=True, sharey=True
    )

    sc = ax1.scatter(
        x, y, c=time_array, cmap="viridis",
        s=10, alpha=0.6, linewidth=0, rasterized=True
    )

    ax1.contour(xx, yy, density, levels=levels, colors="black",
                linewidths=1.2, alpha=0.85)

    if centroid is not None:
        ax1.scatter(
            centroid[0], centroid[1],
            marker="X", s=120,
            color="black", edgecolor="white", linewidth=0.8
        )

    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    ax1.set_xlabel("PC1")
    ax1.set_ylabel("PC2")

    cbar1 = fig.colorbar(sc, ax=ax1, fraction=0.046, pad=0.04)
    cbar1.set_label("Frame", fontsize=14.4)

    density_masked = np.ma.masked_where(density < levels[0], density)

    im = ax2.imshow(
        density_masked,
        extent=[xmin, xmax, ymin, ymax],
        origin='lower',
        cmap="turbo_r",
        alpha=0.75,
        aspect="auto"
    )

    ax2.contour(xx, yy, density, levels=levels, colors="black",
                linewidths=1.2, alpha=0.85)

    if centroid is not None:
        ax2.scatter(
            centroid[0], centroid[1],
            marker="X", s=120,
            color="black", edgecolor="white", linewidth=0.8
        )

    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)
    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")

    cbar2 = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
    cbar2.set_label("Density (blue = most populated)", fontsize=14.4)

    outpath = os.path.join(plot_dir, f"{name}_combined.png")
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)

def main():
    args = parse_args()
    base_folder = args.base_folder
    ref_pdb = args.ref_pdb
    out_folder = args.out_folder
    MAX_RMSD = args.max_rmsd
    N_SAMPLES = args.samples_per_system

    os.makedirs(out_folder, exist_ok=True)

    ref = md.load(ref_pdb)
    bb_idx = ref.topology.select("backbone")

    folders = sorted(f for f in glob.glob(os.path.join(base_folder,"MD*")) if os.path.isdir(f))

    pca2 = IncrementalPCA(n_components=2)
    pca20 = IncrementalPCA(n_components=20)

    for folder in folders:
        traj_files = sorted(glob.glob(os.path.join(folder,"*.h5")))
        if not traj_files:
            continue
        _, traj_bb, _, coords = prepare_trajectory(traj_files, ref, bb_idx)
        sample_idx = np.linspace(0, coords.shape[0]-1, N_SAMPLES).astype(int)
        batch = coords[sample_idx]
        pca2.partial_fit(batch)
        pca20.partial_fit(batch)

    x_all = []
    y_all = []
    for folder in folders:
        traj_files = sorted(glob.glob(os.path.join(folder,"*.h5")))
        if not traj_files:
            continue
        _, traj_bb, _, coords = prepare_trajectory(traj_files, ref, bb_idx)
        sample_idx = np.linspace(0, coords.shape[0]-1, N_SAMPLES).astype(int)
        batch = coords[sample_idx]
        X2_global_batch = pca2.transform(batch)
        x_all.append(X2_global_batch[:,0])
        y_all.append(X2_global_batch[:,1])

    x_all = np.concatenate(x_all)
    y_all = np.concatenate(y_all)

    pad_x = 0.05 * (x_all.max() - x_all.min())
    pad_y = 0.05 * (y_all.max() - y_all.min())

    GLOBAL_XMIN = x_all.min() - pad_x
    GLOBAL_XMAX = x_all.max() + pad_x
    GLOBAL_YMIN = y_all.min() - pad_y
    GLOBAL_YMAX = y_all.max() + pad_y

    results = []

    for folder in folders:
        name = os.path.basename(folder)
        print("Processing:", name)
        traj_files = sorted(glob.glob(os.path.join(folder,"*.h5")))
        if not traj_files:
            continue

        plot_dir = os.path.join(out_folder, f"{name}_plots")
        os.makedirs(plot_dir, exist_ok=True)

        traj, traj_bb, time_array, coords = prepare_trajectory(traj_files, ref, bb_idx)

        X2 = pca2.transform(coords)
        X = pca20.transform(coords)

        km1 = MiniBatchKMeans(n_clusters=1, batch_size=5000, max_iter=100, random_state=42).fit(X)
        km2 = MiniBatchKMeans(n_clusters=2, batch_size=5000, max_iter=100, random_state=42).fit(X)

        center_frame = np.argmin(np.linalg.norm(X - km1.cluster_centers_[0], axis=1))
        traj[center_frame].save(os.path.join(out_folder, f"{name}.pdb"))

        centroid_2d = X2[center_frame]
        plot_pca(X2, time_array, name, plot_dir, centroid_2d)

        reference = traj_bb[center_frame]
        max_rmsd = md.rmsd(traj_bb, reference, 0).max()
        k1_ok = "Yes" if max_rmsd <= MAX_RMSD else "No"

        results.append({
            "System": name,
            "Max_RMSD_to_k1_center_Å": round(max_rmsd,3),
            "k=1_acceptable": k1_ok,
        })

    df = pd.DataFrame(results).sort_values("System")
    df.to_excel(os.path.join(out_folder,"clustering_summary.xlsx"), index=False)

    with open(os.path.join(out_folder,"clustering_summary_readme.txt"),"w") as f:
        f.write(
            f"k=1_acceptable = max RMSD ≤ {MAX_RMSD} Å\n"
            "Global PCA used for all systems.\n"
        )

if __name__ == "__main__":
    main()
