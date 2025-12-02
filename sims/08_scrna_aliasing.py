#!/usr/bin/env python3
"""
Topological Aliasing in scRNA-seq: Real-World Application
==========================================================

Applies the Paper 2 falsifiability metrics to real single-cell data.

THE QUESTION:
When we project high-dimensional gene expression to UMAP (D_obs = 2),
how much topological aliasing do we introduce?

METRICS:
1. D_sys: Intrinsic dimensionality via participation ratio
2. D_obs: UMAP embedding dimension (2)
3. Aliasing rate: How often do UMAP neighbors differ from high-D neighbors?
4. Coverage: What fraction of high-D space is actually sampled?

Dataset: GSE120575 (Sade-Feldman melanoma, ~16k cells, ~55k genes)

Papers: "The Geometry of Biological Shadows" (Paper 2) & "The Limits of Falsifiability" (Paper 1)
Author: Ian Todd
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.manifold import TSNE
import os
import sys

# Add parent directory for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../23_immune_cooperation'))

os.makedirs('../figures', exist_ok=True)
os.makedirs('figures', exist_ok=True)

# Set random seed
np.random.seed(42)


def participation_ratio(data, n_components=100):
    """
    Compute effective dimensionality using PCA Participation Ratio.

    PR = (sum(lambda_i))^2 / sum(lambda_i^2)

    This is D_sys - the intrinsic dimensionality of the data manifold.
    """
    # Subsample for speed if too large
    if data.shape[0] > 3000:
        idx = np.random.choice(data.shape[0], 3000, replace=False)
        data = data[idx]

    # Center the data
    data_centered = data - np.mean(data, axis=0)

    # PCA
    n_comp = min(n_components, min(data.shape) - 1)
    pca = PCA(n_components=n_comp)
    pca.fit(data_centered)

    eigenvalues = pca.explained_variance_
    eigenvalues = eigenvalues[eigenvalues > 1e-10]

    pr = np.sum(eigenvalues)**2 / np.sum(eigenvalues**2)
    return pr, eigenvalues


def compute_topological_aliasing(data_high_d, data_low_d, k=10):
    """
    Measure topological aliasing: how often do low-D neighbors
    differ from high-D neighbors?

    For each point, we find k nearest neighbors in both spaces
    and measure the Jaccard distance between neighbor sets.

    Aliasing rate = 1 - (average Jaccard similarity)

    If UMAP perfectly preserved topology: aliasing = 0
    If UMAP completely scrambles neighbors: aliasing → 1
    """
    n_samples = min(data_high_d.shape[0], 5000)  # Subsample for speed
    idx = np.random.choice(data_high_d.shape[0], n_samples, replace=False)

    high_d = data_high_d[idx]
    low_d = data_low_d[idx]

    # Find k-nearest neighbors in each space
    nn_high = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree')
    nn_low = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree')

    nn_high.fit(high_d)
    nn_low.fit(low_d)

    # Get neighbor indices (excluding self)
    _, neighbors_high = nn_high.kneighbors(high_d)
    _, neighbors_low = nn_low.kneighbors(low_d)

    neighbors_high = neighbors_high[:, 1:]  # Exclude self
    neighbors_low = neighbors_low[:, 1:]

    # Compute Jaccard similarity for each point
    jaccard_sims = []
    for i in range(n_samples):
        set_high = set(neighbors_high[i])
        set_low = set(neighbors_low[i])
        intersection = len(set_high & set_low)
        union = len(set_high | set_low)
        jaccard_sims.append(intersection / union)

    aliasing_rate = 1 - np.mean(jaccard_sims)
    return aliasing_rate, np.array(jaccard_sims)


def compute_coverage(data, n_bins=3):
    """
    Measure coverage of high-dimensional space.

    Same metric as 07_sample_complexity.py but applied to real data.
    Uses PCA-reduced data (first 20 components) to make it tractable.
    """
    # Reduce to manageable dimensions for coverage computation
    n_dim = min(20, data.shape[1])
    pca = PCA(n_components=n_dim)
    data_reduced = pca.fit_transform(data)

    # Normalize to [0, 1] for binning
    data_norm = (data_reduced - data_reduced.min(axis=0)) / (data_reduced.max(axis=0) - data_reduced.min(axis=0) + 1e-10)

    # Discretize into bins
    binned = np.floor(data_norm * n_bins).astype(int)
    binned = np.clip(binned, 0, n_bins - 1)

    # Count unique cells occupied
    cells = set(tuple(b) for b in binned)
    n_occupied = len(cells)

    # Total possible cells
    n_total_cells = n_bins ** n_dim

    coverage = n_occupied / n_total_cells
    return coverage, n_occupied, n_total_cells, n_dim


def generate_synthetic_baseline(n_samples, n_genes, n_clusters=5):
    """Generate synthetic scRNA-seq-like data for comparison."""
    # Simulate cluster structure
    data = []
    for _ in range(n_samples):
        cluster = np.random.randint(n_clusters)
        center = np.zeros(n_genes)
        center[cluster * (n_genes // n_clusters):(cluster + 1) * (n_genes // n_clusters)] = 5
        noise = np.random.exponential(0.5, n_genes)
        data.append(center + noise)
    return np.array(data)


def run_simulation():
    """Main simulation with real or synthetic data."""

    print("=" * 60)
    print("TOPOLOGICAL ALIASING IN scRNA-seq DATA")
    print("Paper 2: The Geometry of Biological Shadows")
    print("=" * 60)

    # Try to load real data
    try:
        from fast_loader import load_sade_feldman_fast
        print("\nLoading Sade-Feldman melanoma data...")
        # Try multiple paths
        data_dir = "../../23_immune_cooperation/sade_feldman_data"
        data_resp, data_non, meta = load_sade_feldman_fast(data_dir)

        if data_resp is not None:
            # Combine for full analysis
            data = np.vstack([data_resp, data_non])
            labels = np.array(['Responder'] * len(data_resp) + ['Non-responder'] * len(data_non))
            data_source = "GSE120575 (Sade-Feldman)"
            print(f"Loaded {data.shape[0]} cells x {data.shape[1]} genes")
        else:
            raise ValueError("Data loading returned None")

    except Exception as e:
        print(f"Could not load real data: {e}")
        print("Using synthetic data for demonstration...")
        data = generate_synthetic_baseline(5000, 1000)
        labels = np.array(['Cluster'] * len(data))
        data_source = "Synthetic"

    # Filter zero-variance genes
    print("\nFiltering low-variance genes...")
    var = np.var(data, axis=0)
    keep = var > 0.01
    data_filt = data[:, keep]
    print(f"Kept {np.sum(keep)} / {len(keep)} genes")

    # =========================================================================
    # 1. COMPUTE D_sys (Intrinsic Dimensionality)
    # =========================================================================
    print("\n" + "-" * 40)
    print("[1] INTRINSIC DIMENSIONALITY (D_sys)")
    print("-" * 40)

    d_sys, eigenvalues = participation_ratio(data_filt)
    print(f"  Participation Ratio: D_sys = {d_sys:.1f}")
    print(f"  (This is the effective dimensionality of the gene expression manifold)")

    # =========================================================================
    # 2. COMPUTE LOW-D EMBEDDING (D_obs = 2)
    # =========================================================================
    print("\n" + "-" * 40)
    print("[2] LOW-DIMENSIONAL EMBEDDING (D_obs = 2)")
    print("-" * 40)

    # First reduce with PCA for speed
    print("  PCA reduction to 50 components...")
    pca = PCA(n_components=min(50, data_filt.shape[1] - 1))
    data_pca = pca.fit_transform(data_filt)

    # Then t-SNE (faster than UMAP without extra dependency)
    print("  t-SNE embedding to 2D...")
    n_subsample = min(5000, len(data_pca))
    idx_sub = np.random.choice(len(data_pca), n_subsample, replace=False)

    tsne = TSNE(n_components=2, perplexity=30, random_state=42, n_iter=1000)
    data_2d = tsne.fit_transform(data_pca[idx_sub])

    print(f"  Embedded {n_subsample} cells into 2D")
    print(f"  D_obs = 2 (the dimension of the shadow)")

    # =========================================================================
    # 3. COMPUTE TOPOLOGICAL ALIASING
    # =========================================================================
    print("\n" + "-" * 40)
    print("[3] TOPOLOGICAL ALIASING")
    print("-" * 40)

    aliasing_rate, jaccard_sims = compute_topological_aliasing(
        data_pca[idx_sub], data_2d, k=10
    )

    print(f"  Aliasing rate (k=10 neighbors): {aliasing_rate:.1%}")
    print(f"  (Fraction of high-D neighbors lost in 2D projection)")
    print(f"  Mean Jaccard similarity: {np.mean(jaccard_sims):.2f}")

    # =========================================================================
    # 4. COMPUTE COVERAGE
    # =========================================================================
    print("\n" + "-" * 40)
    print("[4] COVERAGE OF HIGH-D SPACE")
    print("-" * 40)

    coverage, n_occ, n_total, n_dim = compute_coverage(data_filt, n_bins=3)

    print(f"  Using first {n_dim} PCs with 3 bins per dimension")
    print(f"  Total possible cells: 3^{n_dim} = {n_total:,}")
    print(f"  Cells occupied: {n_occ:,}")
    print(f"  Coverage: {coverage:.2e} ({coverage*100:.4f}%)")

    # =========================================================================
    # PLOTTING
    # =========================================================================
    print("\n" + "-" * 40)
    print("Generating figures...")
    print("-" * 40)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel A: Eigenvalue spectrum
    ax = axes[0, 0]
    ax.semilogy(eigenvalues[:50], 'o-', color='#2A9D8F', markersize=4)
    ax.axhline(y=eigenvalues[int(d_sys)], color='red', linestyle='--',
               label=f'D_sys = {d_sys:.0f}')
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Explained Variance')
    ax.set_title(f'A. Eigenvalue Spectrum\n(D_sys = {d_sys:.1f})', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    # Panel B: t-SNE embedding with "Hairball of Truth"
    ax = axes[0, 1]
    ax.scatter(data_2d[:, 0], data_2d[:, 1], c='#457B9D', s=1, alpha=0.3, zorder=1)

    # Draw the "Hairball of Truth" - lines from random cells to their TRUE high-D neighbors
    np.random.seed(42)
    n_demo_cells = 15  # Show a few example cells
    demo_idx = np.random.choice(len(data_2d), n_demo_cells, replace=False)

    # Get high-D neighbors for demo cells
    nn_high = NearestNeighbors(n_neighbors=6, algorithm='ball_tree')
    nn_high.fit(data_pca[idx_sub])
    _, neighbors_high = nn_high.kneighbors(data_pca[idx_sub][demo_idx])

    # Draw lines to true neighbors (they'll cross the "empty" space)
    for i, cell_idx in enumerate(demo_idx):
        cell_pos = data_2d[cell_idx]
        for neighbor_idx in neighbors_high[i, 1:]:  # Skip self
            neighbor_pos = data_2d[neighbor_idx]
            ax.plot([cell_pos[0], neighbor_pos[0]], [cell_pos[1], neighbor_pos[1]],
                   'r-', alpha=0.4, linewidth=0.5, zorder=2)
        ax.scatter(cell_pos[0], cell_pos[1], c='red', s=20, zorder=3, edgecolor='black')

    ax.set_xlabel('t-SNE 1')
    ax.set_ylabel('t-SNE 2')
    ax.set_title(f'B. The Shadow Lie\n(red lines = TRUE neighbors in high-D)', fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])

    # Panel C: Jaccard similarity distribution
    ax = axes[1, 0]
    ax.hist(jaccard_sims, bins=30, color='#E63946', alpha=0.7, edgecolor='black')
    ax.axvline(x=np.mean(jaccard_sims), color='black', linestyle='--', linewidth=2,
               label=f'Mean = {np.mean(jaccard_sims):.2f}')
    ax.set_xlabel('Jaccard Similarity (high-D vs 2D neighbors)')
    ax.set_ylabel('Count')
    ax.set_title(f'C. Neighbor Preservation\n(Aliasing = {aliasing_rate:.1%})', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)

    # Panel D: Summary
    ax = axes[1, 1]
    ax.axis('off')

    summary_text = f"""
    TOPOLOGICAL ALIASING IN scRNA-seq
    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    Dataset: {data_source}
    Cells: {data.shape[0]:,}
    Genes: {data.shape[1]:,} → {np.sum(keep):,} (filtered)

    SYSTEM (D_sys):
      Participation Ratio = {d_sys:.1f}
      → The gene expression manifold has
        ~{d_sys:.0f} effective dimensions

    SHADOW (D_obs):
      t-SNE dimension = 2
      → We project to a 2D visualization

    ALIASING:
      Rate = {aliasing_rate:.1%}
      → {aliasing_rate:.0%} of high-D neighbors are
        "lost" in the 2D projection

    COVERAGE:
      {coverage:.2e} of 3^{n_dim} cells occupied
      → The data covers a vanishing fraction
        of the high-D space

    ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    IMPLICATION:

    When D_sys >> D_obs, the 2D "shadow"
    introduces systematic distortions.
    Clusters that appear close may be
    far apart in reality, and vice versa.

    This is topological aliasing applied
    to real biological data.
    """

    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
            fontsize=9, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='#f5f5f5', alpha=0.9))
    ax.set_title('D. The Epistemological Point', fontweight='bold')

    plt.tight_layout()

    # Save
    for path in ['figures/fig_scrna_aliasing.pdf', '../figures/fig_scrna_aliasing.pdf']:
        try:
            plt.savefig(path, dpi=150, bbox_inches='tight')
            print(f"Saved: {path}")
            break
        except Exception as e:
            continue

    plt.show()

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("\n" + "=" * 60)
    print("SUMMARY: FALSIFIABILITY METRICS ON REAL DATA")
    print("=" * 60)
    print(f"""
    D_sys (intrinsic dimensionality):  {d_sys:.1f}
    D_obs (observation dimensionality): 2
    Ratio D_sys / D_obs:               {d_sys/2:.1f}x

    Topological aliasing rate:         {aliasing_rate:.1%}
    Coverage of state space:           {coverage:.2e}

    KEY INSIGHT:

    The gene expression manifold has ~{d_sys:.0f} effective dimensions,
    but we routinely visualize it in 2D. This ~{d_sys/2:.0f}x compression
    introduces {aliasing_rate:.0%} aliasing: cells that are neighbors in
    the 2D plot often weren't neighbors in the original space.

    This is not a failure of t-SNE/UMAP—it's a geometric inevitability
    when D_sys >> D_obs. The "structure" we see in 2D plots is partially
    real and partially artifact.

    The falsifiability implication: binary classifications based on
    2D cluster assignments have a ~{aliasing_rate:.0%} baseline error rate
    due to topological aliasing alone.
    """)


if __name__ == '__main__':
    run_simulation()
