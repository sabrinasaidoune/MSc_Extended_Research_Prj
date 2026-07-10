import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# PATHWAY CLUSTERING HEATMAP FUNCTIONS FOR MULTI-OMICS DATA
# ============================================================

def load_multiomics_data(file_path):
    """Load multi-omics data from CSV file"""
    df = pd.read_csv(file_path)
    print(f"✓ Loaded {file_path} with {len(df)} genes/proteins")
    return df


def extract_expression_matrix(df, timepoint_cols=None):
    """Extract expression values and create numeric matrix"""
    if timepoint_cols is None:
        timepoint_cols = ['mean_log2_D14', 'mean_log2_D21', 'mean_log2_M03',
                          'mean_log2_M06', 'mean_log2_M12', 'mean_log2_M18']

    available_cols = [col for col in timepoint_cols if col in df.columns]

    if not available_cols:
        raise ValueError(f"No timepoint columns found. Available columns: {df.columns.tolist()}")

    expr_data = df[available_cols].copy()

    if 'gene_symbol' in df.columns:
        gene_names = df['gene_symbol'].values
    else:
        gene_names = df.iloc[:, 0].values

    expr_matrix = expr_data.apply(pd.to_numeric, errors='coerce').values
    timepoint_labels = [col.replace('mean_log2_', '') for col in available_cols]

    print(f"  Expression matrix shape: {expr_matrix.shape}")
    print(f"  Timepoints: {timepoint_labels}")

    return expr_matrix, gene_names, timepoint_labels


def normalize_expression_data(expr_matrix):
    """Z-score normalize expression data across genes"""
    print(f"  Normalizing expression data...")

    valid_idx = ~np.all(np.isnan(expr_matrix), axis=1)
    expr_matrix_valid = expr_matrix[valid_idx].copy()

    for i in range(expr_matrix_valid.shape[0]):
        row = expr_matrix_valid[i, :]
        row_mean = np.nanmean(row)
        if not np.isnan(row_mean):
            expr_matrix_valid[i, np.isnan(row)] = row_mean

    expr_normalized = np.zeros_like(expr_matrix_valid)
    for i in range(expr_matrix_valid.shape[0]):
        row = expr_matrix_valid[i, :]
        row_mean = np.mean(row)
        row_std = np.std(row)
        if row_std > 0:
            expr_normalized[i, :] = (row - row_mean) / row_std
        else:
            expr_normalized[i, :] = 0

    print(f"  Retained {expr_normalized.shape[0]} genes with valid data")

    return expr_normalized, valid_idx


def cluster_genes_by_pattern(expr_matrix, gene_names, n_clusters=5):
    """Cluster genes into groups based on expression patterns"""
    print(f"  Clustering {len(gene_names)} genes into {n_clusters} clusters...")

    linkage_matrix = linkage(expr_matrix, method='ward')
    clusters = fcluster(linkage_matrix, n_clusters, criterion='maxclust')

    cluster_df = pd.DataFrame({
        'Gene': gene_names,
        'Cluster': clusters
    })

    print(f"  Cluster distribution:")
    print(cluster_df['Cluster'].value_counts().sort_index())

    return cluster_df


def interpret_clusters(expr_df, cluster_df, original_df,
                      go_bp_col='go_biological_process',
                      go_mf_col='go_molecular_function',
                      kegg_col='kegg_pathway',
                      timepoint_labels=None):
    """
    Analyze and interpret what each cluster represents based on GO terms,
    biological processes, and expression patterns

    Args:
        expr_df: DataFrame with normalized expression values
        cluster_df: DataFrame with cluster assignments
        original_df: Original dataframe with GO/KEGG annotations
        go_bp_col: Column name for GO biological process
        go_mf_col: Column name for GO molecular function
        kegg_col: Column name for KEGG pathway
        timepoint_labels: List of timepoint labels

    Returns:
        Dictionary with cluster interpretations
    """
    print(f"\n  Interpreting cluster functions and pathways...")

    # Merge cluster assignments with original data
    merged = pd.merge(cluster_df, original_df, left_on='Gene', right_on='gene_symbol', how='left')

    cluster_interpretations = {}

    for cluster_id in sorted(merged['Cluster'].unique()):
        cluster_genes = merged[merged['Cluster'] == cluster_id]
        n_genes = len(cluster_genes)

        print(f"\n  {'='*70}")
        print(f"  CLUSTER {cluster_id} ({n_genes} genes/proteins)")
        print(f"  {'='*70}")

        # Get genes in cluster
        gene_list = cluster_genes['Gene'].tolist()
        print(f"  Members: {', '.join(gene_list[:5])}" +
              (f", ... +{len(gene_list)-5} more" if len(gene_list) > 5 else ""))

        # Extract and analyze GO terms
        go_bp_terms = []
        for terms_str in cluster_genes[go_bp_col].dropna():
            if isinstance(terms_str, str):
                terms = [t.strip() for t in terms_str.split(';') if t.strip()]
                go_bp_terms.extend(terms)

        go_mf_terms = []
        for terms_str in cluster_genes[go_mf_col].dropna():
            if isinstance(terms_str, str):
                terms = [t.strip() for t in terms_str.split(';') if t.strip()]
                go_mf_terms.extend(terms)

        kegg_terms = []
        for terms_str in cluster_genes[kegg_col].dropna():
            if isinstance(terms_str, str):
                terms = [t.strip() for t in terms_str.split(';') if t.strip()]
                kegg_terms.extend(terms)

        # Get top terms
        top_bp = Counter(go_bp_terms).most_common(5)
        top_mf = Counter(go_mf_terms).most_common(5)
        top_kegg = Counter(kegg_terms).most_common(3)

        # Expression pattern analysis
        cluster_expr = expr_df[expr_df.index.isin(gene_list)]
        cluster_mean = cluster_expr.mean()

        # Find peak expression timepoint
        peak_tp = cluster_mean.idxmax()
        peak_value = cluster_mean.max()

        # Find lowest expression timepoint
        min_tp = cluster_mean.idxmin()
        min_value = cluster_mean.min()

        print(f"\n  🎯 Expression Pattern:")
        print(f"     Peak expression: {peak_tp} (Z-score: {peak_value:.2f})")
        print(f"     Lowest expression: {min_tp} (Z-score: {min_value:.2f})")
        print(f"     Pattern: {'Early→Late activation' if peak_value > 0 and cluster_mean[0] < cluster_mean[-1] else 'Varies'}")

        print(f"\n  📚 Top Biological Processes (GO:BP):")
        for term, count in top_bp:
            if term:
                print(f"     • {term} ({count} proteins)")

        if not top_bp:
            print(f"     • No GO:BP terms annotated")

        print(f"\n  🔬 Top Molecular Functions (GO:MF):")
        for term, count in top_mf:
            if term:
                print(f"     • {term} ({count} proteins)")

        if not top_mf:
            print(f"     • No GO:MF terms annotated")

        print(f"\n  🛣️  KEGG Pathways:")
        for term, count in top_kegg:
            if term:
                print(f"     • {term} ({count} proteins)")

        if not top_kegg:
            print(f"     • No KEGG pathway annotations")

        # Generate cluster label
        if top_bp:
            main_process = top_bp[0][0]
            cluster_label = f"Cluster {cluster_id}: {main_process}"
        else:
            cluster_label = f"Cluster {cluster_id}: Uncharacterized"

        # Store interpretation
        cluster_interpretations[cluster_id] = {
            'n_genes': n_genes,
            'genes': gene_list,
            'go_bp': top_bp,
            'go_mf': top_mf,
            'kegg': top_kegg,
            'peak_timepoint': peak_tp,
            'peak_value': peak_value,
            'label': cluster_label,
            'expression_pattern': f"Peak at {peak_tp}"
        }

    return cluster_interpretations


def plot_cluster_interpretation_summary(cluster_interpretations, figsize=(16, 10), save_path=None):
    """
    Create a summary visualization of cluster interpretations
    """
    print(f"\n  Creating cluster interpretation summary...")

    n_clusters = len(cluster_interpretations)
    fig, axes = plt.subplots(n_clusters, 1, figsize=figsize)

    if n_clusters == 1:
        axes = [axes]

    for idx, (cluster_id, info) in enumerate(sorted(cluster_interpretations.items())):
        ax = axes[idx]

        # Get top processes
        processes = info['go_bp'][:5]
        labels = [p[0][:50] if p[0] else "Uncharacterized" for p in processes]
        counts = [p[1] for p in processes]

        # Create horizontal bar chart
        y_pos = np.arange(len(labels))
        colors = plt.cm.Set3(np.linspace(0, 1, len(labels)))

        ax.barh(y_pos, counts, color=colors, edgecolor='black', alpha=0.8)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel('Protein Count', fontsize=10, fontweight='bold')

        # Add title with cluster info
        title = f"{info['label']}\n({info['n_genes']} proteins, Peak: {info['peak_timepoint']}, Z={info['peak_value']:.2f})"
        ax.set_title(title, fontsize=11, fontweight='bold', pad=10)

        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved to {save_path}")

    plt.show()


def save_cluster_annotations(cluster_interpretations, output_path):
    """
    Save cluster interpretations to a CSV file
    """
    print(f"\n  Saving cluster annotations...")

    data = []
    for cluster_id, info in sorted(cluster_interpretations.items()):
        top_bp = '; '.join([p[0] for p in info['go_bp']]) if info['go_bp'] else "N/A"
        top_mf = '; '.join([p[0] for p in info['go_mf']]) if info['go_mf'] else "N/A"
        top_kegg = '; '.join([p[0] for p in info['kegg']]) if info['kegg'] else "N/A"

        data.append({
            'Cluster': cluster_id,
            'N_Genes': info['n_genes'],
            'Peak_Timepoint': info['peak_timepoint'],
            'Peak_Z_Score': info['peak_value'],
            'Top_GO_BP': top_bp,
            'Top_GO_MF': top_mf,
            'Top_KEGG': top_kegg,
            'Cluster_Label': info['label']
        })

    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
    print(f"  ✓ Saved to {output_path}")

    return df


def plot_clustering_heatmap(expr_matrix, gene_names, timepoint_labels,
                           title="Clustering Heatmap with Annotations",
                           figsize=(12, 14), save_path=None):
    """Create hierarchical clustering heatmap"""
    print(f"  Creating clustering heatmap ({expr_matrix.shape[0]} genes)...")

    hm_df = pd.DataFrame(
        expr_matrix,
        columns=timepoint_labels,
        index=gene_names
    )

    g = sns.clustermap(
        hm_df,
        method='ward',
        metric='euclidean',
        cmap='RdBu_r',
        center=0,
        vmin=-2.5,
        vmax=2.5,
        col_cluster=False,
        row_cluster=True,
        cbar_kws={'label': 'Z-score (Expression)'},
        figsize=figsize,
        dendrogram_ratio=(0.15, 0.15),
        linewidths=0,
        linecolor='white'
    )

    g.ax_heatmap.set_xlabel('Developmental Stage', fontsize=12, fontweight='bold')
    g.ax_heatmap.set_ylabel('Genes/Proteins', fontsize=12, fontweight='bold')
    g.ax_heatmap.set_xticklabels(timepoint_labels, rotation=45, ha='right')
    plt.suptitle(title, fontsize=14, fontweight='bold', y=0.995)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved to {save_path}")

    plt.show()

    return g


def plot_cluster_averages_heatmap(expr_matrix, gene_names, timepoint_labels,
                                  cluster_df, cluster_interpretations=None,
                                  title="Cluster Averages",
                                  figsize=(10, 8), save_path=None):
    """Create heatmap of average expression per cluster"""
    print(f"  Creating cluster averages heatmap...")

    expr_df = pd.DataFrame(expr_matrix, columns=timepoint_labels, index=gene_names)

    cluster_mapping = dict(zip(cluster_df['Gene'], cluster_df['Cluster']))
    expr_df['Cluster'] = expr_df.index.map(cluster_mapping)

    cluster_avg = expr_df.groupby('Cluster')[timepoint_labels].mean()

    # Create cluster labels if interpretations provided
    if cluster_interpretations:
        cluster_labels = [cluster_interpretations[c]['label'].split(': ')[1][:30]
                         for c in cluster_avg.index]
        y_labels = [f"C{c}\n{cluster_interpretations[c]['label'].split(': ')[1][:25]}"
                   for c in cluster_avg.index]
    else:
        y_labels = [f"Cluster {c}" for c in cluster_avg.index]

    plt.figure(figsize=figsize)
    sns.heatmap(
        cluster_avg,
        annot=True,
        fmt='.2f',
        cmap='RdBu_r',
        center=0,
        vmin=-2.5,
        vmax=2.5,
        cbar_kws={'label': 'Mean Z-score'},
        linewidths=1,
        linecolor='black',
        yticklabels=y_labels
    )

    plt.xlabel('Developmental Stage', fontsize=12, fontweight='bold')
    plt.ylabel('Cluster', fontsize=12, fontweight='bold')
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved to {save_path}")

    plt.show()

    return cluster_avg


def plot_all_clusters_together(expr_matrix, gene_names, timepoint_labels,
                               cluster_df, cluster_interpretations=None,
                               title="All Clusters - Temporal Patterns",
                               figsize=(12, 8), save_path=None):
    """Plot all cluster trajectories on same graph"""
    print(f"  Creating combined cluster comparison plot...")

    expr_df = pd.DataFrame(expr_matrix, columns=timepoint_labels, index=gene_names)

    cluster_mapping = dict(zip(cluster_df['Gene'], cluster_df['Cluster']))
    expr_df['Cluster'] = expr_df.index.map(cluster_mapping)

    cluster_avg = expr_df.groupby('Cluster')[timepoint_labels].mean()

    n_clusters = cluster_avg.shape[0]
    colors = plt.cm.tab10(np.linspace(0, 1, n_clusters))

    fig, ax = plt.subplots(figsize=figsize)

    x_pos = np.arange(len(timepoint_labels))

    for cluster_id in range(1, n_clusters + 1):
        cluster_mean = cluster_avg.loc[cluster_id].values
        cluster_genes = expr_df[expr_df['Cluster'] == cluster_id]
        n_genes = len(cluster_genes)

        if cluster_interpretations and cluster_id in cluster_interpretations:
            label = f"C{cluster_id}: {cluster_interpretations[cluster_id]['label'].split(': ')[1][:30]} (n={n_genes})"
        else:
            label = f"Cluster {cluster_id} (n={n_genes})"

        ax.plot(x_pos, cluster_mean, color=colors[cluster_id - 1],
               linewidth=3, marker='o', markersize=10, label=label)

    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(timepoint_labels, fontsize=11, fontweight='bold')
    ax.set_ylabel('Z-score (Expression)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Developmental Stage', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10, framealpha=0.95)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.set_ylim(-2, 2)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved to {save_path}")

    plt.show()


# ============================================================
# MAIN ANALYSIS PIPELINE
# ============================================================

if __name__ == '__main__':

    import os
    os.makedirs('results', exist_ok=True)

    print("\n" + "="*70)
    print("MULTI-OMICS PATHWAY CLUSTERING ANALYSIS WITH INTERPRETATION")
    print("="*70)

    print("\n>>> LOADING DATA <<<")
    cellular_df = load_multiomics_data('multi_omics_significant_cellular.csv')
    soluble_df = load_multiomics_data('multi_omics_significant_soluble.csv')

    # ---- CELLULAR ANALYSIS ----
    print("\n" + "="*70)
    print("CELLULAR PROTEOMICS ANALYSIS")
    print("="*70)

    expr_matrix_cell, genes_cell, timepoints_cell = extract_expression_matrix(cellular_df)
    expr_norm_cell, valid_idx_cell = normalize_expression_data(expr_matrix_cell)
    valid_genes_cell = genes_cell[valid_idx_cell]
    valid_cellular_df = cellular_df[valid_idx_cell].reset_index(drop=True)
    valid_cellular_df['gene_symbol'] = valid_genes_cell

    print("\n>>> Clustering Heatmap <<<")
    plot_clustering_heatmap(
        expr_norm_cell, valid_genes_cell, timepoints_cell,
        title="Cellular Proteins - Clustering Heatmap",
        save_path='results/heatmap_cellular_clustering.png'
    )

    print("\n>>> Gene Clustering <<<")
    cluster_df_cell = cluster_genes_by_pattern(expr_norm_cell, valid_genes_cell, n_clusters=6)

    # Create expression dataframe for interpretation
    expr_df_cell = pd.DataFrame(expr_norm_cell, columns=timepoints_cell, index=valid_genes_cell)

    print("\n>>> Interpreting Clusters <<<")
    cluster_interp_cell = interpret_clusters(
        expr_df_cell, cluster_df_cell, valid_cellular_df,
        go_bp_col='go_biological_process',
        go_mf_col='go_molecular_function',
        kegg_col='kegg_pathway',
        timepoint_labels=timepoints_cell
    )

    print("\n>>> Cluster Averages Heatmap <<<")
    plot_cluster_averages_heatmap(
        expr_norm_cell, valid_genes_cell, timepoints_cell, cluster_df_cell,
        cluster_interpretations=cluster_interp_cell,
        title="Cellular Proteins - Cluster Averages",
        save_path='results/heatmap_cellular_cluster_averages.png'
    )

    print("\n>>> All Clusters Comparison <<<")
    plot_all_clusters_together(
        expr_norm_cell, valid_genes_cell, timepoints_cell, cluster_df_cell,
        cluster_interpretations=cluster_interp_cell,
        title="Cellular Proteins - All Clusters Temporal Comparison",
        save_path='results/cluster_comparison_cellular.png'
    )

    print("\n>>> Cluster Interpretation Summary <<<")
    plot_cluster_interpretation_summary(
        cluster_interp_cell,
        save_path='results/cluster_interpretation_cellular.png'
    )

    print("\n>>> Saving Cluster Annotations <<<")
    cluster_anno_cell = save_cluster_annotations(
        cluster_interp_cell,
        'results/cluster_annotations_cellular.csv'
    )

    # ---- SOLUBLE ANALYSIS ----
    print("\n" + "="*70)
    print("SOLUBLE PROTEOMICS ANALYSIS")
    print("="*70)

    expr_matrix_sol, genes_sol, timepoints_sol = extract_expression_matrix(soluble_df)
    expr_norm_sol, valid_idx_sol = normalize_expression_data(expr_matrix_sol)
    valid_genes_sol = genes_sol[valid_idx_sol]
    valid_soluble_df = soluble_df[valid_idx_sol].reset_index(drop=True)
    valid_soluble_df['gene_symbol'] = valid_genes_sol

    print("\n>>> Clustering Heatmap <<<")
    plot_clustering_heatmap(
        expr_norm_sol, valid_genes_sol, timepoints_sol,
        title="Soluble Proteins - Clustering Heatmap",
        save_path='results/heatmap_soluble_clustering.png'
    )

    print("\n>>> Gene Clustering <<<")
    cluster_df_sol = cluster_genes_by_pattern(expr_norm_sol, valid_genes_sol, n_clusters=6)

    expr_df_sol = pd.DataFrame(expr_norm_sol, columns=timepoints_sol, index=valid_genes_sol)

    print("\n>>> Interpreting Clusters <<<")
    cluster_interp_sol = interpret_clusters(
        expr_df_sol, cluster_df_sol, valid_soluble_df,
        go_bp_col='go_biological_process',
        go_mf_col='go_molecular_function',
        kegg_col='kegg_pathway',
        timepoint_labels=timepoints_sol
    )

    print("\n>>> Cluster Averages Heatmap <<<")
    plot_cluster_averages_heatmap(
        expr_norm_sol, valid_genes_sol, timepoints_sol, cluster_df_sol,
        cluster_interpretations=cluster_interp_sol,
        title="Soluble Proteins - Cluster Averages",
        save_path='results/heatmap_soluble_cluster_averages.png'
    )

    print("\n>>> All Clusters Comparison <<<")
    plot_all_clusters_together(
        expr_norm_sol, valid_genes_sol, timepoints_sol, cluster_df_sol,
        cluster_interpretations=cluster_interp_sol,
        title="Soluble Proteins - All Clusters Temporal Comparison",
        save_path='results/cluster_comparison_soluble.png'
    )

    print("\n>>> Cluster Interpretation Summary <<<")
    plot_cluster_interpretation_summary(
        cluster_interp_sol,
        save_path='results/cluster_interpretation_soluble.png'
    )

    print("\n>>> Saving Cluster Annotations <<<")
    cluster_anno_sol = save_cluster_annotations(
        cluster_interp_sol,
        'results/cluster_annotations_soluble.csv'
    )

    print("\n" + "="*70)
    print("✓ ANALYSIS COMPLETE!")
    print("="*70)
    print("\nOutput files:")
    print("  • Heatmaps: heatmap_*_clustering.png")
    print("  • Cluster patterns: cluster_*_soluble/cellular.png")
    print("  • Interpretations: cluster_interpretation_*.png")
    print("  • Annotations: cluster_annotations_*.csv")
    print("\nAll files saved to 'results/' directory\n")
