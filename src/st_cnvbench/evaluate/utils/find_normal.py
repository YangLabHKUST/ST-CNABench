import pandas as pd
import scanpy as sc
import numpy as np
from sklearn.cluster import KMeans

def find_normal_tumor(ref_df, obs_df, mode='ref', baseline=1.0):
    """
    This function performs binary classification into 'normal' and 'tumor' using KMeans clustering 
    in PCA space based on provided reference and observation CNV profile DataFrames.
    mode: ref means both ref and obs are provided; no_ref means only obs is provided.
    """
    # Standardize orientation: Ensure cells are rows, genes are columns
    # Assuming input is Gene x Cell based on standard InferCNV format
    if ref_df is not None and mode == 'ref':
        df_ref = ref_df.T
        df_obs = obs_df.T

        common_genes = df_ref.columns.intersection(df_obs.columns)
        df_ref = df_ref[common_genes].copy()
        df_obs = df_obs[common_genes].copy()

        df_ref['source'] = 'ref'
        df_obs['source'] = 'obs'
        df_all = pd.concat([df_ref, df_obs])
    else:
        df_obs = obs_df.T
        df_obs['source'] = 'obs'
        df_all = df_obs.copy()

    # Construct AnnData for computation
    adata = sc.AnnData(df_all.drop(columns='source').values.astype(np.float32))
    adata.obs_names = df_all.index
    adata.obs['source'] = df_all['source'].values

    # Dimensionality reduction
    sc.tl.pca(adata, svd_solver='arpack')

    # Binary clustering
    kmeans = KMeans(n_clusters=2, random_state=42)
    clusters = kmeans.fit_predict(adata.obsm['X_pca'])
    adata.obs['cluster'] = clusters.astype(str)

    # Determine which cluster is normal based on cluster cnv profile
    cluster_score = {}
    for cluster in adata.obs['cluster'].unique():
        cells_in_cluster = adata.obs_names[adata.obs['cluster'] == cluster]
        cluster_data = adata[cells_in_cluster].X
        diff = np.mean((cluster_data - baseline)**2)
        cluster_score[cluster] = diff
    normal_cluster = min(cluster_score, key=cluster_score.get)
    print(f"Identified normal cluster: {normal_cluster} with score {cluster_score}")
    adata.obs['label'] = adata.obs['cluster'].apply(
        lambda x: 'normal' if x == normal_cluster else 'tumor'
    )
    if mode == 'ref':
        preds_series = adata.obs.loc[adata.obs['source'] == 'obs', 'label']
    else:
        preds_series = adata.obs['label']
    return preds_series