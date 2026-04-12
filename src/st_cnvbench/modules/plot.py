import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def create_spatial_expression_plot(count_matrix_file, coordinate_file, output_file='spatial_expression_plot.png', markersize=50):
    """
    Create a spatial plot showing total gene expression for each spot.
    
    Parameters:
    -----------
    count_matrix_file : str
        Path to the count matrix file (CSV or TSV format)
    coordinate_file : str
        Path to the coordinate file containing barcode spatial information
    output_file : str
        Output filename for the saved plot (default: 'spatial_expression_plot.png')
    
    Returns:
    --------
    matplotlib.figure.Figure
        The generated figure object
    """
    
    # Read count matrix file
    # Assuming the file is tab-separated with genes as rows and barcodes as columns
    print("Reading count matrix file...")
    count_df = pd.read_csv(count_matrix_file, sep='\t', index_col=0)

    # Transpose to have barcodes as rows and genes as columns
    count_df = count_df.T
    
    # Calculate total expression per spot (sum of all genes for each barcode)
    print("Calculating total expression per spot...")
    total_expression = count_df.sum(axis=1)
    
    # Read coordinate file
    # Assuming the file is space or tab separated with columns: barcode, imagerow, imagecol
    print("Reading coordinate file...")
    coord_df = pd.read_csv(coordinate_file, sep='\t')
    
    # Merge expression data with coordinate data
    print("Merging expression and coordinate data...")
    # Create a DataFrame from total expression
    expression_df = pd.DataFrame({
        'barcode': total_expression.index,
        'total_expression': total_expression.values
    })
    
    # Merge with coordinates
    merged_df = pd.merge(expression_df, coord_df, on='barcode', how='inner')
    
    # Create the spatial expression plot
    print("Creating spatial expression plot...")
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create scatter plot where color represents total expression
    scatter = ax.scatter(
        merged_df['imagecol'],           # x-coordinates (image columns)
        merged_df['imagerow'],           # y-coordinates (image rows)
        c=merged_df['total_expression'], # Color by total expression
        cmap='viridis',                  # Color map for expression levels
        s=markersize,                            # Marker size
    )
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
    cbar.set_label('Total Gene Expression', fontsize=12, rotation=270, labelpad=15)
    
    # Customize plot appearance
    ax.set_xlabel('Image Column Coordinate', fontsize=12)
    ax.set_ylabel('Image Row Coordinate', fontsize=12)
    ax.set_title('Spatial Distribution of Total Gene Expression', fontsize=14, pad=20)
    plt.grid(True, alpha=0.3)
    
    # Invert y-axis to match typical image coordinates (origin at top-left)
    ax.invert_yaxis()
    
    # Customize ticks
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    plt.tight_layout()
    
    # Save the plot
    print(f"Saving plot to {output_file}...")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    
    print("Plot created successfully!")
    
    # Display some basic statistics
    print(f"\nDataset Statistics:")
    print(f"Number of spots: {len(merged_df)}")
    print(f"Number of genes: {count_df.shape[1]}")
    print(f"Total expression range: {merged_df['total_expression'].min():.0f} - {merged_df['total_expression'].max():.0f}")
    print(f"Average expression per spot: {merged_df['total_expression'].mean():.2f}")
    
    return fig

def plot_tumor_spatial_distribution(coord_file, tumor_file, output_plot='tumor_spatial_plot.png', marker_size=60):
    """
    Plot spatial distribution of barcodes colored by tumor/normal classification.
    
    Parameters:
    -----------
    coord_file : str
        Path to file containing barcode coordinates with columns: barcode, imagerow, imagecol
    tumor_file : str
        Path to file containing tumor classification with columns: Barcode, tumor_normal
    output_plot : str, optional
        Output filename for the plot (default: 'tumor_spatial_plot.png')
    
    Returns:
    --------
    matplotlib.figure.Figure
        The generated figure object
    """
    
    # Read coordinate data from the first file
    # Assuming the data is space-separated with specific column structure
    coord_df = pd.read_csv(coord_file, sep=r'\s+', header=0, 
                          names=['barcode', 'imagerow', 'imagecol'])
    
    # Read tumor classification data from the second file  
    tumor_df = pd.read_csv(tumor_file, sep=r'\s+', header=0,
                          names=['Barcode', 'tumor_normal'])
    
    # Merge the two dataframes on barcode
    merged_df = pd.merge(coord_df, tumor_df, 
                        left_on='barcode', right_on='Barcode', 
                        how='inner')
    
    # Create the plot
    plt.figure(figsize=(10, 10))
    
    # Separate tumor and normal samples for plotting
    tumor_samples = merged_df[merged_df['tumor_normal'] == 'Tumor']
    normal_samples = merged_df[merged_df['tumor_normal'] == 'Normal']
    
    # Plot tumor samples in red
    plt.scatter(tumor_samples['imagecol'], tumor_samples['imagerow'], 
               c='red', s=marker_size, label='Tumor')
    
    # Plot normal samples in blue
    plt.scatter(normal_samples['imagecol'], normal_samples['imagerow'], 
               c='blue', s=marker_size, label='Normal')
    
    # Customize the plot
    plt.ylabel('Image Row Coordinate', fontsize=12)
    plt.xlabel('Image Column Coordinate', fontsize=12)
    plt.title('Spatial Distribution of Tumor vs Normal Samples', fontsize=14, pad=20)
    plt.legend(title='Sample Type', title_fontsize=11, fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.gca().invert_yaxis()
    
    # Add statistics to the plot
    tumor_count = len(tumor_samples)
    normal_count = len(normal_samples)
    total_count = tumor_count + normal_count
    tumor_percentage = (tumor_count / total_count) * 100
    
    # Add text box with statistics
    stats_text = f'Total samples: {total_count}\nTumor: {tumor_count} ({tumor_percentage:.1f}%)\nNormal: {normal_count}'
    plt.annotate(stats_text, xy=(0.02, 0.98), xycoords='axes fraction',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                fontsize=10, verticalalignment='top')
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Plot saved as {output_plot}")
    print(f"Dataset statistics: {tumor_count} tumor samples, {normal_count} normal samples")
    
    return plt.gcf()



def plot_clone_clusters(cluster_coord_file, barcode_coord_file, barcode_position_file, output_file=None, marker_size=70):
    """
    Plot spatial transcriptomics clone clusters using absolute positions.
    Non-tumor spots (not in clusters) are shown in gray.
    
    Parameters:
    -----------
    cluster_coord_file : str
        Path to file mapping coordinates to cluster IDs (format: xxxtimesyyy,cluster)
    barcode_coord_file : str
        Path to file mapping barcodes to coordinate pairs (format: barcode xxxtimesyyy)
    barcode_position_file : str
        Path to file mapping barcodes to absolute positions (format: barcode imagerow imagecol value)
    output_file : str, optional
        Path to save the plot. If None, plot is displayed.
    """
    
    # Load and process cluster-coordinate mapping
    # Format: "50x102,0" - coordinate pair and cluster ID
    cluster_data = pd.read_csv(cluster_coord_file, names=['coord', 'cluster'], skiprows=1)
    cluster_data['cluster'] = cluster_data['cluster'].astype(int)
    
    # Load barcode-coordinate mapping
    # Format: barcodes and corresponding coordinate pairs like "50x102"
    barcode_coord = pd.read_csv(barcode_coord_file, sep='\t', header=0)
    barcode_coord.columns = ['barcode', 'coord']
    barcode_coord['barcode'] = barcode_coord['barcode'].str.replace('_U1$', '', regex=True)
    
    # Load barcode absolute position data
    # Format: barcode, imagerow, imagecol, and numeric value
    barcode_pos = pd.read_csv(barcode_position_file, sep=r'\s+', header=0)
    
    # Merge barcode-coordinate with cluster data
    barcode_cluster = pd.merge(barcode_coord, cluster_data[['coord', 'cluster']], 
                              on='coord', how='inner')
    
    # Merge with absolute position data to get cluster spots
    cluster_spots_data = pd.merge(barcode_cluster, barcode_pos[['barcode', 'imagerow', 'imagecol']],
                                 on='barcode', how='inner')
    
    # Get all spots from barcode position data (including non-tumor spots)
    all_spots_data = barcode_pos[['barcode', 'imagerow', 'imagecol']].copy()
    
    # Identify non-tumor spots (spots not in cluster_spots_data)
    non_tumor_spots = all_spots_data[~all_spots_data['barcode'].isin(cluster_spots_data['barcode'])]
    
    # Create the plot
    plt.figure(figsize=(12, 12))
    
    # First plot non-tumor spots in gray (background)
    if len(non_tumor_spots) > 0:
        plt.scatter(non_tumor_spots['imagecol'], non_tumor_spots['imagerow'],
                   c='lightgray', alpha=0.7, label='Non-tumor',
                   s=marker_size, edgecolors='none')
    
    # Then plot cluster spots with colors
    unique_clusters = cluster_spots_data['cluster'].unique()
    colors = plt.cm.Set3(np.linspace(0, 1, len(unique_clusters)))
    color_map = dict(zip(unique_clusters, colors))
    
    # Plot each cluster with specific coloring
    for cluster_id in unique_clusters:
        cluster_spots = cluster_spots_data[cluster_spots_data['cluster'] == cluster_id]
        plt.scatter(cluster_spots['imagecol'], cluster_spots['imagerow'],
                   c=[color_map[cluster_id]], 
                   label=f'Cluster {cluster_id}',
                   s=marker_size)
    
    # Customize plot appearance
    plt.xlabel('Image Column Coordinate', fontsize=12)
    plt.ylabel('Image Row Coordinate', fontsize=12)
    plt.title('Spatial Transcriptomics Clone Clusters', 
              fontsize=14, pad=20)
    plt.legend(title='Spots', loc='upper left')
    plt.grid(True, alpha=0.3)

    plt.gca().invert_yaxis()
    
    # Save or display the plot
    if output_file:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()
    
    # Return both cluster spots and information about non-tumor spots
    result_info = {
        'cluster_spots': cluster_spots_data,
        'non_tumor_spots': non_tumor_spots,
        'total_spots': len(all_spots_data),
        'cluster_spot_count': len(cluster_spots_data),
        'non_tumor_spot_count': len(non_tumor_spots)
    }
    
    return result_info

def plot_cnv_heatmap(filename, output_path=None, dpi=300):
    """
    Read CNV data from file and plot a heatmap.
    
    Parameters:
    filename (str): Path to the CSV file containing CNV data
    output_path (str, optional): Path to save the heatmap image. If None, only display the plot.
    dpi (int): Resolution for saving the image (default: 300)
    
    Returns:
    pd.DataFrame: The processed CNV data
    """
    # Read data from file
    df = pd.read_csv(filename)
    
    # Set first column (gene names) as index
    df.set_index(df.columns[0], inplace=True)
    
    # Automatically detect clone names from column headers (all columns except the first)
    clone_columns = df.columns.tolist()
    
    # Transpose the dataframe to make the heatmap horizontal
    df_transposed = df.T
    
    # Create the heatmap with horizontal orientation
    plt.figure(figsize=(10, 6))  # Adjusted dimensions for horizontal layout
    
    # Define custom colormap
    cmap = sns.color_palette(["blue", "white", "red"], as_cmap=True)
    
    # Plot heatmap without colorbar
    heatmap = sns.heatmap(df_transposed, 
                         cmap=cmap,
                         vmin=0, 
                         vmax=2,
                         cbar=False,  # Disable default colorbar
                         fmt='d')
    
    # Customize plot
    plt.title('CNV Heatmap', fontsize=14, pad=20)
    plt.xlabel('Genes', fontsize=12)
    plt.ylabel('Clone Clusters', fontsize=12)
    
    # Create custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue', label='Deletion (0)'),
        Patch(facecolor='white', label='Normal (1)'),
        Patch(facecolor='red', label='Amplification (2)')
    ]
    
    # Add legend to the plot
    plt.legend(handles=legend_elements, 
               title='CNV Status',
               loc='upper right',
               bbox_to_anchor=(1.3, 1),  # Position legend to the right of heatmap
               frameon=True)
    
    # Add border around the heatmap
    for _, spine in heatmap.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_edgecolor('black')
    
    plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    # Save the plot if output_path is provided
    if output_path:
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"Heatmap saved to: {output_path}")
    
    plt.show()
    return df