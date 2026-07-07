import os
import logging
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patheffects as path_effects
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import gaussian_kde
from scipy.spatial import cKDTree


def clean_label(val):
    if pd.isna(val) or str(val).lower() == 'nan':
        return 'nan'
    s = str(val).strip()
    return s[:-2] if s.endswith('.0') else s


def get_density_peak(x_coords, y_coords):
    if len(x_coords) == 0:
        return np.nan, np.nan
    if len(x_coords) < 5:
        cx, cy = np.mean(x_coords), np.mean(y_coords)
        distances = (x_coords - cx) ** 2 + (y_coords - cy) ** 2
        min_idx = np.argmin(distances.values)
        return x_coords.iloc[min_idx], y_coords.iloc[min_idx]
    try:
        xy = np.vstack([x_coords, y_coords])
        kde = gaussian_kde(xy)
        density = kde.evaluate(xy)
        max_idx = np.argmax(density)
        return x_coords.iloc[max_idx], y_coords.iloc[max_idx]
    except Exception:
        cx, cy = np.mean(x_coords), np.mean(y_coords)
        distances = (x_coords - cx) ** 2 + (y_coords - cy) ** 2
        min_idx = np.argmin(distances.values)
        return x_coords.iloc[min_idx], y_coords.iloc[min_idx]


def get_root_node(G):
    roots = [n for n in G.nodes() if G.in_degree(n) == 0]
    if not roots:
        return None
    return roots[0]


def scale_size(count, min_count, max_count, min_size, max_size):
    if max_count <= min_count:
        return (min_size + max_size) / 2.0
    norm = (np.sqrt(count) - np.sqrt(min_count)) / (np.sqrt(max_count) - np.sqrt(min_count))
    return min_size + norm * (max_size - min_size)


def estimate_spot_size(plot_data, base_area=36.0, ref_dist=1.0, min_area=12.0, max_area=80.0):
    coords = plot_data[['x', 'y']].dropna().values
    if len(coords) < 2:
        return base_area
    try:
        tree = cKDTree(coords)
        dists, _ = tree.query(coords, k=2)
        nn = dists[:, 1]
        nn = nn[np.isfinite(nn)]
        if len(nn) == 0:
            return base_area
        med = float(np.median(nn))
        if not np.isfinite(med) or med <= 0:
            return base_area
        area = base_area * (med / ref_dist) ** 2
        return float(np.clip(area, min_area, max_area))
    except Exception:
        return base_area



def get_slice_id_from_barcode(barcode):
    if pd.isna(barcode):
        return None
    s = str(barcode)
    if '_' in s:
        parts = s.split('_')
        if len(parts) >= 2:
            return f"{parts[0]}_{parts[1]}"
        return parts[0]
    if '-' in s:
        return s.split('-')[0]
    return s


def ensure_slice_id(plot_data):
    if 'Barcodes' not in plot_data.columns:
        return plot_data
    plot_data = plot_data.copy()
    if 'slice_id' in plot_data.columns:
        missing = plot_data['slice_id'].isna()
        if missing.any():
            plot_data.loc[missing, 'slice_id'] = plot_data.loc[missing, 'Barcodes'].apply(get_slice_id_from_barcode)
        return plot_data
    plot_data['slice_id'] = plot_data['Barcodes'].apply(get_slice_id_from_barcode)
    return plot_data


def compute_spatial_anchors(
    G,
    plot_data,
    pred_col,
    root_name=None,
    internal_offset_frac=0.03,
    blend_to_children=0.3,
    internal_lift_frac=0.12,
    internal_left_frac=0.08,
    depth_left_scale=0.3,
):
    """
    Estimate internal-node positions with CalicoST-style Gaussian conditioning.
    Then pull each internal node slightly toward the mean of its children and
    shift the tree toward the upper-left to keep ancestor branches readable.
    """
    pos = {}

    root_name = root_name or get_root_node(G)
    if root_name is None:
        return pos

    global_mean_x = plot_data['x'].mean()
    global_mean_y = plot_data['y'].mean()
    x_span = plot_data['x'].max() - plot_data['x'].min()
    y_span = plot_data['y'].max() - plot_data['y'].min()
    internal_offset = y_span * internal_offset_frac
    lift_offset = y_span * internal_lift_frac
    left_offset = x_span * internal_left_frac

    observed_nodes = [n for n in G.nodes() if G.nodes[n].get('is_observed', False)]
    internal_nodes = [n for n in G.nodes() if n not in observed_nodes]

    known_nodes = []
    known_coords = []

    # Coordinates for observed leaf nodes.
    for node in observed_nodes:
        subset = plot_data[plot_data['clean_pred'] == clean_label(node)]
        if len(subset) > 0:
            mx, my = get_density_peak(subset['x'], subset['y'])
            known_nodes.append(node)
            known_coords.append([mx, my])
        else:
            known_nodes.append(node)
            known_coords.append([global_mean_x, global_mean_y])

    for n, coord in zip(known_nodes, known_coords):
        pos[n] = tuple(coord)

    if not internal_nodes:
        return pos

    # =====================================================================
    # CalicoST-style Gaussian conditioning for internal nodes
    # =====================================================================
    try:
        dist_from_root = nx.single_source_dijkstra_path_length(G, root_name, weight='weight')

        all_nodes = observed_nodes + internal_nodes
        N_nodes = len(all_nodes)

        Sigma_square = np.zeros((N_nodes, N_nodes))
        max_abs_x = np.max(np.abs(plot_data['x']))
        max_abs_y = np.max(np.abs(plot_data['y']))
        base_var = max(max_abs_x, max_abs_y)

        for i, n1 in enumerate(all_nodes):
            for j, n2 in enumerate(all_nodes):
                if i == j:
                    Sigma_square[i, j] = base_var + dist_from_root.get(n1, 0)
                else:
                    try:
                        lca = nx.lowest_common_ancestor(G, n1, n2)
                        if lca == root_name:
                            Sigma_square[i, j] = base_var
                        else:
                            Sigma_square[i, j] = base_var + dist_from_root.get(lca, 0)
                    except Exception:
                        Sigma_square[i, j] = base_var

        k = len(observed_nodes)
        Sigma_11 = Sigma_square[:k, :k]
        Sigma_12 = Sigma_square[:k, k:]

        obs_1 = np.array(known_coords)
        mu_1 = np.zeros((k, 2))
        mu_2 = np.zeros((len(internal_nodes), 2))

        Sigma_11 += np.eye(k) * 1e-6

        expected_internal = mu_2 + Sigma_12.T @ (np.linalg.inv(Sigma_11) @ (obs_1 - mu_1))

        for idx, n in enumerate(internal_nodes):
            pos[n] = (expected_internal[idx, 0], expected_internal[idx, 1])

    except Exception as e:
        logging.error(f"Gaussian conditioning failed: {e}. Falling back to child mean.")
        for node in reversed(list(nx.topological_sort(G))):
            if node not in pos:
                children = list(G.successors(node))
                valid_coords = [pos[c] for c in children if c in pos and not np.isnan(pos[c][0])]
                if valid_coords:
                    pos[node] = (np.mean([c[0] for c in valid_coords]), np.mean([c[1] for c in valid_coords]))
                else:
                    pos[node] = (global_mean_x, global_mean_y)

    # Pull internal nodes toward their children and shift the tree left/up.
    depths = nx.single_source_shortest_path_length(G, root_name)
    for node in internal_nodes:
        children = list(G.successors(node))
        child_coords = [pos[c] for c in children if c in pos and not np.isnan(pos[c][0])]
        depth = depths.get(node, 1)
        depth_scale = 1.0 + depth_left_scale * depth

        if child_coords:
            cx = float(np.mean([c[0] for c in child_coords]))
            cy = float(np.mean([c[1] for c in child_coords])) - internal_offset
            px, py = pos.get(node, (cx, cy))
            pos[node] = (
                px * (1 - blend_to_children) + cx * blend_to_children - left_offset * depth_scale,
                (py * (1 - blend_to_children) + cy * blend_to_children) - lift_offset
            )
        else:
            px, py = pos.get(node, (global_mean_x, global_mean_y))
            pos[node] = (px - left_offset * depth_scale, py - lift_offset)

    return pos


def _tree_layout_positions(G, root_name):
    if root_name is None or root_name not in G:
        return {}

    children = {n: list(G.successors(n)) for n in G.nodes()}
    for n in children:
        children[n] = sorted(children[n], key=lambda x: str(x))

    # X positions by cumulative branch length
    x_pos = {root_name: 0.0}
    queue = [root_name]
    while queue:
        node = queue.pop(0)
        for child in children.get(node, []):
            w = G.edges[node, child].get('weight', 1.0)
            if not np.isfinite(w):
                w = 1.0
            x_pos[child] = x_pos.get(node, 0.0) + float(w)
            queue.append(child)

    # Y positions by leaf order
    y_pos = {}
    counter = [0]

    def assign_y(node):
        kids = children.get(node, [])
        if not kids:
            y_pos[node] = float(counter[0])
            counter[0] += 1
            return y_pos[node]
        vals = [assign_y(k) for k in kids]
        y_pos[node] = float(np.mean(vals))
        return y_pos[node]

    assign_y(root_name)

    pos = {}
    for node in G.nodes():
        pos[node] = (x_pos.get(node, 0.0), y_pos.get(node, 0.0))
    return pos


def _draw_tree_inset(ax, G, root_name, clone_to_label, color_dict, clone_counts, min_count, max_count):
    if root_name is None or root_name not in G:
        return

    pos = _tree_layout_positions(G, root_name)
    if not pos:
        return

    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_span = max(x_max - x_min, 1.0)
    y_span = max(y_max - y_min, 1.0)

    # Normalize to [0,1] for inset plotting
    norm_pos = {
        n: ((p[0] - x_min) / x_span, (p[1] - y_min) / y_span)
        for n, p in pos.items()
    }

    ax_inset = ax.inset_axes([0.02, 0.62, 0.34, 0.34])
    for u, v in G.edges():
        if u not in norm_pos or v not in norm_pos:
            continue
        x0, y0 = norm_pos[u]
        x1, y1 = norm_pos[v]
        ax_inset.plot([x0, x1], [y0, y1], color='black', linewidth=1.0)

        labels = G.edges[u, v].get('edge_labels', [])
        if labels:
            lx = 0.5 * (x0 + x1)
            ly = 0.5 * (y0 + y1)
            text = "\n".join(labels)
            ax_inset.text(
                lx, ly, text,
                ha='center',
                va='center',
                fontsize=6,
                bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.7),
                zorder=5
            )

    for node, (x, y) in norm_pos.items():
        is_observed = G.nodes[node].get('is_observed', G.nodes[node].get('is_leaf', False))
        if is_observed:
            clean_node_name = clean_label(node)
            label = clone_to_label.get(clean_node_name)
            node_color = color_dict.get(label, 'red')
            count_val = clone_counts.get(clean_node_name, 1)
            size = scale_size(count_val, min_count, max_count, 160, 420)
            ax_inset.scatter([x], [y], s=size, c=[node_color], edgecolors='black', linewidths=0.8, zorder=3)
            if label is not None:
                ax_inset.text(
                    x, y,
                    label,
                    ha='center',
                    va='center',
                    fontsize=7,
                    fontweight='bold',
                    color='white',
                    path_effects=[
                        path_effects.Stroke(linewidth=1.2, foreground='black'),
                        path_effects.Normal()
                    ],
                    zorder=4
                )
        else:
            ax_inset.scatter([x], [y], s=120, c='darkgrey', marker='D', edgecolors='black', linewidths=0.6, zorder=2)

    ax_inset.set_xticks([])
    ax_inset.set_yticks([])
    ax_inset.set_xlim(-0.08, 1.08)
    ax_inset.set_ylim(-0.08, 1.08)
    ax_inset.set_title('Clone Tree (Distance)', fontsize=8)


def run_plot_phylogeography(df_all, trees, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    if df_all is None or df_all.empty or 'x' not in df_all.columns:
        return

    for model_key, G in trees.items():
        pred_col = model_key
        display_name = model_key.replace('_Preds', '')

        if pred_col not in df_all.columns:
            continue

        logging.info(f"Plotting Spatial Phylogeography for {display_name}...")

        plot_data = df_all.dropna(subset=[pred_col, 'x', 'y']).copy()
        plot_data['clean_pred'] = plot_data[pred_col].apply(clean_label)
        plot_data = ensure_slice_id(plot_data)

        mapping = {
            n: clean_label(n)
            for n in G.nodes()
            if G.nodes[n].get('is_observed', False)
        }
        nx.relabel_nodes(G, mapping, copy=False)

        root_name = get_root_node(G)

        clones = sorted([c for c in plot_data['clean_pred'].unique() if c != 'nan'])
        clone_to_label = {clone: str(i + 1) for i, clone in enumerate(clones)}
        plot_data['display_pred'] = plot_data['clean_pred'].map(clone_to_label)

        palette_colors = sns.color_palette('Set1', n_colors=max(9, len(clones)))
        color_dict = {clone_to_label[c]: col for c, col in zip(clones, palette_colors)}

        hue_order = sorted(color_dict.keys(), key=lambda x: int(x) if str(x).isdigit() else str(x))

        pos = compute_spatial_anchors(
            G,
            plot_data,
            pred_col,
            root_name=root_name,
            internal_offset_frac=0.03,
            blend_to_children=0.3,
            internal_lift_frac=0.12,
            internal_left_frac=0.08,
            depth_left_scale=0.3
        )

        fig, ax = plt.subplots(figsize=(14, 12))

        # Adapt spot size based on spatial grid spacing (stable across slice subsets)
        spot_size = estimate_spot_size(plot_data)

        sns.scatterplot(
            data=plot_data,
            x='x', y='y',
            hue='display_pred',
            hue_order=hue_order,
            palette=color_dict,
            s=spot_size,
            edgecolor='none',
            alpha=1.0,
            ax=ax,
            legend=True
        )

        x_span = plot_data['x'].max() - plot_data['x'].min()
        jump_threshold = x_span * 0.3

        for u, v in G.edges():
            if u not in pos or v not in pos:
                continue
            ux, uy = pos[u]
            vx, vy = pos[v]

            if np.isnan(ux) or np.isnan(vx):
                continue

            distance = np.sqrt((ux - vx) ** 2 + (uy - vy) ** 2)
            linestyle = "dashed" if distance > jump_threshold else "solid"
            alpha_edge = 0.8 if distance > jump_threshold else 1.0

            arrow = patches.FancyArrowPatch(
                (ux, uy), (vx, vy),
                connectionstyle="arc3,rad=0.15",
                arrowstyle="->",
                color='black',
                linewidth=3.0,
                linestyle=linestyle,
                alpha=alpha_edge,
                mutation_scale=25,
                zorder=2
            )
            ax.add_patch(arrow)

        # Scale node size by the total number of clones.
        clone_counts = plot_data['clean_pred'].value_counts()
        min_count = clone_counts.min() if len(clone_counts) > 0 else 1
        max_count = clone_counts.max() if len(clone_counts) > 0 else 1

        for node in G.nodes():
            if node not in pos:
                continue

            nx_coord, ny_coord = pos[node]
            if np.isnan(nx_coord):
                continue

            is_observed = G.nodes[node].get('is_observed', G.nodes[node].get('is_leaf', False))
            clean_node_name = clean_label(node)

            if is_observed:
                display_label = clone_to_label.get(clean_node_name, clean_node_name)
                node_color = color_dict.get(display_label, 'red')
                marker = 'o'
                count_val = clone_counts.get(clean_node_name, 1)
                size = scale_size(count_val, min_count, max_count, 160, 420)
            else:
                node_color = 'darkgrey'
                marker, size = 'D', 120
                display_label = None

            ax.scatter(
                nx_coord, ny_coord,
                c=[node_color],
                edgecolors='black',
                linewidths=2.0,
                s=size,
                marker=marker,
                zorder=3
            )

            if is_observed and display_label is not None:
                ax.text(
                    nx_coord, ny_coord,
                    display_label,
                    ha='center',
                    va='center',
                    fontsize=12,
                    fontweight='bold',
                    color='white',
                    path_effects=[
                        path_effects.Stroke(linewidth=2.0, foreground='black'),
                        path_effects.Normal()
                    ],
                    zorder=4
                )

        ax.set_title(f"Spatial Phylogeography: {display_name}", fontsize=20, fontweight='bold')
        ax.set_xlabel("Spatial X", fontsize=14)
        ax.set_ylabel("Spatial Y", fontsize=14)

        # Draw inset phylogram with branch lengths
        _draw_tree_inset(ax, G, root_name, clone_to_label, color_dict, clone_counts, min_count, max_count)

        ax.axis('equal')

        # Auto-center y-limits based on current data range
        y_values = []
        if 'y' in plot_data.columns:
            y_values.extend(plot_data['y'].dropna().tolist())
        for node, (nx_coord, ny_coord) in pos.items():
            if np.isfinite(ny_coord):
                y_values.append(ny_coord)
        if len(y_values) > 0:
            y_min, y_max = float(min(y_values)), float(max(y_values))
            span = max(y_max - y_min, 1.0)
            center = 0.5 * (y_min + y_max)
            pad = span * 0.08
            ax.set_ylim(center - span / 2.0 - pad, center + span / 2.0 + pad)

        handles, labels = ax.get_legend_handles_labels()
        clean_handles, clean_labels = [], []
        for h, l in zip(handles, labels):
            if str(l) != 'display_pred':
                clean_handles.append(h)
                clean_labels.append(l)

        if clean_handles:
            ax.legend(clean_handles, clean_labels, title="Subclone",
                      bbox_to_anchor=(1.02, 1), loc='upper left', markerscale=3, fontsize=12)

        plt.tight_layout()
        out_path = os.path.join(out_dir, f"{display_name}_spatial_phylogeography.png")
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close()
