import logging
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.mixture import GaussianMixture
from ...utils.constants import HG38_INFO

def _build_cost_matrix():
    # state indices: 0=loss, 1=neutral, 2=gain
    # Hard irreversible: allow 0->loss/gain, forbid reversals
    inf = np.inf
    cm = np.full((3, 3), inf, dtype=float)
    for s in range(3):
        cm[s, s] = 0.0
    cm[1, 0] = 1.0
    cm[1, 2] = 1.0
    return cm




def discretize_matrix(df_matrix):
    """
    Use GMM to discretize CNV signals into Loss(-1), Neutral(0), Gain(1).
    """
    mat_values = df_matrix.values
    flat_values = mat_values.flatten().reshape(-1, 1)

    if np.var(flat_values) < 1e-6:
        logging.warning("Variance across all clones is near zero. Assuming no CNV events.")
        return pd.DataFrame(0, index=df_matrix.index, columns=df_matrix.columns)

    try:
        gmm = GaussianMixture(n_components=3, covariance_type='full', random_state=42, reg_covar=1e-4)
        gmm.fit(flat_values)
        means = gmm.means_.flatten()
        sorted_indices = np.argsort(means)
        mapping_dict = {
            sorted_indices[0]: -1,
            sorted_indices[1]: 0,
            sorted_indices[2]: 1,
        }
        full_preds = gmm.predict(flat_values)
        discrete_values = np.vectorize(mapping_dict.get)(full_preds)
    except Exception as e:
        logging.error(f"GMM fitting failed: {e}. Falling back to Z-score method.")
        mean_val = np.mean(flat_values)
        std_val = np.std(flat_values)
        discrete_values = np.zeros_like(flat_values)
        discrete_values[flat_values > mean_val + 0.5 * std_val] = 1
        discrete_values[flat_values < mean_val - 0.5 * std_val] = -1

    return pd.DataFrame(
        discrete_values.reshape(mat_values.shape),
        index=df_matrix.index,
        columns=df_matrix.columns,
    )


def _events_from_discrete(discrete_matrix):
    events = {}
    for clone in discrete_matrix.columns:
        series = discrete_matrix[clone]
        nonzero = series[series != 0]
        if len(nonzero) == 0:
            events[clone] = set()
        else:
            keys = nonzero.index.astype(str) + ":" + nonzero.astype(int).astype(str)
            events[clone] = set(keys.tolist())
    return events


def _manhattan_discrete(discrete_matrix, a, b):
    return float(np.abs(discrete_matrix[a].values - discrete_matrix[b].values).sum())


def _subset_violation_ratio(parent_events, child_events):
    if not parent_events:
        return 0.0
    missing = len(parent_events - child_events)
    return missing / max(len(parent_events), 1)


def _ancestor_name_from_leaves(leaf_names):
    cleaned = [str(x) for x in leaf_names if pd.notna(x)]
    cleaned = [c.replace('clone', '').replace('Clone', '') for c in cleaned]
    cleaned = sorted(cleaned)
    return "ancestor_" + "_".join(cleaned)


def _build_subset_tree(discrete_matrix, violation_tolerance=0.05):
    clones = discrete_matrix.columns.tolist()
    G = nx.DiGraph()

    if len(clones) == 0:
        return G

    if violation_tolerance is None:
        violation_tolerance = 0.05
    else:
        violation_tolerance = float(violation_tolerance)

    for c in clones:
        G.add_node(c, is_observed=True)

    events = _events_from_discrete(discrete_matrix)
    event_counts = {c: len(events[c]) for c in clones}
    sorted_clones = sorted(clones, key=lambda c: (event_counts[c], str(c)))

    parent_map = {}

    for idx, clone in enumerate(sorted_clones):
        candidates = []
        for cand in sorted_clones[:idx]:
            if events[cand] == events[clone]:
                continue
            violation = _subset_violation_ratio(events[cand], events[clone])
            if violation <= violation_tolerance:
                dist = _manhattan_discrete(discrete_matrix, cand, clone)
                candidates.append((violation, -event_counts[cand], dist, cand))

        if candidates:
            candidates.sort(key=lambda x: (x[0], x[1], x[2], str(x[3])))
            parent = candidates[0][3]
            parent_map[clone] = parent

    for child, parent in parent_map.items():
        w = _manhattan_discrete(discrete_matrix, parent, child)
        G.add_edge(parent, child, weight=w)

    roots = [c for c in clones if c not in parent_map]
    if len(roots) > 1:
        root_name = _ancestor_name_from_leaves(clones)
        if root_name not in G:
            G.add_node(root_name, is_observed=False)
        for r in roots:
            w = float(np.abs(discrete_matrix[r].values).sum())
            G.add_edge(root_name, r, weight=w)

    for node in G.nodes():
        if 'is_observed' not in G.nodes[node]:
            G.nodes[node]['is_observed'] = False
        G.nodes[node]['is_leaf'] = G.out_degree(node) == 0

    return G


def _encode_states(discrete_matrix):
    states = {}
    for clone in discrete_matrix.columns:
        arr = discrete_matrix[clone].values
        states[clone] = np.where(arr < 0, 0, np.where(arr > 0, 2, 1)).astype(np.int8)
    return states



def _parsimony_score(adj, leaf_states, root, cost_matrix):
    parent = {root: None}
    order = [root]
    for node in order:
        for nb in adj[node]:
            if nb == parent[node]:
                continue
            parent[nb] = node
            order.append(nb)

    costs = {}
    n_bins = len(next(iter(leaf_states.values())))
    for node in reversed(order):
        if node in leaf_states:
            cost = np.full((3, n_bins), np.inf, dtype=float)
            idx = leaf_states[node]
            cost[idx, np.arange(n_bins)] = 0.0
        else:
            children = [ch for ch in adj[node] if ch != parent[node]]
            if not children:
                cost = np.zeros((3, n_bins), dtype=float)
            else:
                cost = np.zeros((3, n_bins), dtype=float)
                for child in children:
                    child_cost = costs[child]
                    mincost = np.empty((3, n_bins), dtype=float)
                    for s in range(3):
                        mincost[s] = np.min(child_cost + cost_matrix[s][:, None], axis=0)
                    cost += mincost
        costs[node] = cost

    root_cost = costs[root]
    # Hard irreversible: force neutral root state (index 1)
    return float(np.sum(root_cost[1]))



def _clone_adj(adj):
    return {k: set(v) for k, v in adj.items()}


def _edge_list(adj):
    edges = []
    seen = set()
    for u in adj:
        for v in adj[u]:
            if (v, u) in seen:
                continue
            edges.append((u, v))
            seen.add((u, v))
    return edges


def _hamming_distance(a, b):
    return int(np.sum(a != b))


def _build_parsimony_tree(discrete_matrix):
    clones = discrete_matrix.columns.tolist()
    adj = {}

    if len(clones) == 0:
        return adj
    if len(clones) == 1:
        adj[clones[0]] = set()
        return adj

    leaf_states = _encode_states(discrete_matrix)
    cost_matrix = _build_cost_matrix()

    # choose initial pair with minimal Hamming distance
    best_pair = (clones[0], clones[1])
    best_dist = None
    for i in range(len(clones)):
        for j in range(i + 1, len(clones)):
            d = _hamming_distance(leaf_states[clones[i]], leaf_states[clones[j]])
            if best_dist is None or d < best_dist:
                best_dist = d
                best_pair = (clones[i], clones[j])

    a, b = best_pair
    adj[a] = {b}
    adj[b] = {a}

    remaining = [c for c in clones if c not in best_pair]
    internal_idx = 1
    root_for_score = clones[0]

    for clone in remaining:
        best_score = None
        best_edge = None

        for u, v in _edge_list(adj):
            temp_adj = _clone_adj(adj)
            temp_internal = f"ancestor_tmp_{internal_idx}"

            temp_adj[u].remove(v)
            temp_adj[v].remove(u)
            temp_adj[temp_internal] = {u, v, clone}
            temp_adj[clone] = {temp_internal}
            temp_adj[u].add(temp_internal)
            temp_adj[v].add(temp_internal)

            score = _parsimony_score(temp_adj, leaf_states, root_for_score, cost_matrix)
            if best_score is None or score < best_score:
                best_score = score
                best_edge = (u, v)

        internal_name = f"ancestor_{internal_idx}"
        internal_idx += 1
        u, v = best_edge
        adj[u].remove(v)
        adj[v].remove(u)
        adj[internal_name] = {u, v, clone}
        adj[clone] = {internal_name}
        adj[u].add(internal_name)
        adj[v].add(internal_name)

    return adj


def _select_root_clone(df_matrix, root_clone=None):
    if root_clone is not None and root_clone in df_matrix.columns:
        return root_clone
    scores = df_matrix.abs().mean(axis=0)
    return scores.idxmin()


def _orient_tree(adj, root, observed, discrete_matrix):
    if root not in adj:
        root = next(iter(observed))

    G = nx.DiGraph()
    for node in adj:
        G.add_node(node, is_observed=node in observed)

    parent = {root: None}
    queue = [root]
    while queue:
        node = queue.pop(0)
        for nb in adj[node]:
            if nb == parent[node]:
                continue
            parent[nb] = node
            queue.append(nb)
            if node in observed and nb in observed:
                w = _manhattan_discrete(discrete_matrix, node, nb)
            else:
                w = 1.0
            G.add_edge(node, nb, weight=w)

    for node in G.nodes():
        G.nodes[node]['is_leaf'] = G.out_degree(node) == 0

    return G


def _normalize_chr(chrom):
    s = str(chrom)
    if s.lower().startswith('chr'):
        s = s[3:]
    return s


def _bin_meta_from_long(df_long):
    if df_long is None or df_long.empty:
        return {}
    cols = ['BinID', 'Chromosome', 'Start', 'End']
    for c in cols:
        if c not in df_long.columns:
            return {}
    meta = df_long[cols].drop_duplicates('BinID')
    out = {}
    for _, row in meta.iterrows():
        out[str(row['BinID'])] = (str(row['Chromosome']), float(row['Start']), float(row['End']))
    return out


def _assign_arm(chrom, start, end):
    chrom_norm = _normalize_chr(chrom)
    info = HG38_INFO.get(chrom_norm)
    if info is None:
        return f"chr{chrom_norm}", None
    cen = info['cen']
    mid = 0.5 * (float(start) + float(end))
    arm = 'p' if mid < cen else 'q'
    return f"chr{chrom_norm}", arm


def _fitch_state_sets(G, root, leaf_states):
    parent = {root: None}
    order = [root]
    for node in order:
        for child in G.successors(node):
            if child == parent.get(node):
                continue
            parent[child] = node
            order.append(child)

    state_sets = {}
    for node in reversed(order):
        if node in leaf_states:
            state_sets[node] = leaf_states[node]
            continue
        children = [state_sets[ch] for ch in G.successors(node)]
        if not children:
            state_sets[node] = np.ones_like(next(iter(leaf_states.values()))) * 2
            continue
        cur = children[0].copy()
        for cs in children[1:]:
            inter = cur & cs
            empty = inter == 0
            cur = np.where(empty, cur | cs, inter)
        state_sets[node] = cur
    return state_sets


def _fitch_assign_states(G, root, state_sets):
    assigned = {}
    queue = [root]
    assigned[root] = _choose_state(state_sets[root], parent_state=None)
    while queue:
        node = queue.pop(0)
        for child in G.successors(node):
            assigned[child] = _choose_state(state_sets[child], parent_state=assigned[node])
            queue.append(child)
    return assigned


def _choose_state(state_set, parent_state=None):
    # state_set is an array of bitmasks (1=loss,2=neutral,4=gain)
    state_set = state_set.astype(np.int8)
    if parent_state is not None:
        parent_state = parent_state.astype(np.int8)
        use_parent = (state_set & parent_state) != 0
    else:
        parent_state = np.zeros_like(state_set, dtype=np.int8)
        use_parent = np.zeros_like(state_set, dtype=bool)

    chosen = np.where(use_parent, parent_state, 0)

    neutral_allowed = (state_set & 2) != 0
    chosen = np.where((chosen == 0) & neutral_allowed, 2, chosen)

    loss_allowed = (state_set & 1) != 0
    gain_allowed = (state_set & 4) != 0
    chosen = np.where((chosen == 0) & loss_allowed, 1, chosen)
    chosen = np.where((chosen == 0) & gain_allowed, 4, chosen)

    chosen = np.where(chosen == 0, 2, chosen)
    return chosen


def _bitmask_to_value(mask):
    if mask == 1:
        return -1
    if mask == 4:
        return 1
    return 0


def _edge_event_labels(G, discrete_matrix, bin_meta, top_k=2):
    leaf_states = {}
    for node in G.nodes():
        if node in discrete_matrix.columns:
            arr = discrete_matrix[node].values
            leaf_states[node] = np.where(arr < 0, 0, np.where(arr > 0, 2, 1)).astype(np.int8)

    root = next((n for n in G.nodes() if G.in_degree(n) == 0), None)
    if root is None or not leaf_states:
        return

    cost_matrix = _build_cost_matrix()
    costs = _sankoff_costs_directed(G, root, leaf_states, cost_matrix)
    assigned = _sankoff_assign_states(G, root, costs, cost_matrix)

    bin_ids = discrete_matrix.index.astype(str).tolist()
    for u, v in G.edges():
        if u not in assigned or v not in assigned:
            continue
        parent_idx = assigned[u]
        child_idx = assigned[v]
        if parent_idx.shape != child_idx.shape:
            continue

        counts = {}
        for idx, bin_id in enumerate(bin_ids):
            meta = bin_meta.get(bin_id)
            if meta is None:
                continue
            chrom, start, end = meta
            chr_name, arm = _assign_arm(chrom, start, end)
            arm_suffix = arm if arm is not None else ''

            p_val = int(parent_idx[idx]) - 1
            c_val = int(child_idx[idx]) - 1
            if c_val == p_val:
                continue
            direction = 'gain' if c_val > p_val else 'loss'
            key = (chr_name, arm_suffix, direction)
            counts[key] = counts.get(key, 0) + 1

        if not counts:
            G.edges[u, v]['edge_labels'] = []
            continue

        items = sorted(counts.items(), key=lambda x: (-x[1], x[0][0], x[0][1], x[0][2]))
        labels = []
        for (chr_name, arm_suffix, direction), _ in items[:top_k]:
            if arm_suffix:
                labels.append(f"{chr_name}{arm_suffix} {direction}")
            else:
                labels.append(f"{chr_name} {direction}")
        G.edges[u, v]['edge_labels'] = labels



def _sankoff_costs_directed(G, root, leaf_states, cost_matrix):
    order = list(nx.topological_sort(G))
    costs = {}
    n_bins = len(next(iter(leaf_states.values())))
    for node in reversed(order):
        if node in leaf_states:
            cost = np.full((3, n_bins), np.inf, dtype=float)
            idx = leaf_states[node]
            cost[idx, np.arange(n_bins)] = 0.0
        else:
            children = list(G.successors(node))
            if not children:
                cost = np.zeros((3, n_bins), dtype=float)
            else:
                cost = np.zeros((3, n_bins), dtype=float)
                for child in children:
                    child_cost = costs[child]
                    mincost = np.empty((3, n_bins), dtype=float)
                    for s in range(3):
                        mincost[s] = np.min(child_cost + cost_matrix[s][:, None], axis=0)
                    cost += mincost
        costs[node] = cost
    return costs


def _pick_states(total_costs, parent_state=None):
    min_cost = np.min(total_costs, axis=0)
    n = min_cost.shape[0]
    chosen = np.full(n, -1, dtype=np.int8)

    if parent_state is not None:
        parent_cost = total_costs[parent_state, np.arange(n)]
        parent_ok = parent_cost == min_cost
        chosen = np.where(parent_ok, parent_state, chosen)

    neutral_ok = total_costs[1] == min_cost
    chosen = np.where((chosen == -1) & neutral_ok, 1, chosen)

    argmin = np.argmin(total_costs, axis=0).astype(np.int8)
    chosen = np.where(chosen == -1, argmin, chosen)
    return chosen


def _sankoff_assign_states(G, root, costs, cost_matrix):
    assigned = {}
    n_bins = costs[root].shape[1]
    assigned[root] = np.ones(n_bins, dtype=np.int8)  # force neutral at root
    queue = [root]
    while queue:
        node = queue.pop(0)
        for child in G.successors(node):
            child_costs = costs[child]
            parent_state = assigned[node]
            total = np.empty_like(child_costs)
            for t in range(3):
                total[t] = child_costs[t] + cost_matrix[parent_state, t]
            assigned[child] = _pick_states(total, parent_state=parent_state)
            queue.append(child)
    return assigned



def run_tree_inference(aligned_results, tree_mode='parsimony', subset_tolerance=0.05, root_clone=None):
    trees = {}

    for model_name, df_long in aligned_results.items():
        if df_long is None or df_long.empty:
            continue

        logging.info(f"--- Inferring Evolution Tree for {model_name} ---")

        df_matrix = df_long.pivot(index='BinID', columns='Clone_ID', values='CN_Score')
        df_matrix = df_matrix.fillna(0)

        # Debug: how many subclones are present for this model
        raw_clones = sorted(df_long['Clone_ID'].dropna().unique().tolist())
        matrix_clones = sorted(df_matrix.columns.tolist())
        logging.info(
            f"Subclones in long table: {len(raw_clones)}; in matrix: {len(matrix_clones)}."
        )

        logging.info("Discretizing CNV signals to mutation events via GMM...")
        discrete_matrix = discretize_matrix(df_matrix)

        try:
            if tree_mode == 'subset':
                logging.info(
                    f"Building subset-based tree (GMM discretization, tolerance={float(subset_tolerance):.2f})..."
                )
                G = _build_subset_tree(discrete_matrix, violation_tolerance=subset_tolerance)
            else:
                logging.info("Building approximate parsimony tree (stepwise addition + irreversible Sankoff scoring)...")
                adj = _build_parsimony_tree(discrete_matrix)
                root = _select_root_clone(df_matrix, root_clone=root_clone)
                G = _orient_tree(adj, root, set(discrete_matrix.columns), discrete_matrix)

            # Attach edge event labels based on parsimony state reconstruction
            bin_meta = _bin_meta_from_long(df_long)
            if bin_meta:
                _edge_event_labels(G, discrete_matrix, bin_meta, top_k=2)

            trees[model_name] = G
            observed = [n for n in G.nodes() if G.nodes[n].get('is_observed', False)]
            logging.info(
                f"Successfully built tree for {model_name} with {G.number_of_nodes()} nodes and "
                f"{G.number_of_edges()} edges; observed subclones in tree: {len(observed)}."
            )
        except Exception as e:
            logging.error(f"Failed to build tree for {model_name}: {e}")

    return trees
