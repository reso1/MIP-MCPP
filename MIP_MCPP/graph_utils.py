import random

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.colors import rgb2hex

from MIP_MCPP.misc import uv_sorted


def nx_graph_read(filepath: str) -> nx.Graph:
    f = open(filepath, 'r')
    lines = f.readlines()
    node_items = lines[0].split('#')

    G = nx.Graph()
    for i, node_str in enumerate(node_items):
        pos_x, pos_y = node_str.strip().split(', ')
        pos = (float(pos_x[1:]), float(pos_y[:-1]))
        G.add_node(i, pos=pos)

    edge_items = lines[1:]
    for edge_str in edge_items:
        s, t, w = edge_str.strip().split('#')
        G.add_edge(int(s), int(t), weight=float(w))

    return G


def graph_plot(
    G: nx.Graph, ax,
    graph_scale=1.0,
    marker_shape='.',
    color_rgba=(0.5, 0.5, 0.5, 1),
    with_node_labels=True,
    with_edge_labels=True
) -> None:

    color_hex = rgb2hex(color_rgba)

    pos_dict = nx.get_node_attributes(G, 'pos')

    if with_node_labels:
        nx.draw(
            G, pos_dict, ax=ax, with_labels=True,
            node_color="w", node_shape='o', edgecolors=color_hex, linewidths=2,
            edge_color=color_hex, width=2
        )
    else:
        nx.draw(
            G, pos_dict, ax=ax, with_labels=with_node_labels,
            node_color=color_hex, node_shape=marker_shape, edgecolors=color_hex, linewidths=1,
            edge_color=color_hex, width=1, node_size=24 * graph_scale
        )
        # nx.draw(
        #     G, pos_dict, ax=ax, with_labels=with_node_labels,
        #     node_color = color_hex, node_shape = '.', edgecolors = 'k', linewidths = 0.5,
        #     edge_color = 'k', width = 1, node_size = 24 * graph_scale
        # )

    if with_edge_labels:
        edge_labels = {
            (u, v): f"{G[u][v]['weight']:.1f}" for u, v in list(G.edges())}
        nx.draw_networkx_edge_labels(
            G, pos=pos_dict, ax=ax, edge_labels=edge_labels)


def create_grid_graph(n_rows: int, n_cols: int) -> nx.Graph:
    # Create an empty graph
    G = nx.Graph()

    def v_index(i, j): return i * n_cols + j

    # Add nodes with random negative weights and positions
    for i in range(n_rows):
        for j in range(n_cols):
            terrain_weight = random.uniform(1, 4)
            G.add_node(v_index(i, j), pos=(i, j),
                       terrain_weight=round(terrain_weight, 3))

    # Add edges with random positive weights
    for i in range(n_rows):
        for j in range(n_cols):
            u = v_index(i, j)
            if i > 0:
                v = v_index(i-1, j)
                G.add_edge(u, v, weight=(
                    G.nodes[u]["terrain_weight"]+G.nodes[v]["terrain_weight"])/2)
            if j > 0:
                v = v_index(i, j-1)
                G.add_edge(u, v, weight=(
                    G.nodes[u]["terrain_weight"]+G.nodes[v]["terrain_weight"])/2)

    return G


def remove_cycles(G: nx.Graph, Ti_edges):
    Ti_edges = set([uv_sorted((u, v)) for u, v in Ti_edges])
    Ti = G.edge_subgraph(Ti_edges).copy()

    while True:
        cycles = nx.cycle_basis(Ti)
        if cycles == []:
            break

        C = cycles[0]
        C_edges = [(C[i], C[(i+1) % len(C)]) for i in range(len(C))]
        u, v = sorted(
            C_edges, key=lambda e: Ti[e[0]][e[1]]["weight"], reverse=True)[0]
        Ti.remove_edge(u, v)
        Ti_edges.remove(uv_sorted((u, v)))

    return Ti, Ti_edges


if __name__ == '__main__':
    import os

    import matplotlib.pyplot as plt
    G = nx_graph_read(os.path.join(
        'graph_data', 'GRID_20x20_UNWEIGHTED_FREE.graph'))
    # find_all_cycles(G)
    C = nx.find_cycle(G)
    print(C)
    fig, ax = plt.subplots()
    graph_plot(G, ax, with_edge_labels=False)
    plt.show()
