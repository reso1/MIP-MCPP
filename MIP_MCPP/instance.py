from __future__ import annotations

import os
import pickle
import random
from itertools import product
from typing import Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from MIP_MCPP.graph_utils import create_grid_graph, graph_plot
from MIP_MCPP.misc import PLT_SHAPES, colormap, uv_sorted


class Instance:

    """ Problem Instance Class """
    DEFAULT_INSTANCE_DIR = os.path.join("data", "instances")
    DEFAULT_SOLUTION_DIR = os.path.join("data", "solutions")
    DEFAULT_IMAGE_DIR = os.path.join("data", "figs")

    def __init__(self, G: nx.Graph, R: list, name: str) -> None:
        self.k = len(R)
        self.n, self.m = G.number_of_nodes(), G.number_of_edges()
        self.name = name

        grid_spec = self.name.split("-")[0]
        self.width, self.height = [int(v) for v in grid_spec.split("x")]

        self.G = G
        self.R = R
        self.I = list(range(self.k))
        self.V = list(range(self.n))
        self.E = list(range(self.m))

        w_dict = nx.get_edge_attributes(G, "weight")
        self.uv2e = {uv_sorted(uv): e for e, uv in enumerate(w_dict.keys())}
        self.e2uv = [uv_sorted(uv) for uv in G.edges()]
        self.w = np.array(list(w_dict.values()))

    def draw_instance(self, ax, marker_scale=1.0, root=True, obstacle=True) -> None:
        pos = nx.get_node_attributes(self.G, "pos")
        pos2v = {p: v for v, p in pos.items()}
        x_lb, y_lb, x_ub, y_ub = self.bounds

        graph_plot(
            self.G, ax,
            graph_scale=marker_scale,
            color_rgba=(0.5, 0.5, 0.5, 1),
            with_node_labels=False, with_edge_labels=False)

        if root:
            for r in self.R:
                ax.scatter(pos[r][0], pos[r][1], 64 *
                           marker_scale, marker='o', color='r')

        if obstacle:
            for p in product(range(x_lb, x_ub+1), range(y_lb, y_ub+1)):
                if p not in pos2v:
                    ax.scatter(p[0], p[1], marker='s',
                               s=200*marker_scale, color='k')

        ax.axis("equal")
        ax.set_xlim(x_lb-1, x_ub+1)
        ax.set_ylim(y_lb-1, y_ub+1)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.axis("off")
        # ax.xaxis.set_ticks(np.linspace(x_lb, x_ub, (x_ub-x_lb+1)//2))
        # ax.yaxis.set_ticks(np.linspace(y_lb, y_ub, (y_ub-y_lb+1)//2))
        # ax.xaxis.set_ticklabels([])
        # ax.yaxis.set_ticklabels([])
        # ax.grid(True, color='gray', linestyle='--', linewidth=0.25)

    def save_instance_image(self, filename, dir: str = DEFAULT_IMAGE_DIR) -> None:
        path = os.path.join(
            dir, self.name + '.png' if filename is None else filename)
        fig, ax = plt.subplots()
        self.draw_instance(ax, obstacle=False)
        fig.set_size_inches(20, 20)
        fig.tight_layout()
        fig.savefig(path, dpi=200)

    def draw_solution(self, edges, ax, graph_scale=3.0) -> None:
        pos_dict = nx.get_node_attributes(self.G, 'pos')
        arrowprops = dict(facecolor='black', shrink=0.05,
                          width=0.1, headwidth=0.1)
        for r in self.R:
            xy = (pos_dict[r][0]+0.1, pos_dict[r][1]+0.1)
            xytext = (pos_dict[r][0]+0.25, pos_dict[r][1]+0.25)
            # ax.annotate(text='R', xy=xy, xytext=xytext, arrowprops=arrowprops)

        cmap = colormap("tab20")
        for i in self.I:
            Ti_e = self.G.edge_subgraph(edges[i]) if edges[i] != [
            ] else self.G.subgraph(self.R[i])
            graph_plot(
                Ti_e, ax,
                graph_scale=graph_scale,
                marker_shape=PLT_SHAPES[i % self.k],
                color_rgba=cmap(i/self.k),
                with_node_labels=False, with_edge_labels=False
            )

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        pos = nx.get_node_attributes(self.G, "pos")
        x_lb = min([p[0] for p in list(pos.values())])
        y_lb = min([p[1] for p in list(pos.values())])
        x_ub = max([p[0] for p in list(pos.values())])
        y_ub = max([p[1] for p in list(pos.values())])
        return x_lb, y_lb, x_ub, y_ub

    @staticmethod
    def varname(k: int, pool_size: int) -> str:
        return f"lambda{k}_{pool_size}"

    @staticmethod
    def read(filename: str, dir: str = DEFAULT_INSTANCE_DIR) -> Instance:
        with open(os.path.join(dir, filename), 'rb') as f:
            istc = pickle.load(f)
        return istc

    def write(self, filename: str = None, dir: str = DEFAULT_INSTANCE_DIR) -> None:
        path = os.path.join(
            dir, self.name+'.istc' if filename is None else filename)
        with open(path, 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        print(f"Wrote instance to {path}")

    @staticmethod
    def write_solution(edges: list, filename: str, dir: str = DEFAULT_SOLUTION_DIR) -> None:
        with open(os.path.join(dir, filename), 'wb') as f:
            pickle.dump(edges, f, pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def read_solution(filename: str, dir: str = DEFAULT_SOLUTION_DIR) -> list:
        with open(os.path.join(dir, filename+'.solu'), 'rb') as f:
            edges = pickle.load(f)
        return edges

    @staticmethod
    def create_random_free(name, R=None) -> Instance:
        items = name.split("-")
        width, height = [int(v) for v in items[0].split("x")]
        k = int(items[-1][1:])

        G = create_grid_graph(width, height)

        if R is None:
            n = G.number_of_nodes()
            R = [random.randint(0, n-1) for _ in range(k)]

        return Instance(G, R, name)

    @staticmethod
    def create_from_binary_map(filename: str) -> Instance:
        OBS = (0, 0, 0)
        ROOT = (1, 0, 0)
        map = plt.imread(filename)
        width, height, _ = map.shape
        def cvt_idx(x, y): return (y, width-x)

        G, R = nx.Graph(), []
        num_nodes = 0
        pos2v = {}
        for x in range(width):
            for y in range(height):
                if tuple(map[x][y]) != OBS:
                    G.add_node(num_nodes, pos=cvt_idx(x, y))
                    pos2v[cvt_idx(x, y)] = num_nodes
                    if tuple(map[x][y]) == ROOT:
                        R.append(num_nodes)
                    num_nodes += 1

        for x in range(width):
            for y in range(height):
                for inc_x, inc_y in [(-1, 0), (1, 0), (0, 1), (0, -1)]:
                    xc, yc = x + inc_x, y + inc_y
                    if cvt_idx(x, y) in pos2v and cvt_idx(xc, yc) in pos2v:
                        G.add_edge(pos2v[cvt_idx(x, y)],
                                   pos2v[cvt_idx(xc, yc)], weight=1)

        istc_name = filename.split("/")[-1]
        return Instance(G, R, istc_name.split(".")[0])

    def __str__(self) -> str:
        return ",\t".join([
            self.name,
            f"% of obs.={1 - self.n/(self.width*self.height):.3f}",
            f"# of verts={self.n}",
            f"# of edges={self.m}",
            f"# of vars.={(self.n + 3*self.m)*self.k + 1}"
        ])

    def draw_covering_nodes(self, ax, graph_scale=1.0) -> None:
        pos = nx.get_node_attributes(self.G, "pos")
        for v in self.G.nodes():
            for xc, yc in [(0.25, -0.25), (0.25, 0.25), (-0.25, 0.25), (-0.25, -0.25)]:
                _nx, _ny = pos[v][0] + xc, pos[v][1] + yc
                ax.plot(_nx, _ny, "o", mfc="w", mec="k",
                        ms=0.5*graph_scale, alpha=0.5)
