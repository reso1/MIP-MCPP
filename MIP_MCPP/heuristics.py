import heapq
from itertools import product
from typing import Dict, Tuple

import networkx as nx
import numpy as np

from MIP_MCPP.instance import Instance


def connectivity_check(G: nx.Graph, r: int, verts_to_remove: set) -> set:
    V_H = set(G.nodes()) - verts_to_remove
    components = list(nx.connected_components(G.subgraph(V_H)))

    if len(components) != 1:
        # record the component index of each v in V_H
        v_in_comp, root_comp_ind = {}, -1
        for v in V_H:
            for i, c in enumerate(components):
                if v in c:
                    if v == r:
                        root_comp_ind = i
                    else:
                        v_in_comp[v] = i
                    break

        # for each component c, get the nearset vertex of c to V_H
        nearest_comp_verts = {i: set() for i in range(len(components))}
        min_path_cost = [float("inf") for i in range(len(components))]
        min_path = [[] for i in range(len(components))]
        for v in verts_to_remove:
            for u in G.neighbors(v):
                if u in V_H and u not in components[root_comp_ind]:
                    comp_ind = v_in_comp[u]
                    if u not in nearest_comp_verts[comp_ind]:
                        nearest_comp_verts[comp_ind].add(u)
                        path = nx.shortest_path(G, r, u, "weight")
                        path_cost = sum([G[path[i]][path[i+1]]["weight"]
                                        for i in range(len(path)-1)])
                        if path_cost < min_path_cost[comp_ind]:
                            min_path_cost[comp_ind], min_path[comp_ind] = path_cost, path

        for path in min_path:
            verts_to_remove = verts_to_remove - set(path)

    return verts_to_remove


""" Parabolic Removal Heuristics """


def get_endpoint_ray(x_lb, x_ub, y_lb, y_ub, pos_s, pos_t, pos_dict) -> Tuple[float, float]:

    px_s, py_s = pos_s
    px_t, py_t = pos_t
    dx, dy = px_t - px_s, py_t - py_s
    def k(): return dy / dx
    def b(): return py_t - k() * px_t

    if (dy >= dx >= 0) or (dy > 0 >= dx):
        for y in np.linspace(y_ub, py_t, num=y_ub-py_t+1, endpoint=True):
            x = round((y - b()) / k()) if dx != 0 else px_t
            x = x_lb if x_lb == x + 1 else (x_ub if x_ub == x - 1 else x)
            if (x, y) in pos_dict:
                return (int(x), int(y))
    elif (dy <= dx <= 0) or (dy < 0 <= dx):
        for y in np.linspace(y_lb, py_t, num=py_t-y_lb+1, endpoint=True):
            x = round((y - b()) / k()) if dx != 0 else px_t
            x = x_lb if x_lb == x + 1 else (x_ub if x_ub == x - 1 else x)
            if (x, y) in pos_dict:
                return (int(x), int(y))
    elif (dx > dy >= 0) or (dx > 0 >= dy):
        for x in np.linspace(x_ub, px_t, num=x_ub-px_t+1, endpoint=True):
            y = round(k() * x + b()) if dy != 0 else py_t
            y = y_lb if y_lb == y + 1 else (y_ub if y_ub == y - 1 else y)
            if (x, y) in pos_dict:
                return (int(x), int(y))
    elif (dx < dy <= 0) or (dx < 0 <= dy):
        for x in np.linspace(x_lb, px_t, num=px_t-x_lb+1, endpoint=True):
            y = round(k() * x + b()) if dy != 0 else py_t
            y = y_lb if y_lb == y + 1 else (y_ub if y_ub == y - 1 else y)
            if (x, y) in pos_dict:
                return (int(x), int(y))

    return (px_t, py_t)


def heur_parabolic(istc: Instance, alpha: float) -> Tuple[dict, dict, dict]:
    assert(alpha >= 0)

    V = {r: set() for r in istc.R}
    Pxr = {r: {} for r in istc.R}
    Pyr = {r: {} for r in istc.R}

    pos = nx.get_node_attributes(istc.G, "pos")
    pos2v = {p: v for v, p in pos.items()}
    x_lb, y_lb, x_ub, y_ub = istc.bounds

    def parabola(x, a, b): return a * (x - b[0])**2 + b[1]
    d2r = {r: nx.shortest_path_length(
        istc.G, source=r, weight="weight") for r in istc.R}
    degrees = nx.degree(istc.G)
    B = [v for v in istc.V if degrees[v] < 4]

    def sigmoid(x): return 1 / (1 + np.exp(-x))

    for ri, rj in product(istc.R, istc.R):

        if ri == rj:
            continue

        d1 = nx.shortest_path_length(istc.G, ri, rj)
        cij = sorted(B, key=lambda v: d2r[ri][v] - d2r[rj][v], reverse=True)[0]
        d2 = nx.shortest_path_length(istc.G, rj, cij)

        verts_to_remove = set()
        a = (alpha * sigmoid(d2/d1))**2
        rot = np.arctan2(pos[ri][1]-pos[rj][1], pos[ri]
                         [0]-pos[rj][0]) + np.pi/2
        for p in product(range(x_lb, x_ub+1), range(y_lb, y_ub+1)):
            if p not in pos2v:
                continue

            v = pos2v[p]
            if v not in istc.V:
                continue

            # Rotate the point around rj
            xr = pos[rj][0] + np.cos(-rot)*(p[0] - pos[rj][0]) - \
                np.sin(-rot)*(p[1] - pos[rj][1])
            yr = pos[rj][1] + np.sin(-rot)*(p[0] - pos[rj][0]) + \
                np.cos(-rot)*(p[1] - pos[rj][1])

            if yr >= parabola(xr, a, pos[rj]):
                verts_to_remove.add(v)

        V[ri] = V[ri].union(verts_to_remove)

        px = np.linspace(x_lb-10, x_ub+10, num=1000)
        py = parabola(px, a, pos[rj])
        pxr = pos[rj][0] + np.cos(rot)*(px - pos[rj][0]) - \
            np.sin(rot)*(py - pos[rj][1])
        pyr = pos[rj][1] + np.sin(rot)*(px - pos[rj][0]) + \
            np.cos(rot)*(py - pos[rj][1])
        Pxr[ri][rj] = pxr
        Pyr[ri][rj] = pyr

    for r in istc.R:
        V[r] = connectivity_check(istc.G, r, V[r])

    return V, Pxr, Pyr


""" Subgraph Removal Heuristics """


def get_separation_graph(G: nx.Graph, ri: int, rj: int) -> nx.Graph:

    S = set()
    dist2i = nx.shortest_path_length(G, source=ri, weight="weight")
    dist2j = nx.shortest_path_length(G, source=rj, weight="weight")
    for v in G.nodes():
        if dist2i[v] > dist2j[v]:
            S.add(v)

    return G.subgraph(S)


def farthest_first_search(G: nx.Graph, c: int, n_max_verts: int, cost2r) -> set:
    visited = set([c])

    Q = [(-cost2r(c), c)]
    while Q and len(visited) <= n_max_verts:
        _, v = heapq.heappop(Q)

        if v in G.nodes():
            for u in G.neighbors(v):
                if u not in visited:
                    visited.add(u)
                    heapq.heappush(Q, (-cost2r(u), u))

    return visited


def heur_subgraph(istc: Instance, beta: float = 1.0) -> Tuple[dict, dict, dict]:

    verts_to_remove = {r: set() for r in istc.R}
    Cs = {r: [] for r in istc.R}
    Ss = {r: {} for r in istc.R}

    def sigmoid(x): return 1 / (1 + np.exp(-x))
    d2r = {r: nx.shortest_path_length(
        istc.G, source=r, weight="weight") for r in istc.R}
    degrees = nx.degree(istc.G)
    B = [v for v in istc.V if degrees[v] < 4]

    for ri, rj in product(istc.R, istc.R):
        if ri == rj:
            continue

        Sij = get_separation_graph(istc.G, ri, rj)

        d1 = nx.shortest_path_length(istc.G, ri, rj)
        cij = sorted(B, key=lambda v: d2r[ri][v] - d2r[rj][v], reverse=True)[0]
        d2 = nx.shortest_path_length(istc.G, rj, cij)

        n_verts_lim = round(len(Sij.nodes()) * beta * sigmoid(d1/(1e-6+d2)))

        M_verts = farthest_first_search(
            Sij, cij, n_verts_lim,
            cost2r=lambda v: nx.shortest_path_length(
                istc.G, v, ri, weight="weight")
        )

        Cs[ri].append(cij)
        Ss[ri][rj] = list(Sij.nodes)
        verts_to_remove[ri] = verts_to_remove[ri].union(M_verts)

    for r in istc.R:
        verts_to_remove[r] = connectivity_check(istc.G, r, verts_to_remove[r])

    return verts_to_remove, Cs, Ss
