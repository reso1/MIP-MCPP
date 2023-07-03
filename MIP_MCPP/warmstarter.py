from collections import defaultdict
from typing import Callable, Tuple

import matplotlib.pyplot as plt
import networkx as nx

from MSTC_Star.mcpp.rtc_planner import RTCPlanner
from MIP_MCPP.graph_utils import graph_plot, remove_cycles
from MIP_MCPP.misc import uv_sorted
from MIP_MCPP.model import Model


def gen_flow(G: nx.Graph, debug=False) -> Tuple[dict, dict]:
    if debug:
        fig, ax = plt.subplots()
    k = G.number_of_nodes()
    f_eu = {uv_sorted(uv): 0 for uv in G.edges()}
    f_ev = {uv_sorted(uv): 0 for uv in G.edges()}
    flow_leftover = {v: 1-1/k for v in G.nodes()}

    assert(nx.cycle_basis(G) == [])

    while G.number_of_edges() != 0:

        for e in G.edges():
            u, v = uv_sorted(e)
            e = (u, v)
            if G.degree(u) == 1 and G.degree(v) != 1:
                G.remove_edge(u, v)
                f_eu[e] = round(flow_leftover[u], 12)
                f_ev[e] = round(1 - f_eu[e], 12)
                assert(f_eu[e]+f_ev[e] == 1)
                flow_leftover[u] -= f_eu[e]
                flow_leftover[v] -= f_ev[e]
            elif G.degree(u) != 1 and G.degree(v) == 1:
                G.remove_edge(u, v)
                f_ev[e] = round(flow_leftover[v], 12)
                f_eu[e] = round(1 - f_ev[e], 12)
                assert(f_eu[e]+f_ev[e] == 1)
                flow_leftover[v] -= f_ev[e]
                flow_leftover[u] -= f_eu[e]
            elif G.degree(u) == 1 and G.degree(v) == 1:
                G.remove_edge(u, v)
                f_eu[e] = round(flow_leftover[u], 12)
                f_ev[e] = round(1-flow_leftover[u], 12)
                assert(f_eu[e]+f_ev[e] == 1)
                flow_leftover[u] -= f_eu[e]
                flow_leftover[v] -= f_ev[e]

            if debug:
                ax.cla()
                graph_plot(G, ax, with_edge_labels=False)
                plt.draw()
                plt.pause(0.05)

    return f_eu, f_ev


class WarmStarter:

    @staticmethod
    def warmstart_RTC(model: Model, varfunc: Callable) -> list:
        istc = model.istc
        pos = nx.get_node_attributes(istc.G, "pos")
        pos2v = {p: v for v, p in pos.items()}

        G = nx.Graph()
        for u, v in istc.G.edges():
            G.add_edge(pos[u], pos[v], weight=istc.G[u][v]["weight"])

        pos_rec = defaultdict(list)
        for i in istc.I:
            pos_rec[pos[istc.R[i]]].append(i)

        r2pos = {}
        for p, ii in pos_rec.items():
            r2pos[ii[0]] = p
            for idx in range(1, len(ii)):
                r_dummy = (p[0]+0.1*idx, p[1]+0.1*idx)
                G.add_edge(r_dummy, p, weight=0)
                r2pos[ii[idx]] = r_dummy
                pos2v[r_dummy] = istc.R[ii[idx]]

        R = [r2pos[i] for i in istc.I]
        assert(len(set(R)) == istc.k)
        rtc_planner = RTCPlanner(G, R, len(R))
        match_tuple, _, _ = rtc_planner.k_tree_cover()

        sol_edges = [[] for _ in istc.I]
        sol_verts = [set([]) for _ in istc.I]

        # get corresponding vars. value
        r2i = {r: i for i, r in enumerate(istc.R)}
        for r, val in match_tuple.items():
            r = pos2v[r]
            i = r2i[r]
            L, S, P = val
            L = [(pos2v[u], pos2v[v]) for u, v in L.edges()]
            S = [(pos2v[u], pos2v[v]) for u, v in S.edges()]
            sol_edges[i] = L + S + \
                [(pos2v[P[idx]], pos2v[P[idx+1]]) for idx in range(len(P)-1)]

            if sol_edges[i] == []:
                varfunc(model.y[model.v_ind[i]+model.v2id[i][r]], 1)
                sol_verts[i].add(r)

            for u, v in sol_edges[i]:
                sol_verts[i].add(u)
                sol_verts[i].add(v)
                if u == v:
                    sol_edges[i].remove((u, v))

            Ti, sol_edges[i] = remove_cycles(istc.G, sol_edges[i])

            assert(len(list(nx.connected_components(Ti))) <= 1)

            f_eu, f_ev = gen_flow(Ti)

            for u, v in sol_edges[i]:
                e = istc.uv2e[uv_sorted((u, v))]
                eid = model.e_ind[i] + model.e2id[i][e]
                uid = model.v_ind[i] + model.v2id[i][u]
                vid = model.v_ind[i] + model.v2id[i][v]
                varfunc(model.x[eid], 1)
                varfunc(model.y[uid], 1)
                varfunc(model.y[vid], 1)
                varfunc(model.fu[eid], f_eu[uv_sorted((u, v))])
                varfunc(model.fv[eid], f_ev[uv_sorted((u, v))])

        V = set()
        for v in sol_verts:
            V = V.union(v)

        assert(len(V) == istc.n)
        return sol_edges

    @staticmethod
    def warmstart_MST(model: Model, H, varfunc: Callable) -> list:
        if H is None:
            H = [model.istc.G for _ in model.istc.I]

        sol_edges = [[] for _ in model.istc.I]

        for i in model.istc.I:
            Mi = nx.minimum_spanning_tree(H[i])
            sol_edges[i] = list(Mi.edges())
            f_eu, f_ev = gen_flow(Mi.copy())
            for u, v in Mi.edges():
                e = model.istc.uv2e[uv_sorted((u, v))]
                eid = model.e_ind[i] + model.e2id[i][e]
                uid = model.v_ind[i] + model.v2id[i][u]
                vid = model.v_ind[i] + model.v2id[i][v]
                varfunc(model.x[eid], 1)
                varfunc(model.y[uid], 1)
                varfunc(model.y[vid], 1)
                varfunc(model.fu[eid], f_eu[uv_sorted((u, v))])
                varfunc(model.fv[eid], f_ev[uv_sorted((u, v))])

        return sol_edges

    @staticmethod
    def apply(model: Model, type="MST", H=None) -> Model:
        def _set_var_start(var, val):
            var.setAttr("start", val)

        if type == "RTC":
            init_sol_edges = WarmStarter.warmstart_RTC(model, _set_var_start)
        elif type == "MST":
            init_sol_edges = WarmStarter.warmstart_MST(
                model, H, _set_var_start)

        return model

    @staticmethod
    def check_feasibility(model: Model, type="MST", H=None) -> None:
        def _fix_var(var, val):
            var.setAttr("lb", val)
            var.setAttr("ub", val)

        if type == "RTC":
            WarmStarter.warmstart_RTC(model, _fix_var)
        elif type == "MST":
            WarmStarter.warmstart_MST(model, H, _fix_var)

        try:
            model.model.computeIIS()
            model.model.write("model.ilp")
        except Exception as e:
            print(e)
