import time
from typing import Tuple, List

import matplotlib
import matplotlib.pyplot as plt
import networkx as nx

from MSTC_Star.mcpp.mfc_planner import MFCPlanner
from MSTC_Star.mcpp.mstc_star_planner import MSTCStarPlanner
from MSTC_Star.mcpp.rtc_planner import RTCPlanner
from MSTC_Star.mcpp.stc_planner import STCPlanner
from MSTC_Star.utils.nx_graph import mst, navigate
from MSTC_Star.utils.robot import Robot
from MIP_MCPP.graph_utils import graph_plot
from MIP_MCPP.instance import Instance


def simulate(
    name: str,
    planner: STCPlanner,
    paths: list,
    weights: list,
    scale: float,
    dt: float,
    obs_graph: nx.Graph,
    is_write: bool = False,
    is_show: bool = False
) -> None:

    k, R = planner.k, planner.R
    color = ['r', 'm', 'b', 'k', 'c', 'g']
    fig = plt.figure()
    fig.set_size_inches(8*scale, 8*scale)
    fig.tight_layout()
    ax = plt.axes()
    # ax.margins(x=0.15, y=0.15)
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    plt.grid(True)
    plt.gcf().canvas.mpl_connect(
        'key_release_event',
        lambda event: [exit(0) if event.key == 'escape' else None])

    robots = [Robot(paths[i], planner.H) for i in range(k)]
    t_finish = [robots[i].T[-1] for i in range(k)]
    t_max = max(t_finish)

    if not is_write and not is_show:
        print(f'Final Max Weights: {max(weights)}')
        return

    lines, markers, texts = [None]*k, [None]*k, [None]*(k+1)
    xs_vec, ys_vec = [None]*k, [None]*k

    def init():
        plt.title(f'{name} (Max Weights={max(weights): .2f})',
                  fontdict={'size': 12*scale})
        texts[-1] = ax.text(
            1, 1, '', va='top', ha='right', transform=ax.transAxes,
            font={'size': 8*scale})

        # MST of spanning graph
        M = mst(planner.G)
        for s, t in M.edges():
            x1, y1 = s
            x2, y2 = t
            ax.plot([x1, x2], [y1, y2], 'ok', mfc='r')
        # covering nodes
        rho = planner.generate_cover_trajectory(R[0], mst(planner.G))
        for cn_x, cn_y in rho:
            ax.plot(cn_x, cn_y, 'o', mec='k', mfc='w', ms=5)
        # obstacle graph
        for s, t in obs_graph.edges():
            x1, y1 = s
            x2, y2 = t
            ax.plot([x1, x2], [y1, y2], '-xk', ms=10, mew=3)

        for i in range(k):
            c = color[i % len(color)]
            line, = ax.plot([], [], '-'+c, alpha=0.35, lw=8)
            marker, = ax.plot([], [], 'o'+c, ms=8)
            # changable texts
            texts[i] = ax.text(
                1, 0.975-i*0.025, '', va='top', ha='right',
                transform=ax.transAxes, font={'size': 8})
            # trajectories and robots
            lines[i], markers[i] = line, marker
            xs_vec[i], ys_vec[i] = zip(*paths[i])
            # depots
            ax.plot(R[i][0], R[i][1], '*k', mfc=c, ms=10)
            # ax.text(R[i][0]+0.1, R[i][1]+0.1, f'R{i}')

        return lines + markers + texts

    # record remaining uncovered nodes
    uncovered = set()
    direction = ['SE', 'NE', 'NW', 'SW']
    for node in planner.G.nodes:
        for sn in [planner.__get_subnode_coords__(node, d) for d in direction]:
            uncovered.add(sn)

    def animate(ti):
        ts = ti * dt
        for i in range(k):
            last_coord_idx, cur_state = robots[i].get_cur_state(ts)
            xs = xs_vec[i][:last_coord_idx+1] + (cur_state.x, )
            ys = ys_vec[i][:last_coord_idx+1] + (cur_state.y, )
            # texts[i].set_text(f'R{i}: ')
            lines[i].set_data(xs, ys)
            markers[i].set_data(cur_state.x, cur_state.y)
            node = (xs_vec[i][last_coord_idx], ys_vec[i][last_coord_idx])
            if node in uncovered:
                uncovered.remove(node)
        texts[-1].set_text(f'T[s]={ts: .2f}, # of uncovered={len(uncovered)}')

        return lines + markers + texts

    anim = matplotlib.animation.FuncAnimation(
        fig, animate, int(t_max/dt), init, interval=5,
        blit=True, repeat=False, cache_frame_data=False)

    if is_write:
        anim.save(f'data/sim_records/{name}.mp4', fps=300, dpi=200,
                  progress_callback=lambda i, n: print(f'{name}: saving frame {i}/{n}'))

    if is_show:
        plt.show()


def mfc_plan(istc: Instance, ax=None, graph_scale=3.0
             ) -> Tuple[MFCPlanner, list, list, float]:

    pos = nx.get_node_attributes(istc.G, "pos")
    tw = nx.get_node_attributes(istc.G, "terrain_weight")
    pos2v = {p: v for v, p in pos.items()}

    G = nx.Graph()
    for u, v in istc.G.edges():
        G.add_edge(pos[u], pos[v], weight=istc.G[u][v]["weight"])

    nx.set_node_attributes(
        G, {pos[v]: w for v, w in tw.items()}, "terrain_weight")
    R = [pos[r] for r in istc.R]
    cap = float("inf")

    planner = MFCPlanner(G, istc.k, R, cap)
    rtc_planner = RTCPlanner(planner.G, planner.R, planner.k)

    ts0 = time.time()
    match_tuple, max_weights, opt_B = rtc_planner.k_tree_cover()
    ts1 = time.time()

    sol_edges = {r:set() for r in R}
    for r, val in match_tuple.items():
        L, S, P = val
        for idx in range(len(P)-1):
            sol_edges[r].add((pos2v[P[idx]], pos2v[P[idx+1]]))
        for u, v in L.edges():
            sol_edges[r].add((pos2v[u], pos2v[v]))
        for u, v in S.edges():
            sol_edges[r].add((pos2v[u], pos2v[v]))

    paths, weights = STC_on_RMMTC_sol(planner, istc, [sol_edges[pos[r]] for r in istc.R])

    for i, t in enumerate(paths):
        tx, ty = zip(*t)
        if ax:
            ax.plot(tx, ty, lw=4*graph_scale, alpha=0.3)

    if ax:
        istc.draw_solution(sol_edges, ax, graph_scale)

    print(f"Planning Time = {ts1-ts0} secs")

    return planner, paths, weights, ts1-ts0


def mstcstar_plan(istc: Instance, ax=None, graph_scale=3.0
                  ) -> Tuple[MSTCStarPlanner, list, list, float]:
    pos = nx.get_node_attributes(istc.G, "pos")
    tw = nx.get_node_attributes(istc.G, "terrain_weight")

    G = nx.Graph()
    for u, v in istc.G.edges():
        G.add_edge(pos[u], pos[v], weight=istc.G[u][v]["weight"])

    nx.set_node_attributes(
        G, {pos[v]: w for v, w in tw.items()}, "terrain_weight")
    R = [pos[r] for r in istc.R]
    cap = float("inf")

    planner = MSTCStarPlanner(G, istc.k, R, cap, True)
    ts0 = time.time()
    plans = planner.allocate()
    ts1 = time.time()

    paths, weights = planner.simulate(plans, False)
    for i, t in enumerate(paths):
        tx, ty = zip(*t)
        if ax:
            ax.plot(tx, ty, lw=4*graph_scale, alpha=0.3)

    if ax:
        M = nx.minimum_spanning_tree(istc.G)
        graph_plot(M, ax, graph_scale=3*graph_scale,
                   with_node_labels=False, with_edge_labels=False)

    print(f"Planning Time = {ts1-ts0} secs")

    return planner, paths, weights, ts1-ts0


class MIPPlanner(STCPlanner):

    def __init__(self, istc: Instance) -> None:
        self.istc = istc

        pos = nx.get_node_attributes(self.istc.G, "pos")
        tw = nx.get_node_attributes(istc.G, "terrain_weight")

        G = nx.Graph()
        for u, v in self.istc.G.edges():
            G.add_edge(pos[u], pos[v], weight=self.istc.G[u][v]["weight"])

        nx.set_node_attributes(
            G, {pos[v]: w for v, w in tw.items()}, "terrain_weight")

        self.R = [pos[r] for r in self.istc.R]
        self.k = istc.k
        self.G = G
        self.H = self.generate_decomposed_graph(G, self.R)
    
    def simulate(self, sol_edges) -> Tuple[List, List]:
        return STC_on_RMMTC_sol(self, self.istc, sol_edges)


def mip_plan(istc: Instance, sol_edges, ax=None, graph_scale=3.0
             ) -> Tuple[MIPPlanner, list, list, float]:
    planner = MIPPlanner(istc)
    paths, weights = planner.simulate(sol_edges)

    for i, t in enumerate(paths):
        tx, ty = zip(*t)
        if ax:
            ax.plot(tx, ty, lw=4*graph_scale, alpha=0.3)

    if ax:
        istc.draw_solution(sol_edges, ax, graph_scale)

    return planner, paths, weights


def STC_on_RMMTC_sol(planner:STCPlanner, istc:Instance, sol_edges):
    pos = nx.get_node_attributes(istc.G, "pos")
    traj, weights = [[] for i in istc.I], [0 for _ in istc.I]
    for i, r in enumerate(istc.R):
        Ti = nx.Graph()
        for u, v in sol_edges[i]:
            Ti.add_edge(pos[u], pos[v], weight=istc.G[u][v]["weight"])

        if Ti.number_of_edges() != 0:
            pi = planner.generate_cover_trajectory(pos[r], Ti)
        else:
            pi = [planner.__get_subnode_coords__(pos[r], d) for d in ["SE", "NE", "NW", "SW"]]

        traj[i] = [pos[r]] + pi + [pos[r]]
        weights[i] = planner.__get_travel_weights__(traj[i])

    return traj, weights
