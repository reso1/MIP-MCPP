import argparse
import os

import networkx as nx

from MIP_MCPP.instance import Instance
from MIP_MCPP.mcpp_planner import mfc_plan, mip_plan, mstcstar_plan, simulate


def plan(
    istc_name: str,
    method: str,
    istc_sol_name: str,
    scale: float,
    dt: float,
    is_write: bool
) -> None:

    istc = Instance.read(f"{istc_name}.istc")

    if method == "MFC":
        planner, paths, weights, runtime = mfc_plan(istc)
    elif method == "MSTC_Star":
        planner, paths, weights, runtime = mstcstar_plan(istc)
    elif method == "MIP":
        if istc_sol_name:
            sol_edges = Instance.read_solution(istc_sol_name)
        else:
            for fn in os.listdir(Instance.DEFAULT_SOLUTION_DIR):
                if fn.startswith(istc_name):
                    sol_edges = Instance.read_solution(fn[:-5])
                    break
        planner, paths, weights = mip_plan(istc, sol_edges)
    else:
        print(f"unsupported method {method}")
        return

    print(f"Method: {method}, Makespan: {max(weights)}")

    simulate(
        name=f"{istc_name}-{method if method != 'MSTC_Star' else 'MSTC*'}",
        planner=planner,
        paths=paths,
        weights=weights,
        scale=scale,
        dt=dt,
        obs_graph=nx.Graph(),
        is_write=is_write,
        is_show=not is_write
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("istc", help="Instance Name")
    parser.add_argument("--method", default="MIP",
                        help="Planner Type from {MFC, MSTC*, MIP}")
    parser.add_argument("--istc_sol_name", default=None,
                        help="MIP Solution Name Stored in 'data/solutions'")
    parser.add_argument("--scale", default=1.0, help="Plot Scaling Factor")
    parser.add_argument("--dt", default=0.02, help="Delta Time of Simulation")
    parser.add_argument("--write", default=False,
                        help="Is Writing the Simulation as MP4")

    args = parser.parse_args()

    plan(
        istc_name=args.istc,
        method=args.method,
        istc_sol_name=args.istc_sol_name,
        scale=float(args.scale),
        dt=float(args.dt),
        is_write=bool(args.write)
    )


if __name__ == "__main__":
    main()
