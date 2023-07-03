import argparse
import os

import yaml

from MIP_MCPP.instance import Instance
from MIP_MCPP.model import Model
from MIP_MCPP.warmstarter import WarmStarter


def solve(
    istc_name: str,
    solver_args: dict,
    alpha: float,
    beta: float,
    warm_start: str
) -> None:

    istc = Instance.read(istc_name + ".istc")

    model = Model(istc)
    H, perc_vars_removed = model.apply_heur(alpha, beta)
    model.wrapup(solver_args, H)
    if warm_start:
        model = WarmStarter.apply(model, warm_start, H)
    sol_edges, sol_verts = model.solve()

    sol_name = "" if alpha is None and beta is None else f"-alpha_{alpha}-beta_{beta}_warmstart{warm_start}"
    res_str = ",\t".join(str(s) for s in [
        istc_name,
        "None" if alpha is None else alpha,
        "None" if beta is None else beta,
        model.num_vars if alpha is None and beta is None else round(
            perc_vars_removed, 3)
    ])

    if sol_edges == []:
        obj_val, mip_gap = "/", "/"
    else:
        obj_val = round(model.model.objVal, 3)
        mip_gap = round(model.model.MIPGap, 3),
        Instance.write_solution(sol_edges, istc_name + sol_name + ".solu")

    res_str += ",\t".join(str(s) for s in [
        "",
        warm_start,
        obj_val,
        mip_gap,
        round(model.model.objBound, 3),
        round(model.model.RunTime, 3)
    ])

    with open(os.path.join("data", "res.txt"), "a") as f:
        f.write(res_str + "\n")

    return sol_edges, sol_verts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("istc", help="Instance Name")
    parser.add_argument("--solver_cfg", default="default.yaml",
                        help="Gurobi Configuration File Path")
    parser.add_argument("--alpha", default=None,
                        help="Parabolic Heuristics Parameter")
    parser.add_argument("--beta", default=None,
                        help="Subgraph Heuristics Parameter")
    parser.add_argument("--warm_start", default=None,
                        help="Warm startup the optimization: [RTC] for non-heuritiscs model and [MST] for heuritsics model")

    args = parser.parse_args()

    with open(os.path.join("data", "cfgs", args.solver_cfg)) as f:
        solver_args = yaml.load(f, yaml.Loader)
        solver_args["OptimalityTol"] = float(solver_args["OptimalityTol"])

    solve(
        istc_name=args.istc,
        solver_args=solver_args,
        alpha=float(args.alpha) if args.alpha else None,
        beta=float(args.beta) if args.beta else None,
        warm_start=args.warm_start if args.warm_start else None
    )


if __name__ == "__main__":
    main()
