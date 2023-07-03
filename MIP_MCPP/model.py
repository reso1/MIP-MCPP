from collections import defaultdict
from typing import List, Tuple

import gurobipy as gp
import networkx as nx
import numpy as np
from gurobipy import GRB

from MIP_MCPP.constraints import BoundType as BT
from MIP_MCPP.constraints import Constraints
from MIP_MCPP.heuristics import heur_parabolic, heur_subgraph
from MIP_MCPP.instance import Instance
from MIP_MCPP.misc import SEP_LITERAL, uv_sorted


class Model:
    """ RMMTC Model """

    def __init__(self, istc:Instance) -> None:
        self.istc = istc

    def _init_var_info(self, H:List[nx.Graph]) -> None:
        self.num_vars = 1
        self.num_v, self.num_e = [], []
        
        for i in self.istc.I:
            self.num_e.append(self.istc.m if H is None else H[i].number_of_edges())
            self.num_v.append(self.istc.n if H is None else H[i].number_of_nodes())
            self.num_vars += self.num_v[-1] + 3 * self.num_e[-1]
        
        self.e_ind = [0] + np.cumsum(self.num_e).tolist()
        self.v_ind = [0] + np.cumsum(self.num_v).tolist()

    def _init_constrs(self, H:List[nx.Graph]) -> None:
        if H is None:
            H = [self.istc.G for i in self.istc.I]

        self.v2id, self.id2v, self.e2id, self.id2e = [], [], [], []
        for i in self.istc.I:
            self.v2id.append({v:idx for idx, v in enumerate(H[i].nodes())})
            self.id2v.append({idx:v for idx, v in enumerate(H[i].nodes())})
            self.e2id.append({self.istc.uv2e[uv_sorted(uv)]:idx for idx, uv in enumerate(H[i].edges())})
            self.id2e.append({idx:self.istc.uv2e[uv_sorted(uv)] for idx, uv in enumerate(H[i].edges())})
        
        self.num_constrs = 3 * self.istc.k + self.istc.n + sum(self.num_e) + sum(self.num_v)
        self.C_makespan  = Constraints('makespan',  self.istc.k,  0, BT.UpperBound, self.num_e, self.num_v)
        self.C_rooted    = Constraints('rooted',    self.istc.k,  1, BT.EqualBound, self.num_e, self.num_v)
        self.C_tree      = Constraints('tree',      self.istc.k,  1, BT.EqualBound, self.num_e, self.num_v)
        self.C_cover     = Constraints('cover',     self.istc.n,  1, BT.LowerBound, self.num_e, self.num_v)
        self.flow        = Constraints('flow',      sum(self.num_e), 0, BT.EqualBound, self.num_e, self.num_v)
        self.C_acyclic   = Constraints('acyclic',   sum(self.num_v), 1-1/self.istc.n, BT.UpperBound, self.num_e, self.num_v)
        self.C_y_defs = []

        # -----------------------------------------------------------------------------
        for i in self.istc.I:
            # CONSTRAINTS(makespan):        sum_e{w_e * x_e^i} - tau <= 0
            self.C_makespan.tau_coeff = -1
            self.C_makespan.x_coeffs[i, self.e_ind[i]:self.e_ind[i+1]] = \
                list(nx.get_edge_attributes(H[i], "weight").values())

            # CONSTRAINTS(rooted):                     y_ri^i = 1
            self.C_rooted.y_coeffs[i, self.v_ind[i] + self.v2id[i][self.istc.R[i]]] = 1

            # CONSTRAINTS(tree): -sum_e{x_e^i} + sum_v{y_v^i} = 1
            self.C_tree.x_coeffs[i, self.e_ind[i]:self.e_ind[i+1]] = -1
            self.C_tree.y_coeffs[i, self.v_ind[i]:self.v_ind[i+1]] = 1

        # -----------------------------------------------------------------------------
        for v in self.istc.V:
            for i in self.istc.I:
                if v in self.v2id[i]:
                    # CONSTRAINTS(cover): sum_i{y_v^i} >= 1
                    self.C_cover.y_coeffs[v, self.v_ind[i]+self.v2id[i][v]] = 1

        # -----------------------------------------------------------------------------
        for i in self.istc.I:
            for e in self.istc.E:
                if e in self.e2id[i]:
                    e_id = self.e_ind[i] + self.e2id[i][e]
                    # CONSTRAINTS(flow): - x_e^i + fu_e^i + fv_e^i = 0
                    self.flow.x_coeffs[e_id, e_id] = -1
                    self.flow.fu_coeffs[e_id, e_id] = 1
                    self.flow.fv_coeffs[e_id, e_id] = 1

        # -----------------------------------------------------------------------------
        v_ngbs_edge = [defaultdict(list) for _ in self.istc.I]
        for i in self.istc.I:
            for v in self.istc.V:
                if H[i].has_node(v):
                    for u in H[i].neighbors(v):
                        uv = uv_sorted((u, v))
                        v_ngbs_edge[i][v].append((uv, self.istc.uv2e[uv]))

        for v in self.istc.V:
            for i in self.istc.I:
                if len(v_ngbs_edge[i][v]) == 0:
                    continue

                self.num_constrs += len(v_ngbs_edge[i][v])
                C_y_def = Constraints(f'y_def_({i},{v})', len(v_ngbs_edge[i][v]), 0, BT.UpperBound, self.num_e, self.num_v)
                for j, val in enumerate(v_ngbs_edge[i][v]):
                    uv, e = val
                    v_id = self.v_ind[i] + self.v2id[i][v]
                    e_id = self.e_ind[i] + self.e2id[i][e]
                    # CONSTRAINTS(acyclic): sum_e{fv_e^i} <= 1 - 1 / |V|
                    if v == uv[0]:
                        self.C_acyclic.fu_coeffs[v_id, e_id] = 1
                    else:
                        self.C_acyclic.fv_coeffs[v_id, e_id] = 1
                    
                    # CONSTRAINTS(y_def): x_e^i - y_v^i <= 0
                    C_y_def.x_coeffs[j, e_id] = 1
                    C_y_def.y_coeffs[j, v_id] = -1

                self.C_y_defs.append(C_y_def)

    def _init_model_params(self, args:dict):
        """ Use GUROBI to solve the RMMTC instance """
        model = gp.Model("RMMTC")
                
        model.params.Threads = int(args["Threads"])
        model.params.OptimalityTol = float(args["OptimalityTol"])
        if args.get("TimeLimit"):
            model.params.TimeLimit = float(args["TimeLimit"])
        if args.get("SoftMemLimit"):
            model.params.SoftMemLimit = float(args["SoftMemLimit"])
        
        print(
            "\n" + SEP_LITERAL + \
            f"Number of trees = {self.istc.k}\n" + \
            f"Number of verts = {self.istc.n}\n" + \
            f"Number of edges = {self.istc.m}\n" + \
            f"Number of variables = {self.num_vars}\n" + \
            f"Number of constraints = {self.num_constrs}\n" + \
            SEP_LITERAL
        )

        return model

    def wrapup(self, args:dict, H:List[nx.Graph]=None) -> None:
        self._init_var_info(H)
        self._init_constrs(H)
        model = self._init_model_params(args)

        """ Variables """
        tau = model.addVar(name="tau", vtype=GRB.CONTINUOUS)
        x   = model.addMVar(name="x",  shape=(sum(self.num_e), 1),  vtype=GRB.BINARY)
        y   = model.addMVar(name="y",  shape=(sum(self.num_v), 1),  vtype=GRB.BINARY)
        fu  = model.addMVar(name="fu", shape=(sum(self.num_e), 1), vtype=GRB.CONTINUOUS)
        fv  = model.addMVar(name="fv", shape=(sum(self.num_e), 1), vtype=GRB.CONTINUOUS)
        
        """ Objective """
        model.setObjective(tau, GRB.MINIMIZE)

        """ Constraints """
        constrs = [self.C_makespan, self.C_rooted, self.C_tree,self.C_cover,
                   self.flow, self.C_acyclic] + self.C_y_defs
        for C in constrs:
            constr, name = C.build(x, y, fu, fv, tau)
            model.addConstr(constr, name = name)
        
        self.model, self.x, self.y, self.fu, self.fv = model, x, y, fu, fv

    def solve(self) -> Tuple[list, list]:

        try:
            self.model.optimize()
        except gp.GurobiError as e:
            print(f"Error code {e.errno}: {e}")
        except AttributeError:
            print('Encountered an attribute error')
        
        sol_edges = [[] for _ in self.istc.I]
        sol_verts = [[] for _ in self.istc.I]

        try:
            for i in self.istc.I:
                for e_id in range(self.e_ind[i], self.e_ind[i+1]):
                    if self.x.X[e_id] > 0.5:
                        e = self.id2e[i][e_id - self.e_ind[i]]
                        sol_edges[i].append(self.istc.e2uv[e])
                for v_id in range(self.v_ind[i], self.v_ind[i+1]):
                    if self.y.X[v_id] > 0.5:
                        v = self.id2v[i][v_id-self.v_ind[i]]
                        sol_verts[i].append(v) 
        except Exception as e:
            print(e)
        
        return sol_edges, sol_verts

    def apply_heur(self, alpha:float, beta:float) -> Tuple[List[nx.Graph], int]:
        if alpha is not None and beta is None:
            Vi, _, _ = heur_parabolic(self.istc, alpha)
        elif alpha is None and beta is not None:
            Vi, _, _ = heur_subgraph(self.istc, beta)
        elif alpha is not None and beta is not None:
            V_prh, _, _ = heur_parabolic(self.istc, alpha)
            V_srh, _, _ = heur_subgraph(self.istc, beta)
            Vi = {r:set.union(V_prh[r], V_srh[r]) for r in self.istc.R}
        else:
            return [self.istc.G for _ in self.istc.I], 0

        H = []
        num_vars_removed = 0
        for i in self.istc.I:
            residual_graph_verts = set(self.istc.V) - Vi[self.istc.R[i]]
            Hi = self.istc.G.subgraph(residual_graph_verts)
            assert(nx.is_connected(Hi))
            H.append(Hi)
            num_vars_removed += H[i].number_of_nodes() + 3 * H[i].number_of_edges()

        return H, 1 - num_vars_removed / (self.istc.k * (self.istc.n + 3 * self.istc.m))
