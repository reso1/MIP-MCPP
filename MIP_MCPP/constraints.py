from enum import Enum
from typing import Dict, List

import numpy as np
import scipy.sparse as sp


class BoundType(Enum):
    LowerBound = 0
    UpperBound = 1
    EqualBound = 2


class Constraints:

    def __init__(
        self,
        name: str,              # constraint group name
        num_constraints: int,   # number of constraints
        bound: float,           # bound value
        bound_type: BoundType,  # bound type
        num_e: List[int],       # list of number of edges foreach i
        num_v: List[int],       # list of number of verts foreach i
        has_z=False           # has z variables
    ) -> None:

        self.k = len(num_e)
        self.name = name
        self.num_constraints = num_constraints
        self.bound_type = bound_type
        self.has_z = has_z

        # [x^1_e1, x^1_e2, ...] + [x^2_e1, x^2_e2, ...] + ... + [x^k_e1, x^k_e2, ...]
        self.x_coeffs = np.zeros((num_constraints, sum(num_e)))
        # [y^1_v1, y^1_v2, ...] + [y^2_v1, y^2_v2, ...] + ... + [y^k_v1, y^k_v2, ...]
        self.y_coeffs = np.zeros((num_constraints, sum(num_v)))
        # [fu^1_e1, fu^1_e2, ...] + [fu^2_e1, fu^2_e2, ...] + ... + [fu^k_e1, fu^k_e2, ...]
        self.fu_coeffs = np.zeros((num_constraints, sum(num_e)))
        # [fv^1_e1, fv^1_e2, ...] + [fv^2_e1, fv^2_e2, ...] + ... + [fv^k_e1, fv^k_e2, ...]
        self.fv_coeffs = np.zeros((num_constraints, sum(num_e)))
        self.tau_coeff = 0

        self.bound = bound

        if self.has_z:
            # [z^1_v1, z^1_v2, ...] + [z^2_v1, z^2_v2, ...] + ... + [z^k_v1, z^k_v2, ...]
            self.z_coeffs = np.zeros((num_constraints, sum(num_v)))

    def build(self, x, y, fu, fv, tau, z=None):

        lhs = sp.csr_matrix(self.x_coeffs)  @ x + \
              sp.csr_matrix(self.y_coeffs)  @ y + \
              sp.csr_matrix(self.fu_coeffs) @ fu + \
              sp.csr_matrix(self.fv_coeffs) @ fv + \
              self.tau_coeff * tau

        if z is not None and self.has_z:
            lhs += sp.csr_matrix(self.z_coeffs) @ z

        bounds = self.bound * np.ones((self.num_constraints, 1))

        if self.bound_type == BoundType.LowerBound:
            return lhs >= bounds, self.name
        elif self.bound_type == BoundType.UpperBound:
            return lhs <= bounds, self.name
        else:
            return lhs == bounds, self.name

    def __str__(self, x_ind, y_ind, fu_ind, fv_ind) -> str:
        def __sub_str(coeff, var_name):
            if coeff > 0:
                return f" +{coeff:.0f}*{var_name}[{i},{j}]"
            elif coeff < 0:
                return f" {coeff:.0f}*{var_name}[{i},{j}]"
            else:
                return ""

        s = f"{self.name}:\n"
        for ci in range(self.num_constraints):
            sub_s = f"({ci}): \t" + \
                f"{self.tau_coeff}*tau" if self.tau_coeff != 0 else ""

            for i in range(self.k):
                for j in range(x_ind[i], x_ind[i+1]):
                    sub_s += __sub_str(self.x_coeffs[ci, j], 'x')
                for j in range(y_ind[i], y_ind[i+1]):
                    sub_s += __sub_str(self.y_coeffs[ci, j], 'y')
                for j in range(fu_ind[i], fu_ind[i+1]):
                    sub_s += __sub_str(self.fu_coeffs[ci, j], 'fu')
                for j in range(fv_ind[i], fv_ind[i+1]):
                    sub_s += __sub_str(self.fv_coeffs[ci, j], 'fv')

            if self.bound_type == BoundType.LowerBound:
                sub_s += f" >= {self.bound}\n"
            elif self.bound_type == BoundType.UpperBound:
                sub_s += f" <= {self.bound}\n"
            else:
                sub_s += f" == {self.bound}\n"

            s += sub_s

        return s
