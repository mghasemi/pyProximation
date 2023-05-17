from .base import Foundation
from .measure import Measure
from .orthsys import OrthSystem
from .collocation import Collocation


class SubRegion(Foundation):
    """
    The `SubRegion` class partitions the region into sub-regions, solve 
    the system of Integro-differential equations on each and glue them 
    together.

    It takes: 
            1) a collocation instance `collsys`;
            2) an optional list of positive integers `num_parts` which shows the number of equal length partitions for each variable.
    """

    def __init__(self, collsys, num_parts=[]):
        from itertools import product
        self.CollSys = collsys
        self.Vars = self.CollSys.Vars
        # find the domain
        self.CollSys.FindDomain()
        self.Domain = self.CollSys.Domain
        self.MiniCollSys = {}
        self.CornerCollPoints = True
        self.Solutions = {}
        if num_parts == []:
            self.Parts = [1 for _ in self.Vars]
        elif len(num_parts) != len(self.Vars):
            raise Exception(
                "Number of variables and dimension of partition does not match.")
        else:
            self.num_vars = len(self.Vars)
            self.Parts = num_parts
            self.breaks = [range(n) for n in self.Parts]
            self.tuples = product(*self.breaks)

    def InitMiniSys(self):
        """
        Initializes a set of collocation systems for all subregions.
        """
        from itertools import product
        self.tuples = product(*self.breaks)
        for tpl in self.tuples:
            # sub-region's collocation instance
            t_Coll = Collocation(
                self.Vars, self.CollSys.uFuncs, self.CollSys.Env)
            # the sub-region
            t_dom = [(self.Domain[idx][0] + tpl[idx] * (self.Domain[idx][1] - self.Domain[idx][0]) / float(self.Parts[idx]),
                      self.Domain[idx][0] + (tpl[idx] + 1) * (self.Domain[idx][1] - self.Domain[idx][0]) / float(self.Parts[idx])) for idx in range(self.num_vars)]
            # configure the orthonormal system of sub-region's collocation
            # instance
            for idx in range(len(self.CollSys.uFuncs)):
                cp_OrthSys = self.CollSys.OrthSys[idx]
                cp_t_dom = [t_dom[self.Vars.index(v)] for v in cp_OrthSys.Vars]
                t_OrthSys = OrthSystem(
                    cp_OrthSys.Vars, cp_t_dom, self.CollSys.Env)
                t_OrthSys.Basis(cp_OrthSys.OriginalBasis)
                t_OrthSys.FormBasis()
                # link
                t_Coll.SetOrthSys(t_OrthSys, self.CollSys.uFuncs[idx])
            # add equations
            t_Coll.Equation(self.CollSys.EQs)
            t_Coll.Verbose = self.CollSys.Verbose
            # append to the Sub-regions
            self.MiniCollSys[tpl] = t_Coll

    def AnalyseConditions(self):
        """
        Associates boundary conditions to each sub-collocation system.
        Starting from one corner, associates boundary conditions to all
        adjacent subregions based on the solution of the current subregion.
        """
        from itertools import product
        # walk through boundary conditions
        for num in range(len(self.CollSys.CndVals)):
            vals = self.CollSys.CndVals[num]
            tmp_tpl = [None for _ in range(self.num_vars)]
            for idx in range(len(vals)):
                if vals[idx] != self.Vars[idx]:
                    if vals[idx] == self.CollSys.Domain[idx][0]:  # an initial condition
                        tmp_tpl[idx] = 0
                    elif vals[idx] == self.CollSys.Domain[idx][1]:  # a final condition
                        tmp_tpl[idx] = self.Parts[idx] - 1
            corr_tpls = product(*self.breaks)
            # find corresponding tuples
            for idx in range(len(tmp_tpl)):
                corr_tpls = filter(lambda x: (x[idx] == tmp_tpl[
                                   idx] or tmp_tpl[idx] == None), corr_tpls)
            # append the condition to the corresponding regions
            for tpl in corr_tpls:
                self.MiniCollSys[tpl].Condition(self.CollSys.Cnds[num], vals)

    def corners(self, tpl):
        """
        Generates end corner points of the region represented `tpl` as 
        collocation points (used internally).
        """
        from itertools import product
        pnts = []
        corr_tpls = product(*self.breaks)
        corr_tpls = filter(lambda x: all(
            [(x[i] == tpl[i] or x[i] == tpl[i] + 1) for i in range(len(x))]), corr_tpls)
        corr_tpls = filter(lambda x: not all(
            [(x[i] == tpl[i] + 1) for i in range(len(x))]), corr_tpls)
        for tp in corr_tpls:
            pnt = []
            for i in range(self.num_vars):
                cord = self.Domain[i][
                    0] + tp[i] * (self.Domain[i][1] - self.Domain[i][0]) / float(self.Parts[i])
                pnt.append(cord)
            pnts.append(tuple(pnt))
        return pnts

    def Solve(self):
        """
        Solves the collocation system for each region and generates boundary conditions 
        for adjacent regions.
        """
        from itertools import product
        from sympy import Eq
        self.InitMiniSys()
        self.AnalyseConditions()
        corr_tpls = product(*self.breaks)
        for tpl in corr_tpls:
            if self.CollSys.Verbose:
                print("Region index:", tpl)
                print("-------------------------------")
            if self.CornerCollPoints:
                self.MiniCollSys[tpl].CollPoints(self.corners(tpl))
            Res = self.MiniCollSys[tpl].Solve()
            self.Solutions[tpl] = Res
            for idx in range(len(tpl)):
                next_tpl = list(tpl)
                if tpl[idx] < self.Parts[idx] - 1:
                    next_tpl[idx] += 1
                    cnd_val = [_ for _ in self.Vars]
                    cnd_val[idx] = self.Domain[idx][0] + next_tpl[idx] * \
                        (self.Domain[idx][1] - self.Domain[idx]
                         [0]) / float(self.Parts[idx])
                    for f in self.CollSys.uFuncs:
                        if f in Res:
                            self.MiniCollSys[tuple(next_tpl)].Condition(
                                Eq(f, Res[f]), cnd_val)

    def ClosedForm(self, func):
        """
        Compiles a piece-wise defined symbolic function from the solution 
        for `func` which is an unknown from the original collocation system.
        """
        from sympy import sign
        c_func = 0
        f_idx = self.CollSys.uFuncs.index(func)
        for idx in self.Solutions:
            mmbfunc = 1
            for v in self.CollSys.OrthSys[f_idx].Vars:
                i = self.CollSys.Vars.index(v)
                er_tol = 0.
                r0 = self.Domain[i][
                    0] + idx[i] * (self.Domain[i][1] - er_tol - self.Domain[i][0]) / float(self.Parts[i])
                r1 = self.Domain[i][
                    0] + (idx[i] + 1) * (self.Domain[i][1] - er_tol - self.Domain[i][0]) / float(self.Parts[i])
                mmbfunc *= (sign(self.Vars[i] - r0) + 1.) * \
                    (sign(r1 - self.Vars[i]) + 1.) / 4.
            if func in self.Solutions[idx]:
                c_func += mmbfunc * self.Solutions[idx][func]
        return c_func
