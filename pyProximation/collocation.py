from .base import Foundation
from .measure import Measure
from .orthsys import OrthSystem


class Collocation(Foundation):
    """
    The ``Collocation`` class tries to approximate the solutions of a system
    of partial differential equations with respect to an orthogonal
    system of functions.

    To initiate an instance of this class one needs to provide two set of parameters:
            1) List of independent symbolic variables `variables`;
            2) List of unknown functions to be found that depend on the independent variables ``ufunc``.
    """

    def __init__(self, variables, ufunc, env='sympy'):
        Env = self.DetSymEnv()
        if Env == []:
            raise Exception("No Symbolic tool is available.")
        if env not in Env:
            raise Exception("The selected symbolic tool is not available.")
        self.Env = env  # The selected tool for symbolic computations
        self.CommonSymFuncs(self.Env)
        self.Vars = variables  # Symbolic variables
        self.num_vars = len(variables)  # Number of symbolic variables
        self.uFuncs = ufunc  # Unknown functions
        self.num_funcs = len(ufunc)  # Number of unknown functions
        # Number of elements in the orthogonal basis
        self.degree = [1 for _ in ufunc]
        self.EQs = []  # Lists of functional equations
        self.Cnds = []  # Storage for initial and boundary conditions
        self.CndVals = []  # Storage for the values of CND
        self.Coeffs = {}
        self.Apprx = {}  # Approximate solutions to uFuncs
        self.Points = []  # Collocation points
        self.OrthSys = {}  # Orthogonal systems of functions corresponding to uFuncs
        self.Solver = 'scipy'  # The solver to find the roots
        self.SolverOption = 'lm'  # Specifies scipy solver
        # Reserved for the domain of variables
        self.Domain = [(None, None) for v in self.Vars]
        self.SampleMeasure = None  # Reserved for the sampling measure
        self.Verbose = False  # Set True to see some messages about the procedure
        # Determines the final status of the solver: `True` for success and
        # `False` for fail
        self.Success = None
        self.init_guess_bnd = 0.1
        # the initial point for solver. Its dimension must be equal to number
        # of unknown coefficients
        self.InitPoint = []
        self.CfSyms = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
                       'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w'
                       'x', 'y', 'z']

    def SetOrthSys(self, obj, func):
        """
        To approximate the solutions of the system of pdes, the class
        requires an orthogonal system of functions ``OrthSystem``.
        This method accepts such a system.
        """
        assert isinstance(
            obj, OrthSystem), "An object of type `OrthSystem` is expected."
        assert func in self.uFuncs, "`func` must be a function symbol given at initiation."
        idx = self.uFuncs.index(func)
        self.OrthSys[idx] = obj
        self.degree[idx] = obj.num_base

    def FindDomain(self):
        """
        Finds the region that all variables are defined.
        """
        dom = [list(tpl) for tpl in self.Domain]
        for v_idx in range(self.num_vars):
            v = self.Vars[v_idx]
            for f in self.OrthSys:
                Orth = self.OrthSys[f]
                if v in Orth.Vars:
                    idx = Orth.Vars.index(v)
                    rng = Orth.Domain[idx]
                    if (dom[v_idx][0] == None) or (rng[0] < dom[v_idx][0]):
                        dom[v_idx][0] = rng[0]
                    if (dom[v_idx][1] == None) or (rng[1] > dom[v_idx][1]):
                        dom[v_idx][1] = rng[1]
        self.Domain = [tuple(lst) for lst in dom]
        # defines the default sampling measure object
        self.SampleMeasure = Measure(self.Domain, 1)

    def setSampleMeasure(self, meas):
        """
        Sets the measure over the domain for sampling collocation points in case of necessity.
        """
        assert isinstance(
            obj, Measure), "The input must be an instance of a `Measure` object."
        self.SampleMeasure = meas

    def Equation(self, eq):
        """
        To enter the system of equations, use this meyhod with a list
        of equations as input.
        """
        if type(eq) is list:
            self.EQs += eq
        else:
            self.EQs.append(eq)

    def Condition(self, eq, val):
        """
        List of initial and boundary conditions.
        """
        # if eq not in self.Cnds:
        self.Cnds.append(eq)
        self.CndVals.append(val)

    def setSolver(self, solver, optn='lm'):
        """
        Currently only two solvers are supported:

                1. the `sage`\`s defult solver for rather simple system of algebraic equations.
                2. the `scipy`\`s `fsolves` to handel more complex and larger systems.
                It also supports the following solvers from scipy:

                        + `hybr`
                        + `lm` (defauls)
                        + `broyden1`
                        + `broyden2`
                        + `anderson`
                        + `krylov`
                        + `df-sane`
        """
        self.Solver = solver.lower()
        self.SolverOption = optn

    def CollPoints(self, pnts):
        """
        Accepts alist of collocation point ``pnts``, to form the algebraic
        system of equations and find the coefficients of the orthogonal
        functions from ``OrthSystem.OrthBase``. Each point must be
        either a list or a tuple.
        """
        self.Points += pnts

    def collocate(self):
        """
        Internal use: generates the system of equations for coefficients
        to be used via collocation points.
        """
        """if self.Env == 'sympy':
			from sympy import Symbol as var
			from sympy import Subs, expand, diff
		elif self.Env == 'sage':
			from sage.all import var, expand, diff
		elif self.Env == 'symengine':
			from symengine import Symbol as var
			from symengine import expand, diff"""
        # symbols for coefficients
        var_syms = self.CfSyms
        # Produce enough symbolic variables
        self.CF = {var_syms[s]: [self.Symbol('%s%d' % (var_syms[s], i)) for i in range(
            self.degree[s])] for s in range(self.num_funcs)}
        self.SymCF = []
        for s in self.CF:
            self.SymCF += self.CF[s]
        self.SR = {}
        self.REq = []
        # loop over unknown functions to be found
        for f_idx in range(self.num_funcs):
            s = var_syms[f_idx]
            T = 0
            # loop over elements of orthogonal basis
            for i in range(self.degree[f_idx]):
                T += self.CF[s][i] * (self.OrthSys[f_idx].OrthBase[i])
            self.SR[s] = T
        # loop over entered equations
        EQ_num = 0
        for eq in self.EQs:
            f_idx = 0
            Teq = eq
            for f in self.uFuncs:
                for v in self.Vars:
                    if self.Env == 'sage':
                        for d_ord in range(1, 6):
                            Teq = Teq.subs({self.diff(f, v, d_ord): self.diff(
                                self.SR[var_syms[f_idx]], v, d_ord)})
                    if self.Env == 'sympy':
                        Teq = (Teq.subs({f: self.SR[var_syms[f_idx]]})).doit()
                    elif self.Env == 'symengine':
                        Teq = (Teq.msubs({f: self.SR[var_syms[f_idx]]}))
                        Teq = (
                            Teq.msubs({self.diff(f, v): self.diff(self.SR[var_syms[f_idx]], v)}))
                        Teq = (Teq.msubs({self.diff(self.diff(f, v), v): self.diff(
                            self.diff(self.SR[var_syms[f_idx]], v), v)}))
                f_idx += 1
            if Teq not in self.REq:
                self.REq.append(Teq)
                EQ_num += 1
                if self.Verbose:
                    print("Equation # %d generated." % (EQ_num))
        # loop over initial and boundary conditions
        cnd_idx = 0
        for eq in self.Cnds:
            f_idx = 0
            Teq = eq
            for f in self.uFuncs:
                for v_idx in range(self.num_vars):
                    v = self.Vars[v_idx]
                    Teq = Teq.subs(
                        {self.diff(f, v): self.diff(self.SR[var_syms[f_idx]], v)})
                    if self.Env == 'sympy':
                        Teq = (Teq.subs({f: self.SR[var_syms[f_idx]]})).doit()
                    elif self.Env == 'symengine':
                        Teq = self.expand(
                            Teq.msubs({f: self.SR[var_syms[f_idx]]}))
                    else:
                        Teq = self.expand(
                            Teq.subs({f: self.SR[var_syms[f_idx]]}))
                f_idx += 1
            # Teq = self.expand(Teq.subs({self.Vars[v]:self.CndVals[cnd_idx][v]
            # for v in range(len(self.CndVals[cnd_idx]))})) # needs more work
            # (index of variables could be off)
            Teq = Teq.subs({self.Vars[v]: self.CndVals[cnd_idx][v]
                            for v in range(len(self.CndVals[cnd_idx]))})
            if Teq not in self.REq:
                self.REq.append(Teq)
            cnd_idx += 1
            if self.Verbose:
                print("Condition # %d added." % (cnd_idx))

    def PlugPoints(self):
        """
        Internal use: plug in collocation points to elliminate independent
        variables and keep the coefficients.
        """
        # Plug in the collocation points to form the algebraic equations

        numeric_eqs = []
        for p in self.Points:
            chg = {self.Vars[i]: p[i] for i in range(self.num_vars)}
            for eq in self.REq:
                tp1 = type(eq)
                Teq = eq.subs(chg)
                tp2 = type(Teq)
                if (tp1 == tp2) and (Teq not in numeric_eqs):
                    numeric_eqs.append(Teq)
                if len(numeric_eqs) >= len(self.SymCF):
                    break
            if len(numeric_eqs) >= len(self.SymCF):
                break

        if len(numeric_eqs) != len(self.SymCF):
            raise Exception(
                "Number of points and equations are not equal! Check the conditions.")
        if self.Verbose:
            print("Solving the system of equations numerically to extract coefficients ...")
        # Solve the algebraic equations
        if self.Solver == 'sage':
            if self.Env != 'sage':
                raise Exception(
                    "Sage solver is not available in selected symbolic environment.")
            sols = solve(numeric_eqs, self.SymCF, solution_dict=True)
            sols = sols[0]
            return sols
        elif self.Solver in ['scipy']:
            from scipy import optimize as opt
            from random import uniform
            if self.Env == 'sympy':
                from sympy import lambdify
                f_ = [lambdify(self.SymCF, (eq.lhs - eq.rhs), "numpy")
                      for eq in numeric_eqs]

                def f(x):
                    z = tuple(float(x.item(i)) for i in range(len(self.SymCF)))
                    return [fn(*z) for fn in f_]
            elif self.Env == 'symengine':
                from symengine import sympify
                from sympy import lambdify
                t_eqs = [sympify(eq) for eq in numeric_eqs]
                f_ = [lambdify(self.SymCF, eq, "numpy") for eq in t_eqs]

                def f(x):
                    z = tuple(float(x.item(i)) for i in range(len(self.SymCF)))
                    return [fn(*z) for fn in f_]
            elif self.Env == 'sage':
                def f(x):
                    chng = {}
                    U = self.SymCF
                    n_var = len(U)
                    chng = {U[i]: float(x.item(i)) for i in range(n_var)}
                    EQs_ = []
                    for eq in numeric_eqs:
                        teq = eq.lhs() - eq.rhs()
                        EQs_.append(teq.subs(chng).n())
                    return EQs_
            nvars = len(self.SymCF)
            if self.Solver == 'scipy':
                if self.InitPoint != []:
                    init_point = tuple(self.InitPoint)
                else:
                    init_point = tuple(
                        uniform(-self.init_guess_bnd, self.init_guess_bnd) for _ in range(nvars))
                sol = opt.root(f, init_point, method=self.SolverOption)
            if sol.success:
                sols = {self.SymCF[i]: list(sol.x)[i] for i in range(nvars)}
                self.Success = True
            else:
                sols = {}
                self.Success = False
            if self.Verbose:
                print(sol.message)
            return sols

    def Solve(self):
        """
        Solves the collocation equations and keep a dictionary of
        coefficients in ``self.Coeffs`` and returns a list of functions
        in the span of orthoginal system.
        """
        if self.Verbose:
            print("Check for collocation points shortage...")
        if self.SampleMeasure is None:
            self.FindDomain()
        num = max(self.degree) - len(self.Points)
        # Check for too many points
        if num < 0:
            raise Exception(
                "Too many points are associated. Reduce at least %d" % (-num))
        cl_points = []
        # Add enough random points to match up for variables
        if num > 0:
            if self.Verbose:
                print("Generating %d new collocation point ..." % (num))
            cl_points = self.SampleMeasure.sample(num)
        # attaching points
        self.CollPoints(cl_points)
        if self.Verbose:
            print("Attaching %d collocation points:" % (len(self.Points)))
            for p in self.Points:
                print(p)
        if self.Verbose:
            print("Generating algebraic equations based on given orthogonal systems of functions ...")
        self.collocate()
        if self.Verbose:
            print("Plug collocation points to extract system of equations ...")
        self.Coeffs = self.PlugPoints()
        if self.Verbose:
            print("Done!")
        if self.Coeffs != {}:
            for fn in self.uFuncs:
                s = self.CfSyms[self.uFuncs.index(fn)]
                self.Apprx[fn] = self.SR[s].subs(self.Coeffs)
        return self.Apprx
