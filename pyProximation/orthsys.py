from base import Foundation
from measure import Measure


class OrthSystem(Foundation):
    """
    ``OrthogonalSystem`` class produces an orthogonal system of functions
    according to a suggested basis of functions and a given measure
    supported on a given region.

    This basically performs a 'Gram-Schmidt' method to extract the
    orthogonal basis. The inner product is obtained by integration of
    the product of functions with respect to the given measure (more
    accurately, the distribution).

    To initiate an instance of this class one should provide a list of
    symbolic variables `variables` and the range of each variable as a
    list of lists ``var_range``.

            To initiate an orthogonal system of functions, one should provide
            a list of symbolic variables ``variables`` and the range of each
            these variables as a list of lists ``var_range``.
    """

    def __init__(self, variables, var_range, env='sympy'):
        """
        To initiate an orthogonal system of functions, one should provide
        a list of symbolic variables ``variables`` and the range of each
        these variables as a list of lists ``var_range``.
        """
        assert (type(variables) is list) and (type(var_range) is list), """The OrthSystem class object
		requires two lists as inputs: (1) list of symbolic variables; (2) range of each variable."""
        self.EnvAvail = self.DetSymEnv()
        if self.EnvAvail == []:
            raise Exception("No Symbolic tool is available.")
        elif (env in self.EnvAvail):
            self.Env = env
        else:
            raise Exception("The selected symbolic tool is not supported.")
        self.Vars = variables
        self.num_vars = len(self.Vars)
        self.Domain = var_range
        self.measure = Measure(self.Domain, 1)
        self.OriginalBasis = []
        self.OrthBase = []
        self.Numerical = False
        self.CommonSymFuncs(self.Env)

    def PolyBasis(self, n):
        """
        Generates a polynomial basis from variables consisting of all
        monomials of degree at most ``n``.
        """
        assert n >= 0, "'n' must be a positive integer."
        from itertools import product
        B = []
        for o in product(range(n + 1), repeat=self.num_vars):
            if sum(o) <= n:
                T_ = 1
                for idx in range(self.num_vars):
                    T_ *= self.Vars[idx]**o[idx]
                B.append(T_)
        return B

    def FourierBasis(self, n):
        """
        Generates a Fourier basis from variables consisting of all
        :math:`sin` & :math:`cos` functions with coefficients at most `n`.
        """
        assert n >= 0, "'n' must be a positive integer."
        from itertools import product
        B = []
        for o in product(range(n + 1), repeat=self.num_vars):
            if sum(o) <= n:
                SinCos = product(range(2), repeat=self.num_vars)
                for ex in SinCos:
                    T_ = 1
                    for idx in range(self.num_vars):
                        period = self.Domain[idx][1] - self.Domain[idx][0]
                        if o[idx] != 0:
                            if ex[idx] == 0:
                                T_ *= self.cos(2 * self.pi *
                                               o[idx] * self.Vars[idx] / period)
                            else:
                                T_ *= self.sin(2 * self.pi *
                                               o[idx] * self.Vars[idx] / period)
                    B.append(T_)
        return list(set(B))

    def TensorPrd(self, Bs):
        """
        Takses a list of symbolic bases, each one a list of symbolic 
        expressions and returns the tensor product of them as a list.
        """
        assert (Bs != []), "Can not compute the tensor product of empty bases."
        from itertools import product
        TP = product(*Bs)
        TBase = []
        for itm in TP:
            t_prd = 1
            for ent in itm:
                t_prd = t_prd * ent
            TBase.append(self.expand(t_prd))
        return TBase

    def SetMeasure(self, M):
        """
        To set the measure which the orthogonal system will be computed,
        simply call this method with the corresponding distribution as
        its parameter `dm`; i.e, the parameter is `d(m)` where `m` is
        the original measure.
        """
        assert isinstance(M, Measure), "The argument must be a `Measure`."
        self.measure = M

    def Basis(self, base_set):
        """
        To specify a particular family of function as a basis, one should
        call this method with a list ``base_set`` of linearly independent
        functions.
        """
        assert type(
            base_set) is list, "A list of symbolic functions is expected."
        self.OriginalBasis = base_set
        self.num_base = len(self.OriginalBasis)

    def inner(self, f, g):
        """
        Computes the inner product of the two parameters with respect to
        the measure ``measure``.
        """
        if self.Env == "sympy":
            from sympy import lambdify
            F = lambdify(self.Vars, f * g, "numpy")
        elif self.Env == "sage":
            from sage.all import fast_callable
            h = f * g + self.Vars[0] * 0
            H = fast_callable(h, vars=self.Vars)
            F = lambda *x: H(*x)
        elif self.Env == 'symengine':
            from symengine import Lambdify
            F = lambda *x: Lambdify(self.Vars, [f * g])(x)[0]
        m = self.measure.integral(F)
        return m

    def project(self, f, g):
        """
        Finds the projection of ``f`` on ``g`` with respect to the inner
        product induced by the measure ``measure``.
        """
        return g * self.inner(f, g) / self.inner(g, g)

    def FormBasis(self):
        """
        Call this method to generate the orthogonal basis corresponding
        to the given basis via ``Basis`` method.
        The result will be stored in a property called ``OrthBase`` which
        is a list of function that are orthogonal to each other with
        respect to the measure ``measure`` over the given range ``Domain``.
        """
        for f in self.OriginalBasis:
            nf = 0
            for u in self.OrthBase:
                nf += self.project(f, u)
            nf = f - nf
            F = self.expand(nf / self.sqrt(self.inner(nf, nf)))
            self.OrthBase.append(F)
        self.num_base = len(self.OrthBase)

    def SetOrthBase(self, base):
        """
        Sets the orthonormal basis to be the given `base`.
        """
        assert (base != []), "Invalid basis."
        self.OrthBase = base
        self.num_base = len(self.OrthBase)

    def Series(self, f):
        """
        Given a function `f`, this method finds and returns the
        coefficients of the	series that approximates `f` as a
        linear combination of the elements of the orthogonal basis.
        """
        cfs = []
        for b in self.OrthBase:
            cfs.append(self.inner(f, b))
        return cfs
