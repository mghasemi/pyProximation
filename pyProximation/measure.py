from base import Foundation


class Measure(Foundation):
    """
    An instance of this class is a measure on a given set `supp`. The support is either
            + a python variable of type `set`, or
            + a list of tuples which represents a box in euclidean space.
    Initializes a measure object according to the inputs:
            + *dom* must be either
                    - a list of 2-tuples
                    - a non-empty dictionary
            + *w* must be a
                    - a function if `dom` defines a region
                    - left blank (None) if `dom` is a dictionary
    """

    def __init__(self, dom, w=None):
        """
        Initializes a measure object according to the inputs:
                +`dom` must be either
                        - a list of 2-tuples
                        - a non-empty dictionary
                + `w` must be a
                        - a function if `dom` defines a region
                        - left blank (None) if `dom` is a dictionary
        """
        self.ErrorMsg = ""
        self.dim = 0
        self.card = 0
        if not self.check(dom, w):
            raise Exception(self.ErrorMsg)

    def boxCheck(self, B):
        """
        Checks the structure of the box *B*.
        Returns `True` id `B` is a list of 2-tuples, otherwise it
        returns `False`.
        """
        flag = True
        for interval in B:
            flag = flag and (type(interval) == tuple) and (len(interval) == 2)
        return flag

    def check(self, dom, w):
        """
        Checks the input types and their consistency, according to the
        *__init__* arguments.
        """
        from types import FunctionType, IntType, LongType, FloatType
        if type(dom) == list:
            if not self.boxCheck(dom):
                self.ErrorMsg = "Each member of the support's list must be a tuple of 2 elements"
                return False
            self.DomType = "box"
            self.dim = len(dom)
            self.supp = dom
            if type(w) not in [FunctionType, IntType, LongType, FloatType]:
                self.ErrorMsg = "Weight must be a `function` defined over the support or a number"
                return False
            elif type(w) == FunctionType:
                self.weight = w
            else:
                self.weight = lambda *x: w
        elif type(dom) == dict:
            if len(dom) == 0:
                self.ErrorMsg = "A measure can not have an empty support"
                return False
            self.DomType = "set"
            self.card = len(dom)
            self.supp = dom.keys()
            self.weight = dom
        else:
            self.ErrorMsg = "First parameter must be either a list of 2-tupels or a dictionary"
            return False
        return True

    def measure(self, S):
        """
        Returns the measure of the set `S`.
        `S` must be a list of 2-tuples.
        """
        m = 0
        if self.DomType == "set":
            if type(S) not in [set, list, tuple]:
                raise Exception(
                    "The type of the given set must be either `set`, `list` or `tuple`")
            for p in S:
                if p in self.supp:
                    m += self.weight[p]
        else:
            if (type(S) != list) or (not self.boxCheck(S)):
                raise Exception("The given set must be a list of 2-tuples")
            from scipy import integrate
            m = integrate.nquad(self.weight, S)[0]
        return m

    def integral(self, f):
        """
        Returns the integral of `f` with respect to the currwnt measure
        over the support.
        """
        from types import FunctionType
        m = 0
        if self.DomType == "set":
            if type(f) not in [dict, FunctionType]:
                raise Exception(
                    "The integrand must be a `function` or a `dict`")
            if type(f) == dict:
                for p in self.supp:
                    if p in f:
                        m += self.weight[p] * f[p]
            else:
                for p in self.supp:
                    try:
                        m += self.weight[p] * f(p)
                    except:
                        pass
        else:
            if type(f) != FunctionType:
                raise Exception("The integrand must be a `function`")
            from scipy import integrate
            fw = lambda *x: f(*x) * self.weight(*x)
            m = integrate.nquad(fw, self.supp)[0]
        return m

    def norm(self, p, f):
        """
        Computes the norm-`p` of the `f` with respect to the current measure.
        """
        from math import pow
        absfp = lambda *x: pow(abs(f(*x)), p)
        return pow(self.integral(absfp), 1. / p)

    def sample(self, num):
        """
        Samples from the support according to the measure.

        """
        assert num >= 1, "Sample size must be a positive integer number."
        if self.DomType == 'box':
            from math import ceil, pow
            from itertools import product
            import random
            from random import uniform
            SubRegions = {}
            NumSample = {}
            points = []
            n = int(ceil(pow(num, (1. / self.dim))))
            delta = [(r[1] - r[0]) / float(n) for r in self.supp]
            for o in product(range(n), repeat=self.dim):
                SubRegions[o] = [(self.supp[i][0] + o[i] * delta[i], self.supp[i]
                                  [0] + (o[i] + 1) * delta[i]) for i in range(self.dim)]
            numpnts = max(num, len(SubRegions))
            muSupp = self.measure(self.supp)
            for o in SubRegions:
                NumSample[o] = ceil(numpnts * self.measure(SubRegions[o]))
            for o in NumSample:
                pnts = []
                while len(pnts) < NumSample[o]:
                    v = []
                    for rng in SubRegions[o]:
                        v.append(uniform(rng[0], rng[1]))
                    v = tuple(v)
                    if v not in pnts:
                        pnts.append(v)
                points += pnts
            return random.sample(set(points), num)
        else:
            from scipy import stats
            TotM = self.measure(self.supp)
            dist = [float(self.weight[p]) / TotM for p in self.supp]
            custm = stats.rv_discrete(name='custm', values=(self.supp, dist))
            return custm.rvs(size=num)
