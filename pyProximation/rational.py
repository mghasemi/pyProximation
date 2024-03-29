from .base import Foundation
from .measure import Measure
from .orthsys import OrthSystem

from sympy import *


class RationalAprox(Foundation):
    """
    ``RationalAprox`` calculates a rational approximation for a given 
    function. It takes one argument `orth` which is an instance of
    ``OrthSystem`` and does all the computations in the scope of this
    object.
    """

    def __init__(self, orth):
        """
        initiate a rational approximation framework.
        """
        self.OrthSys = orth

    def RatLSQ(self, m, n, f):
        """
        Calculates a rational function :math:`\frac{p}{q}` where 
        :math:`p, q` consists of the first :math:`m, n` elements of 
        the orthonormal basis :math:`\|f-\frac{p}{q}\|_2` is minimized.
        """
        from numpy import diag, array, zeros, vstack, hstack
        from scipy.linalg import solve
        assert max(
            m, n) <= self.OrthSys.num_base, "The size of either numerator or denominator exceeds the size of basis."
        self.M = diag([1 for _ in range(m + 1)])
        self.Z = zeros((m + 1, n))
        self.S = zeros((n, n))
        self.F = array([self.OrthSys.inner(self.OrthSys.OrthBase[j], f)
                        for j in range(m + 1)])
        self.G = array(
            [-self.OrthSys.inner(self.OrthSys.OrthBase[k], f**2) for k in range(1, n + 1)])
        for j in range(m + 1):
            for k in range(1, n + 1):
                self.Z[j][
                    k - 1] = - self.OrthSys.inner(self.OrthSys.OrthBase[j], f * self.OrthSys.OrthBase[k])
        for j in range(1, n + 1):
            for k in range(1, n + 1):
                self.S[j - 1][k - 1] = self.OrthSys.inner(
                    self.OrthSys.OrthBase[j], f**2 * self.OrthSys.OrthBase[k])
        H = vstack((hstack((self.M, self.Z)), hstack((self.Z.T, self.S))))
        r = hstack((self.F, self.G))
        ab = solve(H, r)
        a, b = ab[:m + 1], ab[m + 1:]
        numer = sum([a[i] * self.OrthSys.OrthBase[i] for i in range(m + 1)])
        denom = 1 + sum([b[i] * self.OrthSys.OrthBase[i + 1]
                         for i in range(n)])
        return numer / denom

