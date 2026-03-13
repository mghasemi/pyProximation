# Interpolation and least square approximation w.r.t. a discrete measure
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
	sys.path.insert(0, str(ROOT))

from sympy import Symbol, Function, exp, sin, pprint
from pyProximation import Measure, OrthSystem

x = Symbol('x')
y = Function('y')(x)
M = Measure({-3:1, -2:1, -1:1, 0:1, 1:1, 2:1, 3:1})
S = OrthSystem([x], [(-3, 3)])
S.SetMeasure(M)

n = 6

B = S.PolyBasis(n)

m = len(B)

S.Basis(B)

S.FormBasis()

f = exp(x*sin(x))
sries = S.Series(f)
intrpl = sum([S.OrthBase[i]*sries[i] for i in range(m)])

pprint(intrpl)
