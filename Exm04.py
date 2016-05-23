from sympy import *
from Measure import *

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
print intrpl
G = plot(f, intrpl, (x, -3.1, 3.1), show=False)
G[0].line_color = 'blue'
G[1].line_color = 'green'
G.save('Exm04-%d.png'%(n))