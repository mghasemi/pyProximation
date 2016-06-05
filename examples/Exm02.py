from sympy import *
from pyProximation import *

x = Symbol('x')
y = Function('y')(x)

n = 10

S = OrthSystem([x], [(0, 2*pi)])
B = S.PolyBasis(n)

S.Basis(B)
S.FormBasis()

EQ1 = Eq(diff(y, x) + 2*y +5*integrate(y, x), 1)

C = Collocation([x], [y])
C.SetOrthSys(S)
C.Equation([EQ1])

C.Condition(Eq(y, 0), [0])
C.Condition(Eq(y, 0), [2*pi])
#C.Condition(y==sin(pi)*exp(1), [1])
#C.CollPoints(pnts)

C.setSolver('scipy')
Apprx = C.Solve()
f = .5*exp(-x)*sin(2*x)

G = Graphics('sympy')
G.Plot2D(f, (x, 0, 2*pi), color='blue', legend='Exact')
G.Plot2D(Apprx[0], (x, 0, 2*pi), color='red', legend='Approximation')
G.save('Exm02-%d.png'%(n))