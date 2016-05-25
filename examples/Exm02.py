from sympy import *
from IntgDiff import *
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

G = plot(f, Apprx[0], (x, 0, 2*pi), show=False)
G[0].line_color = 'blue'
G[1].line_color = 'green'
G.save('Exm02-%d.png'%(n))