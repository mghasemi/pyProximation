from sympy import *
from Measure import *

x = Symbol('x')
t = Symbol('t')
y = Function('y')(t, x)

n = 6

S = OrthSystem([t, x], [(0, 1), (0, 1)])
B = S.PolyBasis(n)
S.Basis(B)
S.FormBasis()

m = int(sqrt(len(S.OrthBase)))
pnts = [(.2, .4), (.4, .6), (.4, .2), (.6, .2), (.2, .8)]

Z = t*sin(pi*x)
R = diff(Z, t) - diff(Z, x)
EQ1 = Eq(diff(y, t) - diff(y, x), R)
#EQ1 = Eq(diff(y, t) - diff(y, x), R)

#sries = S.Series(Z)
#ChAprx = sum([S.OrthBase[i]*sries[i] for i in range(m)])
#G = plot(ChAprx, (x, -1, 1), color='green')

C = Collocation([t, x], [y])
C.SetOrthSys(S)
C.Equation([EQ1])

C.Condition(Eq(y, 0), [0, 0])
C.Condition(Eq(y, 0), [0, .3])
C.Condition(Eq(y, 0), [1, 1])
C.Condition(Eq(y, 1), [1, .5])
C.Condition(Eq(y, 0), [0, .7])


#C.Condition(Eq(y, 0), [0, 0])
#C.Condition(Eq(y, 0), [0, .3])
#C.Condition(Eq(y, 0), [1, 1])
#C.Condition(Eq(y, 1), [1, .5])
#C.Condition(Eq(y, 0), [0, .7])

#C.CollPoints(pnts)
C.setSolver('scipy')

Apprx = C.Solve()
from sympy.plotting import plot3d
G = plot3d(Apprx[0], (t, 0, 1), (x, 0, 1), adaptive=True, show=False)
#G += plot3d(Z._sympy_(), (ts, 0, 1), (xs, 0, 1))
G.save('Exm3-%d.png'%(n))