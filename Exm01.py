from sympy import *
from Measure import *
x = Symbol('x')
y = Function('y')(x)

n = 8
M = Measure([(-1, 1)], lambda x:1./sqrt(1.-x**2))
S = OrthSystem([x], [(-1, 1)])
#S.SetMeasure(M)
B = S.PolyBasis(n)
#B = S.FourierBasis(n)
S.Basis(B)
S.FormBasis()

m = len(S.OrthBase)
pnts = [[-1 + i*2./m] for i in range(m)]

Z = sin(pi*x)*exp(x)
R = diff(Z, x, x) - diff(Z, x)
EQ1 = (Eq(diff(y, x, x) - diff(y, x), R))

sries = S.Series(Z)
ChAprx = sum([S.OrthBase[i]*sries[i] for i in range(m)])

C = Collocation([x], [y])
C.SetOrthSys(S)
C.Equation([EQ1])

C.Condition(Eq(y, 0), [0])
C.Condition(Eq(y, sin(-pi)*exp(-1)), [-1])
C.Condition(Eq(y, sin(pi)*exp(1)), [1])
#C.CollPoints(pnts)

C.setSolver('scipy')
Apprx = C.Solve()
print Apprx[0]
G = plot(Z, ChAprx, Apprx[0], (x, -1, 1), show=False)
G[0].line_color = 'blue'
G[1].line_color = 'green'
G[2].line_color = 'red'
G.save('Exm01-%d.png'%(n))