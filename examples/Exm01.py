Symbolic = 'sympy'

if Symbolic == 'sympy':
	from sympy import *
	x = Symbol('x')
	y = Function('y')(x)
elif Symbolic == 'sage':
	from sage.all import *
	x = var('x')
	y = function('y')(x)

from ApproxPy import *

n = 6
M = Measure([(-1, 1)], lambda x:1./sqrt(1.-x**2))
S = OrthSystem([x], [(-1, 1)], Symbolic)
S.SetMeasure(M)
B = S.PolyBasis(n)
#B = S.FourierBasis(n)
S.Basis(B)
S.FormBasis()

m = len(S.OrthBase)
pnts = [[-1 + i*2./m] for i in range(m)]

Z = sin(pi*x)*exp(x)
R = diff(Z, x, x) - diff(Z, x)

if Symbolic == 'sympy':
	EQ1 = (Eq(diff(y, x, x) - diff(y, x), R))
elif Symbolic == 'sage':
	EQ1 = (diff(y, x, x) - diff(y, x) == R)

sries = S.Series(Z)
ChAprx = sum([S.OrthBase[i]*sries[i] for i in range(m)])

C = Collocation([x], [y], Symbolic)
C.SetOrthSys(S)
C.Equation([EQ1])

if Symbolic == 'sympy':
	C.Condition(Eq(y, 0), [0])
	C.Condition(Eq(y, sin(-pi)*exp(-1)), [-1])
	C.Condition(Eq(y, sin(pi)*exp(1)), [1])
elif Symbolic == 'sage':
	C.Condition(y == 0, [0])
	C.Condition(y == sin(-pi)*exp(-1), [-1])
	C.Condition(y == sin(pi)*exp(1), [1])

#C.CollPoints(pnts)

C.setSolver('scipy')
Apprx = C.Solve()
print Apprx[0]
G = Graphics(Symbolic)
G.Plot2D(Z, (x, -1, 1), color='blue', legend='Exact')
G.Plot2D(ChAprx, (x, -1, 1), color='green', legend='Orth Apprx')
G.Plot2D(Apprx[0], (x, -1, 1), color='red', legend='Colloc Apprx')

G.save('Exm01-%d.png'%(n))