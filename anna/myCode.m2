--R = QQ[x,y]
--I = ideal(x^2 + y^2 -1, x^2 + x*y + y^2 - 1)
--B = basis(R/I)
R = QQ[x,y]
J = ideal(x^2+y^2-1,x^2+x*y+y^2-1)
B = matrix{{y^2,y,x,1}}
a = x
getH0(a,B,J)


FF = coefficientRing R
S = newRing(R)
F = sub(J,S)
G = gb(F, ChangeMatrix => true)
aS = sub(a,S)
BS = sub(B,S)
coefficients(BS%F)
P = last coefficients(BS%F)
Bhat = rsort lift(basis(S/F),S)
Ta = last coefficients((aS*Bhat)%F, Monomials => Bhat)
V = aS * BS - lift(aS * sub(BS * (inverse P), S/F), S)*P

-*
getHO

input:
ring element a
matrix B
ideal J

output:
matrix H0

*-
