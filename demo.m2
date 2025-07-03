restart

-- Not sure if this comment below is still relevant:
----- The ehrhart method doesn't give a warning when it's applied to a non-integral rational polytope, and sometimes its give the wrong polynomial.

needsPackage("RationalPolytopes", FileName => "Durham_2024/RationalPolytopes.m2")

------------ START examples for Friday's demo

---- (I) EHRHART THEORY from package RationalPolytopes

--- (1) Quasipolynomial
-- a polynomial whose coefficients are periodic functions
M1 = matrix{{1,2,3},{1,4,5}}
Q1 = quasiPolynomial M1
period Q1
M2 = matrix{{1,1},{1,2},{1,1},{1,2}}
Q2 = quasiPolynomial M2
period Q2

--- (2) Ehrhart Polynomial
-- a way to count the number of lattice points in a (rational) polytope
-- a (quasi)polynomial in T where T is the level of dilation of the polytope
P = convexHull transpose matrix {{0,0},{1,0},{0,1/2}}
vertices P
#latticePoints(P) -- might be slow on big polytopes
#latticePoints(4*P)

QP = ehrhartQP P -- compute the Ehrhart Polynomial of P 
period QP
QP(0)
QP(4)
QP(4) == #latticePoints(4*P)

-- another example
P2 = convexHull transpose matrix {{-1,0}, {1,0}, {0,1/2}, {0,-1/2}}
vertices P2
QP2 = ehrhartQP P2 -- note period collapse
period QP2 
QP2(7)

--- (3) h*-polynomial and Ehrhart Series
-- Ehrhart series is the generating function for the Ehrhart polynomial
-- The series is equal to a fraction with some "h* polynomial" as numerator
-- and the expression (1-t^period)^(k+1) on the denominator (period = 1 if not rational)
hStarPolynomial P2
hStarPolynomial(P2, Strategy => M2) -- two strategies but both giving same answer

ehrhartSeries P2
ehrhartSeries(P2, Strategy => M2)

--------------------------------

---- (II) EQUIVARIANT EHRHART THEORY from package EquivariantEhrhart
