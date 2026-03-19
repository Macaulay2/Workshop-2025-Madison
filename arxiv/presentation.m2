--Example 1 (general workflow and copyTemplate)
restart
needsPackage "EliminationTemplates"
R=QQ[x,y]
I=ideal(x^2+y^2-1,x^2+y^3+x*y-2)
B=basis(R/I)
E=eliminationTemplate(x+4*y,I)
getTemplate(E)
getEigenMatrix(E)
sols = templateSolve(E)
assert(all(sols, x -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{x})))
use R
J=ideal(x^2+y^2-2,x^2+y^3+3*x*y-5)
F=copyTemplate(E,J)
getEigenMatrix(F)
sols = templateSolve(F)
assert(all(sols, x -> 1e-6 > norm sub(sub(gens J, CC[gens R]), matrix{x})))

--Example 2 (Help page for package, datatype, and different strategy input)
viewHelp EliminationTemplates
R = QQ[x,y,z]
J = ideal(x^3+y^3+z^3-4,x^2-y-z-1,x-y^2+z-3)
E1 = eliminationTemplate(x, J)
E2 = eliminationTemplate(x, J)
getTemplateMatrix(E1, Strategy => "Greedy")
getTemplateMatrix(E2)
getActionMatrix E1
getActionMatrix E2
eigenvalues getActionMatrix E1
eigenvalues getActionMatrix E2

--Example 3 (Essential Matrix)
--Future direction(s): other minimal problems in computer vision
R = QQ[x,y,z]
Es = apply(4, i -> random(QQ^3, QQ^3))
E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3
I = ideal(E*transpose E * E - (1/2) * trace(E * transpose E) * E)
l = random(1, R)
sols=templateSolve(l, I)
assert(all(sols, x -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{x})))
