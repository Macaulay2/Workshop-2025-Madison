path = prepend("../", path)
needsPackage "EliminationTemplates"
R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
B = lift(basis(R/J), R)
getH0(x, B, J)
(sh, mp) = getTemplate(x, B, J)
M = getTemplateMatrix(x, B, J)
Ma = getActionMatrix(x, mp, M)
assert(all(sort eigenvalues Ma, {-2,0,1}, (e1, e2) -> abs(e1-e2) < 1e-4))
