path = prepend("../", path)
needsPackage "EliminationTemplates"

R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
B = lift(basis(R/J), R)

-- Greedy should run and produce a valid template matrix.
MDefault = getTemplateMatrix(x, B, J)
MGreedy = getTemplateMatrix(x, B, J, Strategy => "Greedy")
assert(numRows MGreedy <= numRows MDefault)

-- Deterministic behavior across runs.
MGreedy2 = getTemplateMatrix(x, B, J, Strategy => "Greedy")
assert(numRows MGreedy == numRows MGreedy2)
assert(numColumns MGreedy == numColumns MGreedy2)

-- Ensure solving still works in the Greedy path.
E = eliminationTemplate(x, J)
A = getActionMatrix(E, Strategy => "Greedy")
assert(numRows A > 0)
