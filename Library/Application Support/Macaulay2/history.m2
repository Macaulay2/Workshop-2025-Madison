-- This is the beginning of your Macaulay2 log stored at /Users/aolongli/Developments/Workshop-2025-Madison/Library/Application Support/Macaulay2/history.m2
-- This is the beginning of your Macaulay2 log stored at /Users/aolongli/Developments/Workshop-2025-Madison/Library/Application Support/Macaulay2/history.m2
print(random 100);
print(random(0,100));
print(random(QQ));
-- This is the beginning of your Macaulay2 log stored at /Users/aolongli/Developments/Workshop-2025-Madison/Library/Application Support/Macaulay2/history.m2
print(random 100);
print(random(0,100));
print(random(QQ));
setRandomSeed 12345;
print(random 100);
setRandomSeed 12345;
print(random 100);
print(random 100);
print(randomSeed);
setRandomSeed 12345;
print(random 100);
setRandomSeed 12345;
print(random 100);
print(random 100);
print(randomSeed);
load "EliminationTemplates.m2"
R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
B = lift(basis(R/J), R)
MGreedy = getTemplateMatrix(x, B, J, Strategy => "Greedy")
print(numRows MGreedy)
print(random 100);
print(randomSeed);
load "EliminationTemplates.m2"
R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
B = lift(basis(R/J), R)
MGreedy = getTemplateMatrix(x, B, J, Strategy => "Greedy")
print(numRows MGreedy)
load "EliminationTemplates.m2"
R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
B = lift(basis(R/J), R)
M1 = getTemplateMatrix(x, B, J, Strategy => "Greedy")
M2m = getTemplateMatrix(x, B, J, Strategy => "Greedy")
print(numRows M1 == numRows M2m)
print(numColumns M1 == numColumns M2m)
J = ideal(x^3 + y^2 - 1, x - y - 1)
B = lift(basis(R/J), R)
M1 = getTemplateMatrix(x, B, J, Strategy => "Greedy")
M2m = getTemplateMatrix(x, B, J, Strategy => "Greedy")
print(numRows M1 == numRows M2m)
print(numColumns M1 == numColumns M2m)
load "EliminationTemplates.m2"
R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
B = lift(basis(R/J), R)
MGreedy = getTemplateMatrix(x, B, J, Strategy => "Greedy")
print(numRows MGreedy)
MGreedy = getTemplateMatrix(x, B, J, Strategy => "Greedy")
print(numRows MGreedy)
print(hash "abc")
print(toString hash "abc")
MGreedy = getTemplateMatrix(x, B, J, Strategy => "Greedy")
print(numRows MGreedy)
MGreedy = getTemplateMatrix(x, B, J, Strategy => "Greedy")
print(numRows MGreedy)
print(hash "abc")
print(toString hash "abc")
load "EliminationTemplates.m2"
R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
B = lift(basis(R/J), R)
MGreedy = getTemplateMatrix(x, B, J, Strategy => "Greedy")
print(numRows MGreedy)
print(hash "abc")
print(toString hash "abc")
load "EliminationTemplates.m2"
R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
B = lift(basis(R/J), R)
MGreedy = getTemplateMatrix(x, B, J, Strategy => "Greedy")
print(numRows MGreedy)
load "EliminationTemplates.m2"
R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
B = lift(basis(R/J), R)
M1 = getTemplateMatrix(x, B, J, Strategy => "Greedy")
M2m = getTemplateMatrix(x, B, J, Strategy => "Greedy")
print(numRows M1 == numRows M2m)
print(numColumns M1 == numColumns M2m)
