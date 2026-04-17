-- Test Strategy => "MatrixHi" via the public API.
path = prepend("./", path)
needsPackage "EliminationTemplates"

setRandomSeed 42
R = QQ[x,y,z]
Es = apply(4, i -> random(QQ^3, QQ^3))
E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3
I = ideal(E * transpose E * E - (1/2) * trace(E * transpose E) * E) + ideal(det E)

l = y
ET = eliminationTemplate(l, I)

M1 = getTemplateMatrix(ET);
<< "Default:   " << numRows M1 << " x " << numColumns M1 << endl;
M2 = getTemplateMatrix(ET, Strategy => "Larsson");
<< "Larsson:   " << numRows M2 << " x " << numColumns M2 << endl;
M3 = getTemplateMatrix(ET, Strategy => "MatrixHi");
<< "MatrixHi:  " << numRows M3 << " x " << numColumns M3 << endl;
<< "Paper target: 10 x 20" << endl;

<< endl << "=== 3-var benchmark ===" << endl;
R2 = QQ[x,y,z];
J2 = ideal(x^3+y^3+z^3-4, x^2-y-z-1, x-y^2+z-3);
ET2 = eliminationTemplate(x, J2);
N1 = getTemplateMatrix(ET2);
<< "Default:   " << numRows N1 << " x " << numColumns N1 << endl;
N2 = getTemplateMatrix(ET2, Strategy => "Larsson");
<< "Larsson:   " << numRows N2 << " x " << numColumns N2 << endl;
N3 = getTemplateMatrix(ET2, Strategy => "MatrixHi");
<< "MatrixHi:  " << numRows N3 << " x " << numColumns N3 << endl;

exit 0
