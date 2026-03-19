ElimProblem = new Type of HashTable
ShiftSet = new Type of List
MonomialPartition = new Type of List

-- length ShiftSet := sh -> sum(sh, numcols)

getH0 = method(Options => {MonomialOrder => null})
getH0 (RingElement, Ideal) := o -> (a, J) -> (
    R := ring J;
    B := basis(R/J);
    getH0(a, B, J, o)
)
getH0 (RingElement, Matrix, Ideal) := o -> (a, B, J) -> (
    R := ring J;
    FF := coefficientRing R;
    MO := if not instance(o.MonomialOrder, Nothing) then o.MonomialOrder else (options R).MonomialOrder;
    S := newRing(R, MonomialOrder => MO);
    F := sub(J, S);
    G := gb(F, ChangeMatrix => true);
    aS := sub(a, S);
    BS := sub(B, S);
    P := last coefficients(BS%F); -- change of basis matrix
    Bhat := rsort lift(basis(S/F), S);
    Ta := last coefficients((aS * Bhat)%F, Monomials => Bhat);
    V := aS * BS - lift(aS * sub(BS * (inverse P), S/F), S) * P;
    V = aS * BS - BS * (inverse P) * Ta * P;
    HVG := V // gens G;
    HGF := getChangeMatrix G;
    assert(gens G * HVG - V == 0);
    assert(gens F * HGF - gens G == 0);
    H0 := HGF * HVG;
    sub(H0, ring J)
    )

shiftPolynomials = (shifts, J) -> (
    assert(length shifts == numgens J);
    apply(shifts, J_*, (m, f) -> f * sub(m, ring J))
    )

getTemplate = method(Options => {MonomialOrder => null})
getTemplate(RingElement, Matrix, Ideal) := o -> (a, B, J) -> (
    H0 := getH0(a, B, J, o);
    shifts := new ShiftSet from apply(numgens J, i -> monomials(H0^{i}));
    allMons := union(set \ flatten \ entries \ monomials \ shiftPolynomials(shifts, J));
    monsB := set flatten entries(lift(B, ring J));
--    monsB := intersect(allMons, set flatten entries(lift(B, ring J)));
    print(B, rsort toList monsB);
    monsR := set flatten entries(a * lift(B, ring J)) - set flatten entries(lift(B, ring J));
    print(a*lift(B, ring J), rsort toList monsR);
    monsE := allMons - union(monsR, monsB);
    monomialPartition := new MonomialPartition from rsort \ toList \ {monsE, monsR, monsB};
    (shifts, monomialPartition)
)


getTemplateMatrix = method(Options => {MonomialOrder => null})
getTemplateMatrix(RingElement, Matrix, Ideal) := o -> (a, B, J) -> (
    (shifts, monomialPartition) := getTemplate(a, B, J, o);
    getTemplateMatrix(shifts, monomialPartition, J, o)
    )
getTemplateMatrix(ShiftSet, MonomialPartition, Ideal) := o -> (shifts, monomialPartition, J) -> (
    allMons := apply(fold(monomialPartition, (a,b) -> a|b), m -> sub(m, ring J));
    sub(tra
    \nspose fold(apply(shiftPolynomials(shifts, J), m -> last coefficients(m, Monomials => allMons)), (a,b) -> a|b), coefficientRing ring J)
    )

needsPackage "NumericalLinearAlgebra"
getActionMatrix = (actVar, mp, M) -> (
    a := length mp#0; -- number of "excessive monomials"
    b := length mp#1; -- number of "reducible monomials"
    c := length mp#2; -- number of "basic monomials"
    (m, n) := (numrows M, numcols M);

    -- eliminate "excessive monomials" w/ LU
    Ma := M_{0..a-1};
    (P, L, U) := LUdecomposition Ma;
    L = L | matrix apply(m, i -> apply(m - a, j -> if i == j + a then 1_RR else 0_RR));
    M1 := inverse(id_(RR^m)_P * L) * M;

    -- extract action matrix from reduced and basic monomials in template
    Mr := M1_{a..a+b-1}^{m-b..m-1};
    Mb := M1_{a+b..n-1}^{m-b..m-1};
    A := -solve(Mr, Mb);

    extraMonomials := toList(mp#2 - set apply(mp#2, p -> numerator(p/actVar)));
    -- b = {y^2, y, x, 1};
    binaryMatrix := matrix apply(extraMonomials, m -> apply(mp#2, n -> if m == n then 1_RR else 0_RR));
    A || binaryMatrix
    )


end--
restart
load "EliminationTemplates.m2"

R = QQ[x,y]
-- Example 1
J = ideal(x^2+y^2-1,x^2+x*y+y^2-1)

-- needsPackage "EigenSolver"
-- zeroDimSolve J

actVar = x
B = lift(basis(R/J), R)
(sh, mp) = getTemplate(actVar, B, J)
M = getTemplateMatrix(sh, mp, J)
Ma = getActionMatrix(actVar, mp, M) 

eigenvalues sub(Ma, RR)

-- Example 2: non-standard basis
B = matrix{{x^2, y, 1}}
J = ideal(x^3 + y^2 - 1, x - y - 1)
getH0(x, B, J)
(sh, mp) = getTemplate(x, B, J)
mp#1
M = getTemplateMatrix(x, B, J)
Mred = reducedRowEchelonForm M
Ma = getActionMatrix(actVar, mp, Mred)


-- essential matrix: template extraction
F = ZZ/32003
FF = frac(F[e_(1,1,1)..e_(4,3,3)])
R = FF[x,y,z]
Es = for i from 1 to 4 list matrix for j from 1 to 3 list for k from 1 to 3 list e_(i,j,k)
E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3
J = ideal(E * transpose E * E - (1/2) * trace(E * transpose E) * E, det E);
R0 = F[x,y,z]--,MonomialOrder=>{Weights => apply(3, i -> random(1, 50))}]
J0 = sub(sub(J, flatten flatten for i from 1 to 4 list for j from 1 to 3 list for k from 1 to 3 list e_(i,j,k) => random F), R0);
(sh, mp) = getTemplate(x, basis(R0/J0), J0)

FF = RR[e_(1,1,1)..e_(4,3,3)]
R = FF[x,y,z]
Es = for i from 1 to 4 list matrix for j from 1 to 3 list for k from 1 to 3 list e_(i,j,k)
E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3
J = ideal(E * transpose E * E - (1/2) * trace(E * transpose E) * E, det E);
elapsedTime T = getTemplateMatrix(sh, mp, J);

first mp
(x0,y0,z0)=(1.0,2.0,3.0)
E0 = matrix{{0,-z0,y0},{z0,0,-x0},{-y0,x0,0}}
Ns = apply(3, i -> first QRDecomposition random(RR^3, RR^3))
Ns = Ns | {E0 - sum apply(Ns, {x0,y0,z0}, (a,b) -> a*b)}

M = sub(T, flatten flatten for i from 1 to 4 list for j from 1 to 3 list for k from 1 to 3 list e_(i,j,k) => (Ns#(i-1))_(j-1,k-1))

-- TODO: read off action matrix
