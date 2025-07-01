-- -*- coding: utf-8 -*-
newPackage(
    "EliminationTemplates",
    Version => "1",
    Date => "July 1, 2025",
    Authors => {
	{Name => "Manav Batavia",
	    Email => "manavbatavia@gmail.com"},
	{Name => "Cheng Chen",
	    Email => "chengchen@math.wisc.edu"},
	{Name => "Wanchun / Rosie Shen", 
	    Email => "wshen@math.harvard.edu"},
	{Name => "Anna Natalie Chlopecki",
	    Email => "achlopec@purdue.edu"},
	{Name => "Tim Duff", 
	    Email => "tduff@missouri.edu"},
	{Name => "Will Huang", 
	    Email => "williamhuang5120@gmail.com"},
	{Name => "Aolong Li", 
	    Email => "lial0921.miu@gmail.com"},
	{Name => "Ikenna Nometa", 
	    Email => "inometa@hawaii.edu"}	
    },
    Headline => "elimination templates",
    Keywords => {"Documentation"},
    DebuggingMode => false
    )

export {
    "getH0",
    "shiftPolynomials",
    "getTemplate",
    "getTemplateMatrix",
    "getActionMatrix"
}

ElimProblem = new Type of HashTable
ShiftSet = new Type of List
MonomialPartition = new Type of List

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
    sub(transpose fold(apply(shiftPolynomials(shifts, J), m -> last coefficients(m, Monomials => allMons)), (a,b) -> a|b), coefficientRing ring J)
    )

needsPackage "NumericalLinearAlgebra"
getActionMatrix = (actVar, mp, M) -> (
    a := length mp#0; -- number of "excessive monomials"
    b := length mp#1; -- number of "reducible monomials"
    c := length mp#2; -- number of "basic monomials"
    (m, n) := (numrows M, numcols M);
    --M1 := reducedRowEchelonForm M;
    
    -- eliminate "excessive monomials" w/ LU
    Ma := M_{0..a-1};
    (P, L, U) := LUdecomposition Ma;
    L = L | matrix apply(m, i -> apply(m - a, j -> if i == j + a then 1_RR else 0_RR));
    M1 := inverse(id_(RR^m)_P * L) * M;

    -- extract action matrix from reduced and basic monomials in template
    Mr := M1_{a..a+b-1}^{m-b..m-1};
    Mb := M1_{a+b..n-1}^{m-b..m-1};
    A := -solve(Mr, Mb);
    
    --A := M1_{a+b..n-1}^{m-b..m-1};
    extraMonomials := rsort toList(mp#2 - set apply(mp#2, p -> numerator(p/actVar)));
    print extraMonomials;
    -- b = {y^2, y, x, 1};
    if #extraMonomials > 0 then (
        binaryMatrix := matrix apply(extraMonomials, m -> apply(mp#2, n -> if m == n then 1_RR else 0_RR));
        A || binaryMatrix
	) else A
    )

beginDocumentation()

doc /// -- TODO
 Node
  Key
   EliminationTemplates
  Headline
     an example Macaulay2 package
  Description
   Text
    {\em EliminationTemplates} is a basic package to be used as an example.
  Caveat
    Still trying to figure this out.
  Subnodes
    firstFunction
 Node
  Key
   (firstFunction,ZZ)
   firstFunction
  Headline
   a silly first function
  Usage
   firstFunction n
  Inputs
   n:
  Outputs
   :
    a silly string, depending on the value of {\tt n}
  Description
   Text
    Here we show an example.
   Example
    firstFunction 1
    firstFunction 0
///

TEST ///
R = QQ[x,y]
J = ideal(x^2+y^2-1,x^2+x*y+y^2-1)
actVar = x
B = lift(basis(R/J), R)
(sh, mp) = getTemplate(actVar, B, J)
M = getTemplateMatrix(sh, mp, J)
Ma = getActionMatrix(actVar, mp, M) 
eigenvalues Ma
-- TODO: add assertion
///


TEST ///
R = QQ[x,y]
B = matrix{{x^2, y, 1}}
J = ideal(x^3 + y^2 - 1, x - y - 1)
getH0(x, B, J)
(sh, mp) = getTemplate(x, B, J)
M = getTemplateMatrix(x, B, J)
Ma = getActionMatrix(actVar, mp, M)
eigenvalues Ma
-- TODO: add assertion
///

end--


-* Development section *-
restart
debug needsPackage "EliminationTemplates"
check "EliminationTemplates"

uninstallPackage "EliminationTemplates"
restart
installPackage "EliminationTemplates"
viewHelp "EliminationTemplates"
