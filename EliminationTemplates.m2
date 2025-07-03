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
    "getActionMatrix",
    "templateSolve",
    "EliminationTemplate",
    "eliminationTemplate",
    "shifts",
    "monomialPartition",
    "templateMatrix",
    "actionVariable",
}

EliminationTemplate = new Type of HashTable
ShiftSet = new Type of List
MonomialPartition = new Type of List

eliminationTemplate = method(Options => {})
eliminationTemplate (RingElement, Ideal) := o -> (aVar, J) -> (
    R := ring J;
--    (sh, mp) := getTemplate(aVar, basis(R/J), J);
--    M := getTemplateMatrix(shifts, monomialPartition, J);
    new EliminationTemplate from {
--	    shifts => sh,
--        monomialPartition => mp,
  --      templateMatrix => M,
        "actionVariable" => aVar,
        ideal => J,
	cache => new CacheTable from {}
    }
)

actionVariable = method()
actionVariable EliminationTemplate := E -> E#"actionVariable"

ideal EliminationTemplate := E -> E#ideal


getH0 = method(Options => {MonomialOrder => null, Strategy => null})
getH0 (RingElement, Ideal) := o -> (a, J) -> (
    R := ring J;
    B := basis(R/J);
    H0 := getH0(a, B, J, o);
    if (o.Strategy == null) then H0 else if (o.Strategy == "Larsson") then (
	H0%ideal(syz(H0))
	) else (error "Strategy not yet implemented.") 
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
getTemplate EliminationTemplate := o -> E -> (
    if (E.cache#?"shifts" and E.cache#?"monomialPartition") then (E.cache#"shifts", E.cache#"monomialPartition") else (
	    aVar := actionVariable E;
	    J := ideal E;
	    R := ring J;
	    (sh, mp) := getTemplate(aVar, basis(R/J), J);
	    E.cache#"shifts" = sh;
	    E.cache#"monomialPartition" = mp;
	    (sh, mp)
	    )
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
getTemplateMatrix(EliminationTemplate) := o -> E -> (
    if E.cache#?"templateMatrix" then E.cache#"templateMatrix" else (
	(shifts, monomialPartition) := getTemplate E;
	J := ideal E;
	ret := getTemplateMatrix(shifts, monomialPartition, J);
	E.cache#"templateMatrix" = ret;
	ret
    )
)
    
net EliminationTemplate := E -> (
    str := " action variable: " | toString(actionVariable E);
    if E.cache#?"templateMatrix" then str = "Template matrix:\n" | net(E.cache#"templateMatrix") | str;
    if E.cache#?"actionMatrix" then str = "Action matrix:\n" | net(E.cache#"actionMatrix") | str;
    str
    )

needsPackage "NumericalLinearAlgebra"
getActionMatrix = method(Options => {MonomialOrder => null})
getActionMatrix(RingElement, MonomialPartition, Matrix) := o -> (actVar, mp, M) -> (
    a := length mp#0; -- number of "excessive monomials"
    b := length mp#1; -- number of "reducible monomials"
    c := length mp#2; -- number of "basic monomials"
    (m, n) := (numrows M, numcols M);
    --M1 := reducedRowEchelonForm M;
    
    -- eliminate "excessive monomials" w/ LU
    Ma := M_{0..a-1};
    (P, L, U) := LUdecomposition Ma;
    L = L | matrix apply(m, i -> apply(m - a, j -> if i == j + a then 1_CC else 0_CC));
    M1 := inverse(id_(CC^m)_P * L) * M;

    -- extract action matrix from reduced and basic monomials in template
    Mr := M1_{a..a+b-1}^{m-b..m-1};
    Mb := M1_{a+b..n-1}^{m-b..m-1};
    A := -solve(Mr, Mb);
    
    --A := M1_{a+b..n-1}^{m-b..m-1};
    extraMonomials := rsort toList(mp#2 - set apply(mp#2, p -> numerator(p/actVar)));
    print extraMonomials;
    -- b = {y^2, y, x, 1};
    if #extraMonomials > 0 then (
        binaryMatrix := matrix apply(extraMonomials, m -> apply(mp#2, n -> if m == n then 1_CC else 0_CC));
        A || binaryMatrix
	) else A
)

getActionMatrix(EliminationTemplate) := o -> E -> (
    if E.cache#?"actionMatrix" then E.cache#"actionMatrix" else (
        actVar := actionVariable E;
        (sh, mp) := getTemplate E;
        templateMatrix := getTemplateMatrix E;
        ret := getActionMatrix(actVar, mp, templateMatrix);
	    E.cache#"actionMatrix" = ret;
	    ret
    )
)


templateSolve = method(Options => {MonomialOrder => null, Tolerance => 1e-10})
--templateSolve(EliminationTemplate) := o -> (template) -> ()
templateSolve(Ideal) := o -> (I) -> (

)
templateSolve(RingElement, Ideal) := o -> (a, J) -> (
    R := ring J;
    K := coefficientRing R;
    ringVars := flatten entries vars R;
    R = K[prepend("s", ringVars), MonomialOrder => Eliminate 1];
    I := sub(J, R) + ideal(R_0 - sub(a, R));
    actvar := R_0;
    
    B := lift(basis(R/I), R);
    (sh, mp) := getTemplate(actvar, B, I);
    M := getTemplateMatrix(sh, mp, I);
    Ma := getActionMatrix(actvar, mp, M);
    (svals, P) := eigenvectors Ma;
    (transpose rsort B, clean_(1e-10) (P * inverse diagonalMatrix(P^{numColumns P - 1})))
)

beginDocumentation()

doc ///
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
///

doc ///
 Node
  Key
   getTemplate
   (getTemplate, RingElement, Matrix, Ideal)
   (getTemplate, EliminationTemplate)
  Headline
     extracts a "sparse" representation of an elimination template
  Usage
    (sh, mp) = getTemplate(a, B, J)
  Inputs
    a:RingElement
      the action polynomial defining a multiplication matrix
    B:Matrix
      a basis for a zero-dimensional quotient ring
    J:Ideal
      a zero-dimensional ideal
    MonomialOrder=>Thing
      the monomial order used on the ambient ring
  Outputs
    shifts:ShiftSet
      A list of matrices, each encoding rows of the template matrix
    monomialPartition:MonomialPartition
      A list of monomials encoding columns of the template matrix
  Description
   Text
     This method builds an elimination template. It returns a Sequence of length two, which can be used to recover the template matrix.

     The elements of this sequence encode the rows and columns of a Macaulay matrix (the template matrix.)
     The last element consists of lists of three monomials supported on equations indexing the rows of the template matrix.
     These are called excessive monomials, reducible monomials, and basic monomials.
   Example
     R = QQ[x,y];
     J = ideal(x^2+y^2-1, x^2+x*y+y^2-1);    
     actVar = x;
     B = lift(basis(R/J), R);
     (sh, mp) = getTemplate(actVar, B, J)
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
J = ideal(x^3 + y^2 - 1, x - y - 1)
B = lift(basis(R/J), R)
getH0(x, B, J)
(sh, mp) = getTemplate(x, B, J)
M = getTemplateMatrix(x, B, J)
Mx = getActionMatrix(x, mp, M)
eigenvalues Mx
assert(all(sort eigenvalues Mx, {-2,0,1}, (e1, e2) -> abs(e1-e2) < 1e-4))
-- TODO: add assertion
///

TEST ///
R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
E = eliminationTemplate(x, J)
getTemplateMatrix E
///

end--

-* Development section *-
restart
loadPackage "EliminationTemplates"
needsPackage "NumericalAlgebraicGeometry"
R=QQ[x,y,z]
J=ideal(x^3+y^3+z^3-4,x^2-y-z-1,x-y^2+z-3)
B=basis(R/J)
templateSolve(z,J)
templateSolve(x,J)
templateSolve(x+2*y+3*z,J)
templateSolve(x,J)
netList solveSystem J_*




restart
debug needsPackage "EliminationTemplates"
check "EliminationTemplates"



uninstallPackage "EliminationTemplates"
restart
installPackage "EliminationTemplates"
help EliminationTemplates
help getTemplate

viewHelp "EliminationTemplates"

R = QQ[x,y]
E = eliminationTemplate(x, ideal(x^3 + y^2 - 1, x - y - 1))
getActionMatrix E
net E
