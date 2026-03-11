-- -*- coding: utf-8 -*-
newPackage(
    "EliminationTemplates",
    Version => "1",
    Date => "July 1, 2025",
    Authors => {
    {Name => "Manav Batavia",
    Email => "manavbatavia@gmail.com",
    HomePage => ""},
    {Name => "Cheng Chen",
    Email => "chengchen@math.wisc.edu",
    HomePage => ""},
    {Name => "Wanchun / Rosie Shen", 
    Email => "wshen@math.harvard.edu",
    HomePage => ""},
    {Name => "Anna Natalie Chlopecki",
    Email => "achlopec@purdue.edu",
    HomePage => ""},
    {Name => "Tim Duff", 
    Email => "tduff@missouri.edu",
    HomePage => "https://timduff35.github.io/timduff35/"},
    {Name => "Will Huang", 
    Email => "williamhuang5120@gmail.com",
    HomePage => ""},
    {Name => "Aolong Li", 
    Email => "lial0921.miu@gmail.com",
    HomePage => ""},
    {Name => "Ikenna Nometa", 
    Email => "inometa@hawaii.edu",
    HomePage => ""}	
    },
    Headline => "elimination templates",
    PackageImports => {"EigenSolver", "NumericalAlgebraicGeometry"},
    Keywords => {"Documentation"},
    HomePage => "",
    DebuggingMode => false,
    AuxiliaryFiles => true
)

-- Greedy-specific helper routines are maintained separately for readability.
load "./EliminationTemplates/GreedyHelpers.m2"

export {
    "getH0",
    "shiftPolynomials",
    "getTemplate",
    "getTemplateMatrix",
    "getActionMatrix",
    "getEigenMatrix",
    "templateSolve",
    "EliminationTemplate",
    "eliminationTemplate",
    "shifts",
    "monomialPartition",
    "templateMatrix",
    "actionVariable",
    "copyTemplate"
}

EliminationTemplate = new Type of HashTable
ShiftSet = new Type of List
MonomialPartition = new Type of List

eliminationTemplate = method(Options => {})
eliminationTemplate (RingElement, Ideal) := o -> (aVar, J) -> (
    R := ring J;
    -- (sh, mp) := getTemplate(aVar, basis(R/J), J);
    -- M := getTemplateMatrix(shifts, monomialPartition, J);
    new EliminationTemplate from {
        -- shifts => sh,
        -- monomialPartition => mp,
        -- templateMatrix => M,
        "actionVariable" => aVar,
        ideal => J,
	      cache => new CacheTable from {}
    }
)

net EliminationTemplate := E -> (
    str := " action variable: " | toString(actionVariable E);
    if E.cache#?"templateMatrix" then str = "Template matrix:\n" | net(E.cache#"templateMatrix") | str;
    if E.cache#?"actionMatrix" then str = "Action matrix:\n" | net(E.cache#"actionMatrix") | str;
    str
)

copyTemplate = method(Options => {})
copyTemplate(EliminationTemplate, Ideal) := o -> (E, J) -> (
    Rnew := ring J;
    
    -- 1. Safely substitute the action variable into the new ring
    aNew := sub(actionVariable E, Rnew);
    F := eliminationTemplate(aNew, J);
    
    -- 2. Promote and copy the basis
    F.cache#basis = sub(basis E, Rnew);
    
    -- 3. Promote and copy the offline structure (if it was computed)
    if E.cache#?"shifts" then 
        F.cache#"shifts" = apply(E.cache#"shifts", sh -> sub(sh, Rnew));
        
    if E.cache#?"monomialPartition" then 
        F.cache#"monomialPartition" = apply(E.cache#"monomialPartition", m -> sub(m, Rnew));
        
    if E.cache#?"graphIdeal" then
        F.cache#"graphIdeal" = sub(E.cache#"graphIdeal", Rnew);

    -- NOTE: We intentionally do NOT copy templateMatrix or actionMatrix.
    -- Those depend on the specific coefficients of J. By copying the shifts,
    -- getTemplateMatrix(F) will instantly build the new matrix without calling getH0.
    
    F
)
-*
copyTemplate(EliminationTemplate, Ideal) := o -> (E, J) -> (
    F := eliminationTemplate(E#"actionVariable", J);
    -- anything else to copy??
    F.cache#basis = basis E;
    F
)
*-

actionVariable = method()
actionVariable EliminationTemplate := E -> E#"actionVariable"
ideal EliminationTemplate := E -> E#ideal
basis EliminationTemplate := o -> E -> E.cache#basis

getH0 = method(Options => {MonomialOrder => null, Strategy => null})
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
    V := aS * BS - lift(aS * sub(BS * (inverse P), S/F), S) * P;
    HVG := V // gens G;
    HGF := getChangeMatrix G;
    assert(gens G * HVG - V == 0);
    assert(gens F * HGF - gens G == 0);
    H0 := HGF * HVG;
    H0 = sub(H0, ring J);

    if (o.Strategy === null) then (
--        print("Using default strategy to compute H0.");
        H0
    )
    else if (o.Strategy == "Greedy") then (
        print("Using Greedy strategy to compute H0.");
        -- compute H = H0^T + Theta*H1, where H1 is the transposed syzygy matrix.
        -- We keep columns aligned with generators of J for downstream greedy helpers.
        H1 := transpose sub(syz(gens J), ring J);
        -- print("H0: " | toString H0);
        -- print("H1: " | toString H1);

        -- create an extension ring of R with the theta variables
        ThetaExt := R[apply(numcols H0 * numrows H1, i -> "t" | toString i)];

        -- coerce H0 and H1 into the extension ring
        toTheta := map(ThetaExt, R);
        H0e := transpose(toTheta H0);
        H1e := toTheta H1;

        -- build Theta over the extension ring
        Theta := genericMatrix(ThetaExt, ThetaExt_0, numcols H0, numrows H1);

        -- now everything is in the same ring
        H := H0e + Theta * H1e;

        data := monomialVectorAndWData(H);
        W := data#"W";
        columnInfo := data#"columnInfo";
        SH := ring H;
        baseR := coefficientRing SH;
        allVars := flatten entries vars SH;
        baseVars := flatten entries vars baseR;
        thetaVars := drop(allVars, #baseVars);
        thetaToZeroMap := map(SH, SH, baseVars | apply(thetaVars, t -> 0_SH));

        rowA := rowWiseGreedyAssignments(W, thetaVars, thetaToZeroMap, SH, baseR);
        excessiveMons := computeExcessiveMonomials(a, B, J, columnInfo, baseR);
        colA := columnWiseGreedyAssignments(W, excessiveMons, columnInfo, J, thetaVars, thetaToZeroMap, SH, baseR);

        rowZero := countZeroColumns(W, rowA, SH);
        colZero := countZeroColumns(W, colA, SH);
        bestA := if colZero > rowZero then colA else rowA;
        bestName := if colZero > rowZero then "Column-wise" else "Row-wise";
        print("Greedy selected: " | bestName | " strategy.");
        -- print("Zero columns in W (row-wise): " | toString rowZero);
        -- print("Zero columns in W (column-wise): " | toString colZero);
        instantiateHWithAssignments(H, bestA, SH, baseR)
    )
    else if (o.Strategy == "Larsson") then (
        print("Using Larsson's strategy to compute H0.");
        -- Ensure the reduction stays in the base ring R
        H0res := H0 % image(syz(gens(J)));
        sub(H0res, ring J)
    )
    else (error "Strategy not yet implemented.") 
)

shiftPolynomials = (shifts, J) -> (
    assert(length shifts == numgens J);
    apply(shifts, J_*, (m, f) -> f * sub(m, ring J))
)

getTemplate = method(Options => {MonomialOrder => null, Strategy => null})
getTemplate(RingElement, Matrix, Ideal) := o -> (a, B, J) -> (
    H0 := getH0(a, B, J, o);
    shifts := new ShiftSet from apply(numgens J, i -> monomials(H0^{i}));
    allMons := union(set \ flatten \ entries \ monomials \ shiftPolynomials(shifts, J));
    if (allMons == set {}) then error "allMons is empty!";
    monsB := set flatten entries(lift(B, ring J));
    monsR := set flatten entries(a * lift(B, ring J)) - set flatten entries(lift(B, ring J));
    monsE := allMons - union(monsR, monsB);
    monomialPartition := new MonomialPartition from rsort \ toList \ {monsE, monsR, monsB};
    (shifts, monomialPartition)
)
getTemplate(EliminationTemplate) := o -> E -> (
    J := ideal E;
    a := actionVariable E;
    R := ring J;
    B := lift(basis(R/J), R);

    -- 1. Compute shifts for the original ideal in R
    (shOrig, mpOrig) := getTemplate(a, B, J, o);

    -- 2. Set up the extended ring Rs
    K := coefficientRing R;
    ringVars := flatten entries vars R;
    MO := if not instance(o.MonomialOrder, Nothing) then o.MonomialOrder else (options R).MonomialOrder;
    Rs := K[prepend("s", ringVars), MonomialOrder => {Eliminate 1, MO}];
    
    aS := sub(a, Rs); -- This is the 'x' that lives in Rs
    Js := sub(J, Rs);
    actVar := Rs_0; -- This is 's'
    Is := Js + ideal(actVar - aS);
    Bs := sub(B, Rs);
    
    sortedBs := rsort flatten entries Bs;
    
    -- 3. Construct shifts for the graph ideal directly
    shiftsGraph := new ShiftSet from (
        apply(shOrig, sh -> sub(sh, Rs)) | {matrix {sortedBs}}
    );

    E.cache#basis = Bs;
    E.cache#"graphIdeal" = Is;

    -- 4. Reconstruct the monomial partition in Rs
    -- WE USE aS HERE instead of a to keep everything in the same ring!
    allMons := union(set \ flatten \ entries \ monomials \ shiftPolynomials(shiftsGraph, Is));
    monsB := set sortedBs;
    monsR := set apply(sortedBs, b -> actVar * b);
    monsE := allMons - union(monsR, monsB);
    mpGraph := new MonomialPartition from rsort \ toList \ {monsE, monsR, monsB};

    (shiftsGraph, mpGraph)
)
-*
getTemplate(EliminationTemplate) := o -> E -> (
    J := ideal E;
    a := actionVariable E;

    -- 1. Compute shifts for the original ideal to avoid graph ideal GB bloat
    R := ring J;
    B := lift(basis(R/J), R);
    (shOrig, mpOrig) := getTemplate(a, B, J, o);

    -- 2. Set up the extended ring with the 's' variable
    K := coefficientRing R;
    ringVars := flatten entries vars R;
    MO := if not instance(o.MonomialOrder, Nothing) then o.MonomialOrder else (options R).MonomialOrder;
    Rs := K[prepend("s", ringVars), MonomialOrder => {Eliminate 1, MO}];
    
    aS := sub(a, Rs);
    Js := sub(J, Rs);
    actVar := Rs_0;
    Is := Js + ideal(actVar - aS);
    Bs := sub(B, Rs);
    
    -- Sort the basis descending to match the partition column order!
    sortedBs := rsort flatten entries Bs;
    
    -- 3. Construct shifts for the graph ideal directly
    shiftsGraph := new ShiftSet from (
        apply(shOrig, sh -> sub(sh, Rs)) | {matrix {sortedBs}}
    );

    E.cache#basis = Bs;
    E.cache#"graphIdeal" = Is;

    -- 4. Reconstruct the monomial partition in the extended ring
    allMons := union(set \ flatten \ entries \ monomials \ shiftPolynomials(shiftsGraph, Is));
    monsB := set sortedBs;
    monsR := set apply(sortedBs, b -> actVar * b);
    monsE := allMons - union(monsR, monsB);
    mpGraph := new MonomialPartition from rsort \ toList \ {monsE, monsR, monsB};

    (shiftsGraph, mpGraph)
)
getTemplate(EliminationTemplate) := o -> E -> (
    J := ideal E;
    a := actionVariable E;
    -- 1. Compute shifts for the original ideal without 's' variable (avoids template matrix bloat)
    R := ring J;
    B := lift(basis(R/J), R);
    (shOrig, mpOrig) := getTemplate(a, B, J, o);
    -- 2. Set up the extended ring with the 's' variable
    K := coefficientRing R;
    ringVars := flatten entries vars R;
    MO := if not instance(o.MonomialOrder, Nothing) then o.MonomialOrder else (options R).MonomialOrder;
    Rs := K[prepend("s", ringVars), MonomialOrder => {Eliminate 1, MO}];
    aS := sub(a, Rs);
    Js := sub(J, Rs);
    actVar := Rs_0;
    Is := Js + ideal(actVar - aS);
    Bs := sub(B, Rs);
    -- 3. Construct shifts for the graph ideal Is directly
    -- The first generators use the original shifts.
    -- The last generator (s - a) is shifted perfectly by the basis elements Bs.
    shiftsGraph := new ShiftSet from (
        apply(shOrig, sh -> sub(sh, Rs)) | {matrix {flatten entries Bs}}
    );
    E.cache#basis = Bs;
    E.cache#"graphIdeal" = Is;
    -- 4. Construct the monomial partition in the extended ring
    allMons := union(set \ flatten \ entries \ monomials \ shiftPolynomials(shiftsGraph, Is));
    monsB := set flatten entries Bs;
    monsR := set flatten entries(actVar * Bs);
    monsE := allMons - union(monsR, monsB);
    mpGraph := new MonomialPartition from rsort \ toList \ {monsE, monsR, monsB};
    (shiftsGraph, mpGraph)
)
getTemplate(EliminationTemplate) := o -> E -> (
    J := ideal E;
    a := actionVariable E;

    R := ring J;
    K := coefficientRing R;
    ringVars := flatten entries vars R;
    MO := if not instance(o.MonomialOrder, Nothing) then o.MonomialOrder else (options R).MonomialOrder;
    R = K[prepend("s", ringVars), MonomialOrder => {Eliminate 1, MO}];
    I := sub(J, R) + ideal(R_0 - sub(a, R));
    actVar := R_0;

    B := lift(basis(R/I), R);
    E.cache#basis = B;
    E.cache#"graphIdeal" = I;
    getTemplate(actVar, B, I, o)
)
*-

getTemplateMatrix = method(Options => {MonomialOrder => null, Strategy => null})
getTemplateMatrix(RingElement, Matrix, Ideal) := o -> (a, B, J) -> (
    (shifts, monomialPartition) := getTemplate(a, B, J, o);
    getTemplateMatrix(shifts, monomialPartition, J, o)
)
getTemplateMatrix(ShiftSet, MonomialPartition, Ideal) := o -> (shifts, monomialPartition, J) -> (
    allMons := apply(fold(monomialPartition, (a,b) -> a|b), m -> sub(m, ring J));
    sub(transpose fold(apply(shiftPolynomials(shifts, J), m -> last coefficients(m, Monomials => allMons)), (a,b) -> a|b), coefficientRing ring J)
)
getTemplateMatrix(EliminationTemplate) := o -> E -> (
    if (E.cache#?"lastTemplateStrategy" === o.Strategy) and (E.cache#?"templateMatrix") then E.cache#"templateMatrix" else (
        (shifts, monomialPartition) := getTemplate(E, o);
        ret := getTemplateMatrix(shifts, monomialPartition, E.cache#"graphIdeal", o);
        E.cache#"templateMatrix" = ret;
        E.cache#"lastTemplateStrategy" = o.Strategy;
        ret
    )
)

getActionMatrix = method(Options => {MonomialOrder => null, Strategy => null})
getActionMatrix(RingElement, MonomialPartition, Matrix) := o -> (actVar, mp, M) -> (
    a := length mp#0; -- number of "excessive monomials"
    b := length mp#1; -- number of "reducible monomials"
    c := length mp#2; -- number of "basic monomials"
    (m, n) := (numrows M, numcols M);
    
    -- eliminate "excessive monomials" w/ LU
    Ma := M_{0..a-1};
    (P, L, U) := LUdecomposition Ma;
    L = L | matrix apply(m, i -> apply(m - a, j -> if i == j + a then 1_CC else 0_CC));
    M1 := solve(id_(CC^m)_P * L, sub(M, CC));

    -- extract action matrix from reduced and basic monomials in template
    Mr := M1_{a..a+b-1}^{m-b..m-1};
    Mb := M1_{a+b..n-1}^{m-b..m-1};
    A := -solve(Mr, Mb);
    
    extraMonomials := rsort toList(mp#2 - set apply(mp#2, p -> numerator(p/actVar)));
    if #extraMonomials > 0 then (
        binaryMatrix := matrix apply(extraMonomials, m -> apply(mp#2, n -> if m == n then 1_CC else 0_CC));
        A || binaryMatrix
	  ) else A
)
getActionMatrix(RingElement, MonomialPartition, Matrix) := o -> (actVar, mp, M) -> (
    numE := length mp#0;
    numR := length mp#1;
    numB := length mp#2;
    m := numrows M;
    n := numcols M;

    -- The bottom numR rows are exactly the shifts of s - a.
    -- The top rows belong to the original ideal J.
    numTop := m - numR;

    if numE == 0 then (
        -- If there are no excessive monomials, the action matrix is just the basic block.
        -M_{numE+numR .. n-1}^{numTop .. m-1}
    ) else (
        -- Slice the top block (original ideal shifts)
        MtopE := M_{0 .. numE-1}^{0 .. numTop-1};
        MtopB := M_{numE+numR .. n-1}^{0 .. numTop-1};

        -- Slice the bottom block (s - a shifts)
        MbotE := M_{0 .. numE-1}^{numTop .. m-1};
        MbotB := M_{numE+numR .. n-1}^{numTop .. m-1};

        -- solve(A,B) finds X such that A*X = B. 
        -- Action Matrix A = MbotE * X - MbotB
        MbotE * solve(MtopE, MtopB) - MbotB
    )
)
getActionMatrix(EliminationTemplate) := o -> E -> (
    if (E.cache#?"lastActionStrategy" === o.Strategy) and E.cache#?"actionMatrix" then E.cache#"actionMatrix" else (
        (sh, mp) := getTemplate(E, o);
        templateMatrix := getTemplateMatrix(E, o);
        
        -- The action variable in the extended ring is always the first variable (Rs_0)
        Rs := ring first first mp;
        actVar := Rs_0;
        
        ret := getActionMatrix(actVar, mp, templateMatrix, o);
        E.cache#"actionMatrix" = ret;
        E.cache#"lastActionStrategy" = o.Strategy;
        ret
    )
)
-*
getActionMatrix(EliminationTemplate) := o -> E -> (
    if (E.cache#?"lastActionStrategy" === o.Strategy) and E.cache#?"actionMatrix" then E.cache#"actionMatrix" else (
        (sh, mp) := getTemplate E;
        templateMatrix := getTemplateMatrix E;
        R := ring first first mp;
        actVar := R_0;
        ret := getActionMatrix(actVar, mp, templateMatrix);
	      E.cache#"actionMatrix" = ret;
        E.cache#"lastActionStrategy" = o.Strategy;
	      ret
    )
)
*-

getEigenMatrix = method(Options => {MonomialOrder => null})
getEigenMatrix(EliminationTemplate) := o -> (E) -> (
    Ma := getActionMatrix(E);
    (svals, P) := eigenvectors Ma;
    cleanEvecs := clean_(1e-10) (P * inverse diagonalMatrix(P^{numColumns P - 1}));

    (transpose rsort basis E, cleanEvecs)
)

getEigenMatrix(Ideal) := o -> (I) -> getEigenMatrix(random(1, ring I), I, o)
getEigenMatrix(RingElement, Ideal) := o -> (a, J) -> (
    E := eliminationTemplate(a, J);
    getEigenMatrix(E, o)
)
-*
getEigenMatrix(RingElement, Ideal) := o -> (a, J) -> (
    R := ring J;
    K := coefficientRing R;
    ringVars := flatten entries vars R;
    MO := if not instance(o.MonomialOrder, Nothing) then o.MonomialOrder else (options R).MonomialOrder;
    R = K[prepend("s", ringVars), MonomialOrder => {Eliminate 1, MO}];
    I := sub(J, R) + ideal(R_0 - sub(a, R));
    actvar := R_0;
    
    B := lift(basis(R/I), R);
    (sh, mp) := getTemplate(actvar, B, I, o);
    M := getTemplateMatrix(sh, mp, I, o);
    Ma := getActionMatrix(actvar, mp, M, o);
    (svals, P) := eigenvectors Ma;
    cleanEvecs := clean_(1e-10) (P * inverse diagonalMatrix(P^{numColumns P - 1}));

    (transpose rsort B, cleanEvecs)
)
*-

templateSolve = method(Options => {MonomialOrder => null})
templateSolve(EliminationTemplate) := o -> (E) -> (
    (B, M) := getEigenMatrix(E, o);
    recoverSolutions(B, M, ideal E)
)
templateSolve(Ideal) := o -> (I) -> templateSolve(random(1,ring I), I, o)
-*
templateSolve(RingElement, Ideal) := o -> (a, J) -> (
    (B, M) := getEigenMatrix(a, J, o);
    recoverSolutions(B, M, J)
)
*-
templateSolve(RingElement, Ideal) := o -> (a, J) -> (
    E := eliminationTemplate(a, J);
    templateSolve(E, o)
)


-- Helper method, don't export
recoverSolutions = method()
recoverSolutions(Matrix, Matrix, Ideal) := (B, M, J) -> (
    basisMons := apply(flatten entries B, m -> sub(m, ring J));
    solutions := {};
    varsList := flatten entries vars ring J;

    for rootIndex from 0 to numColumns M - 1 do (
        monomialValues := new MutableHashTable;
        for i from 0 to #basisMons - 1 do (
            m := basisMons#i;
            monomialValues#m = M_(i, rootIndex);
        );

        -- TODO: check excessive polynomials first before falling back to Groebner basis

        root := {};
        for v in varsList do (
            -- If the variable is in the basis, use directly
            if monomialValues#?v then (
              root = append(root, monomialValues#v);
            )
            else (
                -- Otherwise, reduce modulo ideal: x_i mod J
                r := sub(v % J, ring J);

                -- write r as linear combination of basis monomials
                coeffs := last coefficients(r, Monomials => basisMons);
                value := 0;
                for i from 0 to #basisMons - 1 do (
                    m := basisMons#i;
                    if monomialValues#?m then (
                      value += sub(coeffs_(i,0), coefficientRing ring J) * monomialValues#m;
                    )
                );
                root = append(root, value);
            );
        );
        solutions = append(solutions, root);
    );
    solutions
)

beginDocumentation()

doc ///
 Node
  Key
    EliminationTemplates
  Headline
     zero-dimensional polynomial solvers based on linear algebra
  Description
   Text
    {\em EliminationTemplates} is a package that supports solvers for the following problem: given a zero-dimensional radical ideal $I \subset R := \mathbb{C} [x_1, \ldots , x_n]$, find approximate values for the isolated solutions $(p_1, \ldots , p_n ) \in V_{\mathbb{C}} (I).$
    The main applications occur when the ideal $I$ occur in a parametric family of problems with similar structure.

    Following the references below, the package is geared twowards implementing a "two-stage" approach, consisting of (1) offline stage and (2) an online stage.

    In the offline stage, the structure of a "template matrix" for $I$ is determined using Groebner basis computations.

    In the online stage, prior knowledge of the template matrix can be used to construct a multiplication matrix for the quotient ring $R/I.$
    From this multiplication matrix, solutions can be extracted using eigenvector methods, such as the ones implemented in the package @TO EigenSolver@.
  Caveat
    This package is a work-in-progress!
  References
    @UL {
	{"Optimizing Elimination Templates by Greedy Parameter Search, Martyushev-Vrablikova-Pajdla", EM "CVPR 2022"},
	{"Efficient solvers for minimal problems by syzygy-based reduction, Larsson-Oskarsson-Astrom", EM "CVPR 2017"}
	}@
///

doc ///
 Node
    Key
        EliminationTemplate
    Headline
        type for elimination template objects
    Description
        Text
            The type `EliminationTemplate` represents objects that store the data required for elimination template computations.
            An `EliminationTemplate` object encodes the action variable, the ideal, and a cache for storing computed template data such as shifts, monomial partitions, and template matrices.
            These objects are constructed using the function `eliminationTemplate` and are used as input to other functions in this package, such as `getTemplateMatrix`, `getActionMatrix`, and `templateSolve`.
        Example
          R = QQ[x,y]
          J = ideal(x^2+y^2-1, x^2+x*y+y^2-1)
          E = eliminationTemplate(x, J)
          E
///

doc ///
 Node
    Key
      eliminationTemplate
      (eliminationTemplate, RingElement, Ideal)
    Headline
      constructor for an EliminationTemplate object
    Usage
      E = eliminationTemplate(a, J)
    Inputs
      a:RingElement
        the action polynomial defining a multiplication matrix
      J:Ideal
        a zero-dimensional ideal
    Outputs
      E:EliminationTemplate
        an EliminationTemplate object encoding the data for elimination template computations
    Description
      Text
        This function constructs an EliminationTemplate object, which stores the action variable and ideal, and provides a cache for storing computed template data.
        The EliminationTemplate object can be used with other functions in this package to compute template matrices, action matrices, and solve polynomial systems.
      Example
        R = QQ[x,y]
        J = ideal(x^2+y^2-1, x^2+x*y+y^2-1)
        E = eliminationTemplate(x, J)
///

--Do not want to export this function potentially?
--doc ///
 --Node
  --Key
    --getTemplate
    --(getTemplate, RingElement, Matrix, Ideal)
    --(getTemplate, EliminationTemplate)
  --Headline
    --extracts a "sparse" representation of an elimination template
  --Usage
    --(sh, mp) = getTemplate(a, B, J)
  --Inputs
    --a:RingElement
      --the action polynomial defining a multiplication matrix
    --B:Matrix
      --a basis for a zero-dimensional quotient ring
    --J:Ideal
      --a zero-dimensional ideal
  --Outputs
    --shifts:ShiftSet
      --A list of matrices, each encoding rows of the template matrix
    --monomialPartition:MonomialPartition
      --A list of monomials encoding columns of the template matrix
  --Description
    --Text
      --This method builds an elimination template. It returns a Sequence of length two, which can be used to recover the template matrix.

      --The elements of this sequence encode the rows and columns of a Macaulay matrix (the template matrix.)
      --The last element consists of lists of three monomials supported on equations indexing the rows of the template matrix.
      --These are called excessive monomials, reducible monomials, and basic monomials.
    --Example
      --R = QQ[x,y];
      --J = ideal(x^2+y^2-1, x^2+x*y+y^2-1);    
      --actVar = x;
      --B = lift(basis(R/J), R);
      --(sh, mp) = getTemplate(actVar, B, J)
--///

--doc ///
 --Node
  --Key
--viewH   --[getTemplate, MonomialOrder]
  --Headline
    --the monomial order used on the ambient ring, 
  --Usage
    --getTemplate(a, B, J, MonomialOrder => Eliminate 1)
  --Description
    --Text
      --The monomial order used on the ambient ring. This is used to determine the ordering of the columns of the template matrix.
      --The default is `Eliminate 1`, which is a monomial order that eliminates the first variable.
      --Other monomial orders can be used, such as `Eliminate 2` or `Eliminate 3`.
      --See the documentation for `Macaulay2` for more information on monomial orders.
--///

doc ///
 Node
  Key
    templateSolve
    (templateSolve, EliminationTemplate)
    (templateSolve, Ideal)
    (templateSolve, RingElement, Ideal)
  Headline
    polynomial system solver using elimination templates
  Usage
    (B, ev) = templateSolve(et)
    (B, ev) = templateSolve(J)
    (B, ev) = templateSolve(a, J)
  Inputs
    a:RingElement
      the action polynomial defining a multiplication matrix
    J:Ideal
      a zero-dimensional ideal
    et:EliminationTemplate
      the elimination template for this problem
    MonomialOrder=>Thing
      the monomial order used on the ambient ring
  Outputs
    B:Matrix
      a column matrix containing the basis for R/J used
    ev:Matrix
      a matrix whose columns are the eigenvectors of the action matrix
  Description
   Text
      In the example below, the ideal $J$ defines a zero-dimensional variety with four points.
      This method finds numerical approximations to these four points by solving an eigenvalue problem, much like the package @TO EigenSolver@.
      The main difference between this package and ours is that, for ours, the internal template matrix may be reused for problems of a "similar structure."
   Example
      R = QQ[x,y]
      J = ideal(x^2+y^2-1,x^2+x*y+y^2-1)
      actVar = x + 2*y
      templateSolve(actVar, J)
///

doc ///
 Node
    Key
      copyTemplate
      (copyTemplate, EliminationTemplate, Ideal)
    Headline
      copies EliminationTemplate object
    Usage
      F = copyTemplate(E, J)
    Inputs
      E:EliminationTemplate
        the elimination template to copy
      J:Ideal
        a zero-dimensional ideal
    Outputs
      F:EliminationTemplate
        an EliminationTemplate object
    --Description
      --Text
      --Example
    SeeAlso
      eliminationTemplate
///

doc ///
 Node
    Key
      getTemplateMatrix
      (getTemplateMatrix, RingElement, Matrix, Ideal)
      (getTemplateMatrix, ShiftSet, MonomialPartition, Ideal)
      (getTemplateMatrix, EliminationTemplate)
    --Headline
    --Usage
    --Inputs
    --Outputs
    --Description
      --Text
      --Example
    --SeeAlso
///

doc ///
 Node
    Key
      actionVariable
    --Headline
    --Usage
    --Inputs
    --Outputs
    --Description
      --Text
      --Example
    --SeeAlso
///

TEST ///
  R = QQ[x,y]
  J = ideal(x^2+y^2-1,x^2+x*y+y^2-1)
  actVar = x
  E = eliminationTemplate(actVar, J)
  M = getTemplateMatrix(E)
  Ma = getActionMatrix(E)
  evals = eigenvalues Ma
  assert(all(sort evals, {-1,0,0,1}, (e1, e2) -> abs(e1 - e2) < 1e-4))
///

TEST ///
  R = QQ[x,y]
  J = ideal(x^3 + y^2 - 1, x - y - 1)
  E = eliminationTemplate(x, J)
  Mx = getActionMatrix(E)
  evals = eigenvalues Mx
  assert(all(sort evals, {-2,0,1}, (e1, e2) -> abs(e1 - e2) < 1e-4))
///

TEST ///
  R = QQ[x,y,z]
  J = ideal(x^3+y^3+z^3-4,x^2-y-z-1,x-y^2+z-3)
  E = eliminationTemplate(x, J)
  --H0 = getH0(x,J,Strategy=>"Larsson")
  --H0 = getH0(x,J,Strategy=> null)
  getTemplateMatrix E
  -- getTemplateMatrix(E, Strategy => "Greedy")
  getActionMatrix E
  eigenvalues getActionMatrix E
///

-*TEST ///
  R = QQ[x,y,z]
  J = ideal(x^3+y^3+z^3-4,x^2-y-z-1,x-y^2+z-3)
  E1 = eliminationTemplate(x, J)
  E2 = eliminationTemplate(x, J)
  E3 = eliminationTemplate(x, J)
  --H0 = getH0(x,J,Strategy=>"Larsson")
  --H0 = getH0(x,J,Strategy=> null)
  getTemplateMatrix(E1, Strategy => null)
  getActionMatrix E1
  eigenvalues getActionMatrix E1
  getTemplateMatrix(E2, Strategy => "Larsson")
  getActionMatrix E2
  eigenvalues getActionMatrix E2
  getTemplateMatrix(E3, Strategy => "Greedy")
  getActionMatrix E3
  eigenvalues getActionMatrix E3
///
*-
TEST ///
  R = QQ[x,y,z]
  J = ideal(x^3+y^3+z^3-4,x^2-y-z-1,x-y^2+z-3)

  -- 3 templates, 3 strategies
  E1 = eliminationTemplate(x, J)
  E2 = eliminationTemplate(x, J)
  E3 = eliminationTemplate(x, J)
  
  -- Test 1: Default Strategy
  M1 = getActionMatrix(E1);
  evals1 = eigenvalues M1;
  assert(#evals1 == 12)
  
  -- Test 2: Larsson Strategy
  M2 = getActionMatrix(E2, Strategy => "Larsson");
  evals2 = eigenvalues M2
  assert(#evals2 == 12)
  
  -- Test 3: Greedy Strategy -- !! this is a good example, but 20s is probably too slow for a test
  M3 = getActionMatrix(E3, Strategy => "Greedy")
  evals3 = eigenvalues M3
  assert(#evals3 == 12)
///

TEST /// -- 5-point essential matrix problem
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;  -- essential matrix
  I = ideal(E*transpose E * E - (1/2) * trace(E * transpose E) * E);  -- Demazure constraints
  l = random(1, R);
  sols = templateSolve(l, I)
  assert(all(sols, x -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{x})))
///

TEST /// -- change of ideals
  R = QQ[x,y]
  I = ideal(x^2+y^2-1,x^2+y^3+x*y-2)
  E = eliminationTemplate(x+4*y,I)
  sols = templateSolve(E)
  assert(all(sols, x -> 1e-6 > norm sub(sub(gens I, QQ[gens R]), matrix{x})))

  J = ideal(x^2+y^2-2,x^2+y^3+3*x*y-5)
  F = copyTemplate(E,J)
  sols = templateSolve(F)
  assert(all(sols, x -> 1e-6 > norm sub(sub(gens J, QQ[gens R]), matrix{x})))
///

end


-- 5-point essential matrix problem: DEBUGGING TEMPLATE SIZE & STRATEGY
restart
path = prepend("./", path)
needsPackage "EliminationTemplates"
check "EliminationTemplates"
R = QQ[x,y,z]
Es = apply(4, i -> random(QQ^3, QQ^3))
E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3  -- essential matrix
I = ideal(E*transpose E * E - (1/2) * trace(E * transpose E) * E, det E);  -- Demazure constraints
--l = random(1, R)
l = y
ET = eliminationTemplate(l, I)
printWidth = 10000
getTemplateMatrix ET
load "Benchmarks.m2";
runBenchmarks()



-* Development section *-
-- basic solve, compare with known solution
restart
debug needsPackage "EliminationTemplates"
needsPackage "NumericalAlgebraicGeometry"
R=QQ[x,y,z]
J=ideal(x^3+y^3+z^3-4,x^2-y-z-1,x-y^2+z-3)
B=basis(R/J)
getEigenMatrix(x,J)
templateSolve(x, J)
templateSolve(x+2*y+3*z,J)
templateSolve(x,J)
netList solveSystem J_*

-- change of basis
restart
debug needsPackage "EliminationTemplates"
R=QQ[x,y]
I=ideal(x^2+y^2-1,x^2+y^3+x*y-2)
J=ideal(x^2+y^2-2,x^2+y^3+3*x*y-5)
B=basis(R/I)
E=eliminationTemplate(x+4*y,I)
getTemplate(E)
getEigenMatrix(E)
sols = templateSolve(E)
assert(all(sols, x -> 1e-6 > norm sub(sub(gens I, QQ[gens R]), matrix{x})))

F=copyTemplate(E,J)
getEigenMatrix(F)
sols = templateSolve(F)
assert(all(sols, x -> 1e-6 > norm sub(sub(gens J, QQ[gens R]), matrix{x})))

restart
debug needsPackage "EliminationTemplates"
R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
templateSolve(x, J)
actVar = x
getEigenMatrix(x, J)

restart
debug needsPackage "EliminationTemplates"
R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
templateSolve(x, J)

restart
debug needsPackage "EliminationTemplates"
R = QQ[x]
J = ideal(x^2-1, x^3-x)
getH0(x, basis(R/J), J, Strategy => "Greedy")

-- Benchmark tests: just run these three lines
restart
load "Benchmarks.m2";
runBenchmarks()
--

uninstallPackage "EliminationTemplates"
restart
installPackage "EliminationTemplates"
viewHelp "EliminationTemplates"
check "EliminationTemplates"

help EliminationTemplates
help getTemplate

viewHelp "EliminationTemplates"

-- 5-point essential matrix problem: DEBUGGING TEMPLATE SIZE & STRATEGY
restart
path = prepend("./", path)
needsPackage "EliminationTemplates"
R = QQ[x,y,z]
Es = apply(4, i -> random(QQ^3, QQ^3))
E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3  -- essential matrix
I = ideal(E*transpose E * E - (1/2) * trace(E * transpose E) * E, det E);  -- Demazure constraints
--l = random(1, R)
(sh, mp) = getTemplate ET
l = y
ET = eliminationTemplate(l, I)
getTemplateMatrix(ET); -- 27 X 44
getTemplateMatrix(ET, Strategy => "Greedy"); -- 15 x 44
getTemplateMatrix(ET, Strategy => "Larsson"); -- 24 x 44

-* 
-- problem! should be 24 x 34
Rosie's proposed solution: 
  1. Store most recently used strategy in cache of ET
  2. If NEW strategy is passed, recompute
*-


-- E+f+k 7pt relative pose
restart
needsPackage "EliminationTemplates"
R = QQ[w,x,y,lambda];
mons = {x^2, y^2, lambda^2, x*y, x*lambda, y*lambda};
coeffs = apply(6, i -> random(QQ));
h = sum(0..#mons-1, i -> coeffs#i * mons#i);  -- random quadratic function
Fs = apply(4, i -> random(QQ^3, QQ^3));
F = x * Fs#0 + y * Fs#1 + lambda * Fs#2 + Fs#3;
Q = diagonalMatrix({1, 1, w});
I = ideal(F * Q * transpose F * Q * F - (1/2) * trace(F * Q * transpose F * Q) * F) + ideal(det F) + ideal(lambda * y - h);
l = random(1, R)
errorDepth = 0 
ET = eliminationTemplate(l, I)
getTemplateMatrix(ET); -- 788 x 530
getTemplateMatrix(ET, Strategy => "Larsson"); -- 256 x 339
-- getTemplateMatrix(ET, Strategy => "Larsson"); -- will exceed runtime limit
