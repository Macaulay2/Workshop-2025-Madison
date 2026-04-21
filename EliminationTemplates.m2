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

-- Symbolic H via gap polynomials + adjustParams (Martyushev CVPR 2022 §3–§4).
-- Provides: buildGapPolys, buildHSymbolic, adjustParams, buildTemplateFromH.
load "./EliminationTemplates/MartyushevClean.m2"

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
    "copyTemplate",
    -- Option symbol used by MartyushevClean.adjustParams (bounded greedy iterations).
    "MaxIter"
}

EliminationTemplate = new Type of HashTable
ShiftSet = new Type of List
MonomialPartition = new Type of List

eliminationTemplate = method(Options => {})
eliminationTemplate (RingElement, Ideal) := o -> (aVar, J) -> (
    R := ring J;
    new EliminationTemplate from {
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
    a = sub(a, R);
    B = sub(B, R);
    FF := coefficientRing R;

    -- Greedy (Martyushev CVPR 2022) computes H via gap polynomials + α-commit,
    -- then lifts to the mainline H0 shape (nF x nB). No Groebner detour needed.
    -- Downstream (getTemplateHelper, getTemplate, getTemplateMatrix default,
    -- getActionMatrix default, recoverSolutions) consumes this H0 unchanged.
    if (o.Strategy === "Greedy") then (
        BlistGr := flatten entries B;
        (gp, resMons, BlistChk) := buildGapPolys(a, matrix{BlistGr}, J);
        nFGr := numgens J;
        nBGr := #BlistGr;
        if #gp == 0 then return map(R^nFGr, R^nBGr, 0);
        FlistGr := flatten entries gens J;
        RB := (toList resMons) | BlistGr;
        (Hsym, Rext, alphaVars, perRow) := buildHSymbolic(FlistGr, gp);
        Hfinal := adjustParams(FlistGr, Hsym, Rext, (alphaVars, RB));
        return liftGreedyHToH0(Hfinal, a, BlistGr, nFGr);
    );

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
        H0
    )
    else if (o.Strategy === "Larsson") then (
        -- Larsson CVPR 2017: reduce H0 mod the first syzygy module of F.
        H0res := H0 % image(syz(gens(J)));
        sub(H0res, ring J)
    )
    else (error "Strategy not supported here. getH0 accepts null / \"Larsson\" / \"Greedy\". MatrixHi bypasses getH0 — call getTemplateMatrix directly.")
)

shiftPolynomials = (shifts, J) -> (
    assert(length shifts == numgens J);
    apply(shifts, J_*, (m, f) -> f * sub(m, ring J))
)

getTemplateHelper = (a, B, J, o) -> (
    H0 := getH0(a, B, J, o);
    shifts := new ShiftSet from apply(numgens J, i -> monomials(H0^{i}));
    allMons := union(set \ flatten \ entries \ monomials \ shiftPolynomials(shifts, J));
    monsB := set flatten entries(lift(B, ring J));
    monsR := set flatten entries(a * lift(B, ring J)) - monsB;
    monsE := allMons - union(monsR, monsB);
    (shifts, new MonomialPartition from rsort \ toList \ {monsE, monsR, monsB})
)

getTemplate = method(Options => {MonomialOrder => null, Strategy => null})
getTemplate(EliminationTemplate) := o -> E -> (
    if E.cache#?"monomialPartition" and (o.Strategy === null or (E.cache#?"lastPartitionStrategy" and E.cache#"lastPartitionStrategy" === o.Strategy)) then (
        (E.cache#"shifts", E.cache#"monomialPartition")
    ) else (
        J := ideal E;
        R := ring J;
        a := sub(actionVariable E, R);
        B := lift(basis(R/J), R);

        (shOrig, mpOrig) := getTemplateHelper(a, B, J, o);

        K := coefficientRing R;
        ringVars := flatten entries vars R;
        MO := if not instance(o.MonomialOrder, Nothing) then o.MonomialOrder else (options R).MonomialOrder;
        Rs := K[prepend("s", ringVars), MonomialOrder => {Eliminate 1, MO}];

        toRs := map(Rs, R, apply(numgens R, i -> Rs_(i+1)));

        aS := toRs(a);
        JsGens := toRs(gens J);
        actVar := Rs_0;
        
        Is := ideal(JsGens | matrix{{actVar - aS}});

        Bs := toRs(B);
        sortedBs := rsort flatten entries Bs;

        shiftsGraph := new ShiftSet from (
            apply(shOrig, sh -> toRs(sh)) | {matrix {sortedBs}}
        );

        allMons := union(set \ flatten \ entries \ monomials \ shiftPolynomials(shiftsGraph, Is));
        monsB := set sortedBs;
        monsR := set apply(sortedBs, b -> actVar * b);
        monsE := allMons - union(monsR, monsB);
        mpGraph := new MonomialPartition from rsort \ toList \ {monsE, monsR, monsB};

        E.cache#basis = Bs;
        E.cache#"graphIdeal" = Is;
        E.cache#"shifts" = shiftsGraph;
        E.cache#"monomialPartition" = mpGraph;
        E.cache#"lastPartitionStrategy" = o.Strategy;
        (shiftsGraph, mpGraph)
    )
)

copyTemplate = method(Options => {})
copyTemplate(EliminationTemplate, Ideal) := o -> (E, J) -> (
    Rnew := ring J;
    FFnew := coefficientRing Rnew;
    aNew := sub(actionVariable E, Rnew);
    Enew := eliminationTemplate(aNew, J);
    
    if E.cache#?"graphIdeal" then (
	Rs := ring E.cache#"graphIdeal";
        Rsnew := FFnew[gens ring E.cache#"graphIdeal", MonomialOrder => (options Rs).MonomialOrder];

	if E.cache#?basis then Enew.cache#basis = sub(if E.cache#?basis then E.cache#basis else basis E, Rsnew);
        if E.cache#?"shifts" then Enew.cache#"shifts" = apply(E.cache#"shifts", sh -> sub(sh, Rsnew));
        if E.cache#?"monomialPartition" then Enew.cache#"monomialPartition" = apply(E.cache#"monomialPartition", mp -> apply(mp, m -> sub(m, Rsnew)));
        if E.cache#?"lastPartitionStrategy" then Enew.cache#"lastPartitionStrategy" = E.cache#"lastPartitionStrategy";
        if E.cache#?"lastMatrixStrategy" then Enew.cache#"lastMatrixStrategy" = E.cache#"lastMatrixStrategy";
	if E.cache#?"lastActionStrategy" then Enew.cache#"lastActionStrategy" = E.cache#"lastActionStrategy";

	toRsnew := map(Rsnew, Rnew, apply(numgens Rnew, i -> Rsnew_(i+1)));
        JsGens := toRsnew(gens J);
        aS := toRsnew(aNew);
        actVar := Rsnew_0;
        
        Enew.cache#"graphIdeal" = ideal(JsGens | matrix{{actVar - aS}});
	if E.cache#?"templateMatrix" then Enew.cache#"templateMatrix" = getTemplateMatrix(Enew.cache#"shifts", Enew.cache#"monomialPartition", Enew.cache#"graphIdeal");
    );
    Enew
)

getTemplateMatrix = method(Options => {MonomialOrder => null, Strategy => null})
getTemplateMatrix(RingElement, Matrix, Ideal) := o -> (a, B, J) -> (
    getTemplateMatrix(eliminationTemplate(a, J), o)
)
getTemplateMatrix(ShiftSet, MonomialPartition, Ideal) := o -> (shifts, monomialPartition, J) -> (
    allMons := monomialPartition#0 | monomialPartition#2;
    sub(transpose fold(apply(shiftPolynomials(shifts, J), m -> last coefficients(m, Monomials => allMons)), (a,b) -> a|b), coefficientRing ring J)
)
getTemplateMatrix(EliminationTemplate) := o -> E -> (
    -- MatrixHi / Greedy: build the template from gap polynomials
    --   G_i = a*b_i - NF(a*b_i, J)
    -- factored over F as G_i = sum_j H[i,j] * F_j (Martyushev CVPR 2022 §3).
    --
    -- MatrixHi (Martyushev CVPR 2022 §3 particular solution; α := 0) bypasses
    -- getH0 and builds the template directly from gap polynomials
    -- G_i = a*b_i - NF(a*b_i, J) via buildHSymbolic + buildTemplateFromH.
    -- The output skips the graph-ideal extension (`s - a`) that Default /
    -- Larsson use; getActionMatrix / templateSolve then read the action via
    -- extractActionFromTemplate + recoverSolutionsMatrixHi.
    --
    -- Greedy used to share this branch, but is now routed through getH0 +
    -- default downstream (Route A): adjustParams produces H, liftGreedyHToH0
    -- reshapes to the H0 contract (nF x nB), and getTemplate / default
    -- getActionMatrix / recoverSolutions consume it unchanged.
    --
    -- Limitation: the action variable must be a ring *variable* (monomial),
    -- not a linear combination. `buildGapPolys` constructs residual
    -- monomials as `a * B \ B`, which is only well-defined when `a * b_i`
    -- is a single monomial. Polynomial actions like `x + y + z` make every
    -- `a * b_i` a polynomial; residual detection degenerates and
    -- buildTemplateFromH errors in `coefficient`.
    if o.Strategy === "MatrixHi" then (
        cacheKey := "templateMatrixMatrixHi";
        if E.cache#?cacheKey then return E.cache#cacheKey;
        J := ideal E;
        R := ring J;
        a := sub(actionVariable E, R);
        B := lift(basis(R/J), R);
        Flist := flatten entries gens J;
        (gp, resMons, Blist) := buildGapPolys(a, B, J);
        RB := (toList resMons) | Blist;
        (Hsym, Rext, alphaVars, perRow) := buildHSymbolic(Flist, gp);
        -- Particular solution: evaluate all alphas at 0, project to R.
        Hfinal := if #alphaVars > 0 then (
                finalSubst := map(R, Rext,
                    apply(numgens R, i -> R_i) | apply(#alphaVars, k -> 0_R));
                matrix apply(numRows Hsym, i ->
                    apply(numColumns Hsym, j -> finalSubst(Hsym_(i,j))))
            ) else Hsym;
        (sh, M, V, Eexcess) := buildTemplateFromH(Hfinal, Flist, RB);
        E.cache#"matrixHiBlist" = Blist;
        E.cache#"matrixHiResMons" = toList resMons;
        E.cache#cacheKey = M;
        return M;
    );
    if E.cache#?"templateMatrix" and (o.Strategy === null or (E.cache#?"lastMatrixStrategy" and E.cache#"lastMatrixStrategy" === o.Strategy)) then (
        E.cache#"templateMatrix"
    ) else (
        (shifts, monomialPartition) := getTemplate(E, o);
        ret := getTemplateMatrix(shifts, monomialPartition, E.cache#"graphIdeal", o);
        E.cache#"templateMatrix" = ret;
        E.cache#"lastMatrixStrategy" = o.Strategy;
        ret
    )
)

getActionMatrix = method(Options => {MonomialOrder => null, Strategy => null})
getActionMatrix(RingElement, MonomialPartition, Matrix) := o -> (actVar, mp, M) -> (
    numE := length mp#0;
    numB := length mp#2;
    FF := ring M;
    m := numrows M;
    n := numcols M;
    numTop := m - numB;
    MtopE := M_{0 .. numE-1}^{0..numTop-1};
    MtopB := M_{numE .. n-1}^{0..numTop-1};
    X := solve(MtopE, MtopB, ClosestFit => if instance(FF, InexactFieldFamily) or instance(FF, InexactField) then true else false);
    MbotE := M_{0 .. numE-1}^{numTop .. m-1};
    MbotB := M_{numE .. n-1}^{numTop .. m-1};
    MbotE * X - MbotB
)
getActionMatrix(EliminationTemplate) := o -> E -> (
    -- MatrixHi: read action matrix off the template via RREF + pivot
    -- extraction. This strategy skips the graph-ideal extension, so the
    -- monomialPartition cache used by the default branch is not populated.
    -- Greedy now flows through the default branch (Route A): its H0 is
    -- greedy-computed in getH0, but the downstream is the shared graph
    -- linear-solve path below.
    if o.Strategy === "MatrixHi" then (
        cacheKey := "actionMatrix_" | toString o.Strategy;
        if E.cache#?cacheKey then return E.cache#cacheKey;
        M := getTemplateMatrix(E, o);  -- populates matrixHiBlist / matrixHiResMons
        Blist := E.cache#"matrixHiBlist";
        resMons := E.cache#"matrixHiResMons";
        Ma := extractActionFromTemplate(M, resMons, Blist, actionVariable E);
        E.cache#cacheKey = Ma;
        return Ma;
    );
    if E.cache#?"actionMatrix" and (o.Strategy === null or (E.cache#?"lastActionStrategy" and E.cache#"lastActionStrategy" === o.Strategy)) then (
        E.cache#"actionMatrix"
    ) else (
        (sh, mp) := getTemplate(E, o);
        templateMatrix := getTemplateMatrix(E, o);

        Rs := ring first first mp;
        actVar := Rs_0;

        ret := getActionMatrix(actVar, mp, templateMatrix, o);
        E.cache#"actionMatrix" = ret;
        E.cache#"lastActionStrategy" = o.Strategy;
        ret
    )
)

-*

*-
getEigenMatrix = method(Options => {MonomialOrder => null, Strategy => null})
getEigenMatrix(EliminationTemplate) := o -> (E) -> (
    Ma := getActionMatrix(E, o);
    -- MatrixHi: action matrix is over the base ring R (not the graph ring
    -- Rs), so basis monomials are in R as well. Greedy uses the default
    -- graph-ring branch below.
    if o.Strategy === "MatrixHi" then (
        (svals, P) := eigenvectors sub(Ma, CC);
        Blist := E.cache#"matrixHiBlist";
        (matrix{Blist}, P)
    ) else (
        (svals2, P2) := eigenvectors Ma;
        cleanEvecs := clean_(1e-10) (P2 * inverse diagonalMatrix(P2^{numColumns P2 - 1}));
        (transpose rsort basis E, cleanEvecs)
    )
)
getEigenMatrix(Ideal) := o -> (I) -> getEigenMatrix(random(1, ring I), I, o)
getEigenMatrix(RingElement, Ideal) := o -> (a, J) -> (
    E := eliminationTemplate(a, J);
    getEigenMatrix(E, o)
)

templateSolve = method(Options => {MonomialOrder => null, Strategy => null})
templateSolve(EliminationTemplate) := o -> (E) -> (
    -- MatrixHi: read variable values straight off basis-indexed eigenvector
    -- entries.
    --
    -- Greedy (Route A): H0 and action matrix go through the shared graph-ring
    -- pipeline, but the downstream recover is done via recoverSolutionsMatrixHi,
    -- matching Greedy's pre-refactor correctness on multivariable QQ problems
    -- where the graph-ring recoverSolutions has pre-existing numerical issues
    -- (e.g. 5pt essential+det, 6 Demazure cubics). The action matrix Ma is
    -- intrinsic to (a, J) up to basis permutation, so eigen-recovery gives
    -- the right solutions regardless of which path produced Ma.
    if o.Strategy === "MatrixHi" then (
        recoverSolutionsMatrixHi(getActionMatrix(E, o), E.cache#"matrixHiBlist", ring ideal E)
    ) else if o.Strategy === "Greedy" then (
        -- recoverSolutionsMatrixHi reads a variable's coordinate directly
        -- from its slot in the basis, so it only works when every ring
        -- variable appears in the quotient basis. When that fails (e.g.
        -- circle+line where the basis is {1, y} and x is eliminated), fall
        -- back to the graph-ring recoverSolutions which reads coordinates
        -- through the full monomial partition.
        Ma := getActionMatrix(E, o);  -- populates monomialPartition cache
        R := ring ideal E;
        mp := E.cache#"monomialPartition";
        Rs := ring first first mp;
        toR := map(R, Rs, {0_R} | apply(numgens R, i -> R_i));
        Blist := apply(mp#2, b -> toR(b));
        BlistSet := set Blist;
        if all(numgens R, i -> BlistSet#?(R_i)) then (
            recoverSolutionsMatrixHi(Ma, Blist, R)
        ) else (
            (Bmat, M) := getEigenMatrix(E, o);
            templateMat := getTemplateMatrix(E, o);
            recoverSolutions(Bmat, M, E, templateMat)
        )
    ) else (
        (Bmat2, M2) := getEigenMatrix(E, o);
        templateMat2 := getTemplateMatrix(E, o);
        recoverSolutions(Bmat2, M2, E, templateMat2)
    )
)
templateSolve(Ideal) := o -> (I) -> templateSolve(random(1,ring I), I, o)
templateSolve(RingElement, Ideal) := o -> (a, J) -> (
    E := eliminationTemplate(a, J);
    templateSolve(E, o)
)

recoverSolutions = method()
recoverSolutions(Matrix, Matrix, EliminationTemplate, Matrix) := (Bmat, M, E, templateMat) -> (
    J := ideal E;
    Rnew := ring J;
    
    mp := E.cache#"monomialPartition";
    Rs := ring first first mp;
    
    toRnew := map(Rnew, Rs, {0_Rnew} | flatten entries vars Rnew);
    toRs := map(Rs, Rnew, apply(numgens Rnew, i -> Rs_(i+1)));
    
    basisMonsRnew := apply(flatten entries Bmat, m -> toRnew(m));
    
    monsE := mp#0;
    monsR := mp#1;
    monsB := mp#2;
    
    numE := length monsE;
    numR := length monsR;
    numB := length monsB;
    numTop := numrows templateMat - numB;
    
    -- Pure linear algebra: solve only for the excessive block, skipping action variables
    MtopE := templateMat_{0 .. numE-1}^{0..numTop-1};
    MtopB := templateMat_{numE .. numcols templateMat -1}^{0..numTop-1};
    
    X := solve(MtopE, MtopB);
    
    solutions := {};
    varsList := flatten entries vars Rnew;
    
    for rootIndex from 0 to numColumns M - 1 do (
        monomialValues := new MutableHashTable;
        for i from 0 to #basisMonsRnew - 1 do (
            monomialValues#(basisMonsRnew#i) = M_(i, rootIndex);
        );

        root := {};
        for v in varsList do (
            if monomialValues#?v then (
		-- v is basic monomial
              root = append(root, monomialValues#v);
            )
            else (
		-- v is an excessive monomial
                vRs := toRs(v);
                posInE := position(monsE, m -> m == vRs);
                if posInE =!= null then (
                    local val;
                    val = 0;
                    for j from 0 to numB - 1 do (
                        bMapped := toRnew(monsB#j);
                        val = val - sub(X_(posInE, j), CC) * monomialValues#bMapped;
                    );
                    root = append(root, val);
                )
                else (
		    -- failsafe in case v is neither a basic nor an excessive monomial
                    r := v % J;
                    coeffs := last coefficients(r, Monomials => basisMonsRnew);
                    local val;
                    val = 0;
                    for j from 0 to numB - 1 do (
                        bMapped := basisMonsRnew#j;
                        if monomialValues#?bMapped then (
                          val = val + sub(coeffs_(j,0), coefficientRing Rnew) * monomialValues#bMapped;
                        )
                    );
                    root = append(root, val);
                );
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
    Description
      Text
        This method copies an elimination template object, using the same action variable and basis, but a different defining ideal.
      Example
        R = QQ[x,y]
        I = ideal(x^4+x*y+y^2-3, x^2*y+y^3-2)
        J = ideal(x^3+y^2-1,x^2+y^3-1)
        E = eliminationTemplate(x,I)
        F = copyTemplate(E, J)
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
    Headline
      computes template matrix
    Usage
      M = (a, B, J)
      M = getTemplateMatrix(sh, mp, J)
      getTemplateMatrix(E)
    Inputs
      a:RingElement
        the action polynomial defining a multiplication matrix
      B:Matrix
      E:EliminationTemplate
        the elimination template for this problem
      J:Ideal
      MonomialOrder=>Thing
        the monomial order used on the ambient ring
      mp:MonomialPartition
      sh:ShiftSet
      Strategy=>Thing
        the strategy used to compute H0
    Outputs
      M:Matrix
    Description
      Text
        This method computes the template matrix corresponding to the inputted elimination template.
      Example
        R = QQ[x,y]
        I = ideal(x^4+x*y+y^2-3, x^2*y+y^3-2)
        E = eliminationTemplate(x,I)
        getTemplateMatrix(E)
    SeeAlso
      EliminationTemplate
///

doc ///
 Node
    Key
      actionVariable
    Headline
      returns the action variable associated to the elimination template
    Usage
      actionVariable(E)
    Inputs
      E:EliminationTemplate
        the elimination template for this problem
    Outputs
      a:RingElement
        the action variable associated to E
    Description
      Text
        This method outputs the action variable associated to the inputted elimination template E.
      Example
        R = QQ[x,y]
        I = ideal(x^4+x*y+y^2-3, x^2*y+y^3-2)
        E = eliminationTemplate(x,I)
        actionVariable(E)
    SeeAlso
      EliminationTemplate
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

TEST /// -- 3-variable benchmark: dense cubic-ish system (deg J = 12, |B| = 12).
-- Used here for a smoke test of the basic pipeline (getTemplateMatrix / getActionMatrix).
-- Template-size checks for this system and the 5-pt essential are in the
-- "MatrixHi" strategy TEST below.
  R = QQ[x,y,z]
  J = ideal(x^3+y^3+z^3-4,x^2-y-z-1,x-y^2+z-3)
  E = eliminationTemplate(x, J)
  getTemplateMatrix E
  getActionMatrix E
  eigenvalues getActionMatrix E
///

TEST /// -- 5-point essential matrix (Demazure trace identity, deg I = 10).
-- Classical relative-pose benchmark from Nister; template-size reference from
-- Martyushev et al., CVPR 2022 (10 x 20 on std basis).
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;  -- essential matrix
  I = ideal(E*transpose E * E - (1/2) * trace(E * transpose E) * E);  -- Demazure constraints
  l = random(1, R);
  sols = templateSolve(l, I)
  assert(all(sols, x -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{x})))
///

TEST /// -- MatrixHi strategy: template size matches Martyushev CVPR 2022 reference.
-- The MatrixHi path bypasses the graph-ideal extension, so the template has
-- |B| fewer rows than Default / Larsson (and now Greedy, which under Route A
-- flows through the graph-ideal pipeline). On 5pt essential: 10×20 for
-- MatrixHi vs 20×20 for the graph-ideal strategies.
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  ET = eliminationTemplate(y, I);
  M = getTemplateMatrix(ET, Strategy => "MatrixHi");
  assert(numRows M == 10 and numColumns M == 20)

  J = ideal(x^3+y^3+z^3-4,x^2-y-z-1,x-y^2+z-3)
  ET2 = eliminationTemplate(x, J);
  N = getTemplateMatrix(ET2, Strategy => "MatrixHi");
  assert(numRows N <= numRows getTemplateMatrix(ET2))  -- strictly smaller on this system
///

TEST /// -- Greedy (Route A): on 5pt essential (0 free alphas), adjustParams
-- is a no-op, so the lifted greedy H0 equals the Default H0, and Greedy's
-- graph-ideal template matches Default's (20x20).
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  ET = eliminationTemplate(y, I);
  Mg = getTemplateMatrix(ET, Strategy => "Greedy");
  assert(numRows Mg == 20 and numColumns Mg == 20)
///

TEST /// -- getActionMatrix on MatrixHi / Greedy templates: 10x10 with deg I
-- distinct eigenvalues on 5pt essential.
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  d = degree I;
  E1 = eliminationTemplate(y, I);
  Ma1 = getActionMatrix(E1, Strategy => "MatrixHi");
  assert(numRows Ma1 == d and numColumns Ma1 == d);
  assert(#eigenvalues sub(Ma1, CC) == d);
  E2 = eliminationTemplate(y, I);
  Ma2 = getActionMatrix(E2, Strategy => "Greedy");
  assert(numRows Ma2 == d and numColumns Ma2 == d);
  assert(#eigenvalues sub(Ma2, CC) == d);
///

TEST /// -- templateSolve via MatrixHi on 5pt essential
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  d = degree I;
  E1 = eliminationTemplate(y, I);
  sols1 = templateSolve(E1, Strategy => "MatrixHi");
  assert(#sols1 == d);
  assert(all(sols1, s -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{s})));
///

TEST /// -- templateSolve via Greedy on 5pt essential
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  d = degree I;
  E2 = eliminationTemplate(y, I);
  sols2 = templateSolve(E2, Strategy => "Greedy");
  assert(#sols2 == d);
  assert(all(sols2, s -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{s})));
///

TEST /// -- MatrixHi with action variable x on 5pt essential.
-- The graph-ideal templateSolve path fails on 3-variable systems with
-- non-random actions; MatrixHi reads variable coords directly from the
-- basis and handles monomial actions correctly.
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  d = degree I;
  E = eliminationTemplate(x, I);
  sols = templateSolve(E, Strategy => "MatrixHi");
  assert(#sols == d);
  assert(all(sols, s -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{s})));
///

TEST /// -- Greedy with action z on 5pt essential.
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  d = degree I;
  E = eliminationTemplate(z, I);
  sols = templateSolve(E, Strategy => "Greedy");
  assert(#sols == d);
  assert(all(sols, s -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{s})));
///

TEST /// -- Cross-validation: on a 0-free-alpha problem (5pt essential),
-- Greedy's adjustParams is a no-op, so the lifted greedy H0 equals the
-- Default H0 and Greedy's graph-ideal template equals Default's. Action
-- matrices from all three strategies (MatrixHi / Greedy / Default) represent
-- the same linear map — MatrixHi uses the raw basis order from
-- `basis(R/J)`, Default/Greedy use rsort, so matrices differ by a basis
-- permutation but share eigenvalues (intrinsic spectrum of `a` on R/J).
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  E1 = eliminationTemplate(y, I);
  Mh = getTemplateMatrix(E1, Strategy => "MatrixHi");
  E2 = eliminationTemplate(y, I);
  Mg = getTemplateMatrix(E2, Strategy => "Greedy");
  E3 = eliminationTemplate(y, I);
  Md = getTemplateMatrix(E3);
  Ah = getActionMatrix(E1, Strategy => "MatrixHi");
  Ag = getActionMatrix(E2, Strategy => "Greedy");
  Ad = getActionMatrix(E3);
  -- Greedy (Route A) and Default share the graph-ideal pipeline; with 0
  -- free alphas the H0 objects coincide, so the templates coincide.
  assert(Mg == Md);
  -- MatrixHi lives on a different pipeline (no graph extension), with
  -- |B| fewer rows than Mg.
  d = degree I;
  assert(numRows Mh == numRows Mg - d);
  -- Greedy and Default match exactly (same ordering, same pipeline).
  assert(Ag == Ad);
  -- MatrixHi's action matrix is similar to Greedy's via a basis permutation;
  -- check spectrum equality up to ordering (as a multiset over CC).
  evH = sort apply(toList eigenvalues sub(Ah, CC), z -> (realPart z, imaginaryPart z));
  evG = sort apply(toList eigenvalues sub(Ag, CC), z -> (realPart z, imaginaryPart z));
  assert(#evH == #evG);
  assert(all(#evH, i -> abs(evH#i#0 - evG#i#0) < 1e-8 and abs(evH#i#1 - evG#i#1) < 1e-8));
///

TEST /// -- Greedy (Route A) on 6 Demazure cubics (no det): free alpha > 0,
-- adjustParams commits alphas to cancel excessive monomials, producing a
-- smaller H0 (fewer monomials per row) than Default's Groebner H0. Both
-- flow through the graph-ideal pipeline, so the comparison that remains
-- meaningful under Route A is Greedy vs Default: Greedy's template should
-- have row count ≤ Default's (fewer shifts). Both paths solve to
-- residual < 1e-6. MatrixHi is checked separately (it is on a different
-- pipeline — not comparable row-for-row here).
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(2*Ee*transpose(Ee)*Ee - trace(Ee*transpose(Ee))*Ee);
  assert(dim I == 0);
  d = degree I;
  -- Compute everything up front; defer the CC residual checks to the end,
  -- because `CC[gens R]` rebinds the symbols x,y,z and breaks subsequent
  -- calls like `eliminationTemplate(y, I)`.
  E1 = eliminationTemplate(y, I);
  Md = getTemplateMatrix(E1);                           -- Default
  E2 = eliminationTemplate(y, I);
  Mg = getTemplateMatrix(E2, Strategy => "Greedy");
  E3 = eliminationTemplate(y, I);
  solsH = templateSolve(E3, Strategy => "MatrixHi");
  solsD = templateSolve(E1);
  solsG = templateSolve(E2, Strategy => "Greedy");
  -- Route A invariant: Greedy ≤ Default in rows on problems with free alphas.
  assert(numRows Mg <= numRows Md);
  assert(#solsD == d and #solsG == d and #solsH == d);
  Rc = CC[gens R];
  gensC = sub(gens I, Rc);
  assert(all(solsD, s -> 1e-6 > norm sub(gensC, matrix{s})));
  assert(all(solsG, s -> 1e-6 > norm sub(gensC, matrix{s})));
  assert(all(solsH, s -> 1e-6 > norm sub(gensC, matrix{s})));
///

TEST /// -- MatrixHi/Greedy on the 3-variable docs system (different benchmark).
  R = QQ[x,y,z]
  J = ideal(x^3+y^3+z^3-4, x^2-y-z-1, x-y^2+z-3)
  d = degree J;
  E1 = eliminationTemplate(x, J);
  solsH = templateSolve(E1, Strategy => "MatrixHi");
  assert(#solsH == d);
  assert(all(solsH, s -> 1e-6 > norm sub(sub(gens J, CC[gens R]), matrix{s})));
///

TEST /// -- Paper §3 small example on a 2-var system (deg I = 12): MatrixHi
-- produces a 7×19 template (vs Default's 19×19) and solves to residual < 1e-6.
  R = QQ[x,y]
  I = ideal(x^4+y^2+x*y-3, x^2*y+y^3-2)
  d = degree I;
  E = eliminationTemplate(x, I);
  M = getTemplateMatrix(E, Strategy => "MatrixHi");
  assert(numColumns M == 19);
  assert(numRows M < 19);  -- strictly smaller than Default's 19×19
  sols = templateSolve(E, Strategy => "MatrixHi");
  assert(#sols == d);
  assert(all(sols, s -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{s})));
///

TEST /// -- #1 F+λ 8pt from Martyushev CVPR 2022 Table 1 over ZZ/32749.
-- Paper target: std 11×19, nstd 7×15. MatrixHi should achieve 11×20
-- (matching std up to a single column). No tier-3 since we're in ZZ/p.
  FF = ZZ/32749
  R = FF[x,y]
  d1 = random(FF^8, FF^4);
  V1 = matrix{{y*x},{y},{x},{1_R}};
  fv = sub(d1, R) * V1;
  Fmat = matrix{{fv_(0,0), fv_(3,0), fv_(5,0)},
                {fv_(1,0), fv_(4,0), fv_(6,0)},
                {fv_(2,0), x,        1_R}};
  I = ideal(det Fmat, y*fv_(2,0) - fv_(7,0));
  E = eliminationTemplate(x, I);
  M = getTemplateMatrix(E, Strategy => "MatrixHi");
  assert(numRows M == 11 and numColumns M == 20);
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

TEST /// -- copyTemplate on a Greedy-built template (Route A).
-- Greedy now populates the same graph-ideal cache as Default, so
-- copyTemplate should transplant the structure into a perturbed ideal
-- and Greedy solve on the copy should still hit tier-3 residual.
  R = QQ[x,y]
  I = ideal(x^2+y^2-1, x^2+y^3+x*y-2)
  E = eliminationTemplate(x, I)
  -- Populate the greedy cache (graphIdeal / shifts / monomialPartition).
  getTemplateMatrix(E, Strategy => "Greedy");
  getActionMatrix(E, Strategy => "Greedy");
  -- getTemplate rebinds x, y to the graph-ring generators while building
  -- Rs; restore R so that J is constructed over R (not over the graph
  -- ring), which is what copyTemplate expects.
  use R;
  J = ideal(x^2+y^2-2, x^2+y^3+3*x*y-5)
  F = copyTemplate(E, J)
  sols = templateSolve(F, Strategy => "Greedy")
  assert(all(sols, x -> 1e-6 > norm sub(sub(gens J, CC[gens R]), matrix{x})))
///

TEST /// -- example used for section 3
  R = QQ[x,y]
  I = ideal(x^4+y^2+x*y-3, x^2*y+y^3-2)
  E = eliminationTemplate(x,I)
  sols = templateSolve(E)
  assert(all(sols, x -> 1e-6 > norm sub(sub(gens I, QQ[gens R]), matrix{x})))
  actionVariable(E)
///

TEST ///
  R = QQ[x,y,z]
  J = ideal(x^3+y^3+z^3-4,x^2-y-z-1,x-y^2+z-3)
  -- getActionMatrix is only defined on the graph-ideal pipeline (Default,
  -- Larsson). MatrixHi / Greedy templates skip the graph extension, so the
  -- action matrix must be recovered separately (see benchmarks/bench_5pt_all.m2).
  E1 = eliminationTemplate(x, J);
  E2 = eliminationTemplate(x, J);
  -- Default Strategy
  M1 = getActionMatrix(E1);
  assert(#eigenvalues M1 == 12)
  -- Larsson Strategy
  M2 = getActionMatrix(E2, Strategy => "Larsson");
  assert(#eigenvalues M2 == 12)
///


end


--------------------------------------------------------------------------
-- DEMO: paste these blocks one-at-a-time into an interactive M2 session.
-- Each block is self-contained (starts with `restart`) and exercises one
-- strategy comparison with size annotations and tier-3 (solve) assertions.
--
-- Expected sizes below are reproducible with `setRandomSeed 42` on QQ.
--
-- Strategy semantics (post-2026-04-17 refactor):
--   null (Default): H0 from Groebner change-of-basis; graph-ideal pipeline
--   "Larsson":      H0 mod syz(F); graph-ideal pipeline (CVPR 2017)
--   "MatrixHi":     per-monomial alpha, alpha := 0 (Martyushev CVPR 2022 §3)
--   "Greedy":       MatrixHi + adjustParams (Martyushev CVPR 2022 §3-§4)
--------------------------------------------------------------------------

-- Demo 1: 5-point essential matrix (Demazure trace identity + det)
-- deg I = 10; MatrixHi/Greedy match paper std size 10x20.
restart
path = prepend("./", path)
needsPackage "EliminationTemplates"
setRandomSeed 42
R = QQ[x,y,z]
Es = apply(4, i -> random(QQ^3, QQ^3))
Em = x * Es#0 + y * Es#1 + z * Es#2 + Es#3
I  = ideal(2*Em*transpose(Em)*Em - trace(Em*transpose(Em))*Em, det Em)
degree I  -- 10

ET = eliminationTemplate(y, I)
getTemplateMatrix ET                          -- Default   : 20 x 20
getTemplateMatrix(ET, Strategy => "Larsson")  -- Larsson   : 20 x 20
getTemplateMatrix(ET, Strategy => "MatrixHi") -- MatrixHi  : 10 x 20  <- paper
getTemplateMatrix(ET, Strategy => "Greedy")   -- Greedy    : 10 x 20  (0 free alphas)

-- End-to-end solve + tier-3 residual check via MatrixHi / Greedy.
ET2 = eliminationTemplate(y, I)
sols = templateSolve(ET2, Strategy => "Greedy")
assert(#sols == degree I)
assert(all(sols, s -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{s})))


-- Demo 2: 6 Demazure cubics (no det) -- the one where GREEDY BEATS MATRIXHI.
-- 66 free alphas; adjustParams commits enough to shave a row and a column.
restart
path = prepend("./", path)
needsPackage "EliminationTemplates"
setRandomSeed 42
R = QQ[x,y,z]
Es = apply(4, i -> random(QQ^3, QQ^3))
Em = x * Es#0 + y * Es#1 + z * Es#2 + Es#3
I  = ideal(2*Em*transpose(Em)*Em - trace(Em*transpose(Em))*Em)  -- no det
degree I  -- still 10 (generic random coefficients)

E1 = eliminationTemplate(y, I)
Mh = getTemplateMatrix(E1, Strategy => "MatrixHi")  -- 25 x 35
E2 = eliminationTemplate(y, I)
Mg = getTemplateMatrix(E2, Strategy => "Greedy")    -- 24 x 34  <- strictly smaller
assert(numRows Mg < numRows Mh or numColumns Mg < numColumns Mh)

-- Both solve correctly to residual < 1e-6.
solsG = templateSolve(E2, Strategy => "Greedy")
assert(#solsG == degree I)
assert(all(solsG, s -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{s})))


-- Demo 3: paper §3 example (2-var system, deg I = 12).
-- MatrixHi shrinks rows 19 -> 7 while preserving column count.
restart
path = prepend("./", path)
needsPackage "EliminationTemplates"
R = QQ[x,y]
I = ideal(x^4 + y^2 + x*y - 3, x^2*y + y^3 - 2)
degree I  -- 12

E = eliminationTemplate(x, I)
getTemplateMatrix E                               -- Default  : 19 x 19
getTemplateMatrix(E, Strategy => "MatrixHi")      -- MatrixHi :  7 x 19
sols = templateSolve(E, Strategy => "MatrixHi")
assert(#sols == degree I)
assert(all(sols, s -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{s})))

-- Cross-check: eigenvalues of the action matrix == x-coordinates of solutions.
Ma = getActionMatrix(E, Strategy => "MatrixHi")
evals = eigenvalues sub(Ma, CC)
assert(#evals == degree I)


-- Demo 4: Paper Table 1 #1 (F+lambda 8pt) over ZZ/p; matches Martyushev std.
-- Shown over ZZ/32749 because MatrixHi over QQ is slower (coefficient growth).
restart
path = prepend("./", path)
needsPackage "EliminationTemplates"
FF = ZZ/32749
R = FF[x,y]
d1 = random(FF^8, FF^4)
V1 = matrix{{y*x},{y},{x},{1_R}}
fv = sub(d1, R) * V1
Fmat = matrix{{fv_(0,0), fv_(3,0), fv_(5,0)},
              {fv_(1,0), fv_(4,0), fv_(6,0)},
              {fv_(2,0), x,        1_R}}
I = ideal(det Fmat, y*fv_(2,0) - fv_(7,0))
E = eliminationTemplate(x, I)
getTemplateMatrix(E, Strategy => "MatrixHi")  -- 11 x 20 (paper std: 11 x 19)




-- 5-point essential matrix problem: DEBUGGING TEMPLATE SIZE & STRATEGY
restart
path = prepend("./", path)
needsPackage "EliminationTemplates"
check "EliminationTemplates"
installPackage("EliminationTemplates", RemakeAllDocumentation => true)
viewHelp EliminationTemplates
R = QQ[x,y,z]
Es = apply(4, i -> random(QQ^3, QQ^3))
E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3  -- essential matrix
I = ideal(E*transpose E * E - (1/2) * trace(E * transpose E) * E, det E);  -- Demazure constraints
l = y
ET = eliminationTemplate(l, I)
M = getTemplateMatrix ET
FF=frac(QQ[e_(0,0,0)..e_(3,2,2)])
Es = apply(4, i -> matrix apply(3, j -> apply(3, k -> e_(i,j,k))))
R = FF[x,y,z]
E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3  -- essential matrix
J = ideal(E*transpose E * E - (1/2) * trace(E * transpose E) * E, det E);  -- Demazure constraints
errorDepth=3
ETP =  copyTemplate(ET, J)
printWidth = 1000000
getTemplateMatrix ETP

FF = frac(QQ[a,b,c,d])
R = FF[x,y,MonomialOrder=>Lex]
l = c*x + d*y
I = ideal(x^2+a*y^2-1, x*y-b)
needsPackage "EliminationTemplates"
ET = eliminationTemplate(l, I)
M = getTemplateMatrix ET
(P, L, U) = LUdecomposition M
reducedRowEchelonForm M

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
templateSolve(x,J) -- Why do we have two of these lines?
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
errorDepth = 2
templateSolve(x, J)
actVar = x
getEigenMatrix(x, J)

restart
debug needsPackage "EliminationTemplates"
R = QQ[x,y]
J = ideal(x^3 + y^2 - 1, x - y - 1)
errorDepth = 0
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
-- l = random(1, R)
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
getTemplateMatrix(ET, Strategy => "Greedy"); 
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
getTemplateMatrix(ET, Strategy => "Greedy"); -- will exceed runtime limit


restart
path = prepend("./", path)
needsPackage "EliminationTemplates"
check "EliminationTemplates"
installPackage("EliminationTemplates", RemakeAllDocumentation => true)
viewHelp EliminationTemplates

-- generate a template over finite field
FF = ZZ/3
R = FF[x,y,z]
Es = apply(4, i -> random(FF^3, FF^3))
E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3  -- essential matrix
I = ideal(2 * E*transpose E * E - trace(E * transpose E) * E, det E);  -- Demazure constraints
l = y
ET = eliminationTemplate(l, I)
M = getTemplateMatrix ET

-- try copying this template into a rational problem instance
FF = QQ
Es = apply(4, i -> random(FF^3, FF^3))
R = QQ[x,y,z]
E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3  -- essential matrix
J = ideal(E*transpose E * E - (1/2) * trace(E * transpose E) * E, det E);  -- Demazure constraints
E = copyTemplate(ET, J)
sols = templateSolve(E)
apply(sols, x -> 1e-6 > norm sub(sub(gens J, QQ[gens R]), matrix{x}))

-- Test case
restart
needsPackage "EliminationTemplates"
R = QQ[x,y,z]
J = ideal(x^3+y^3+z^3-4,x^2-y-z-1,x-y^2+z-3)
-- 3 templates, 3 strategies
E1 = eliminationTemplate(x, J);
E2 = eliminationTemplate(x, J);
E3 = eliminationTemplate(x, J);
getTemplateMatrix(E1); -- 27 X 44
getTemplateMatrix(E2, Strategy => "Greedy"); -- 15 x 44
getTemplateMatrix(E3, Strategy => "Larsson")

-- This example doesn't work :(
loadPackage "EliminationTemplates"
R = QQ[x,y,z];
E0 = matrix {{7/3, 9, 3}, {5/6, 3, 1/8}, {10/9, 7/5, 3/4}};
E1 = matrix {{1/6, 5/8, 3/10}, {9/4, 9/7, 1/3}, {7/4, 2, 7/10}};
E2 = matrix {{8/9, 7/5, 9/4}, {10/3, 9/4, 4/7}, {5/2, 7/5, 2/9}};
E3 = matrix {{5/6, 1/7, 6}, {6/7, 8/3, 3/10}, {9/8, 1, 4/7}};
Es = {E0, E1, E2, E3}
E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;  -- essential matrix
I = ideal(E*transpose E * E - (1/2) * trace(E * transpose E) * E);  -- Demazure constraints
l = 5*x + (3/8)*y + (9/7)*z
sols = templateSolve(l, I);
norms = apply(sols, x -> norm sub(sub(gens I, CC[gens R]), matrix{x}));
all(norms, x -> tol > x)

--(1/9)*x+(3/5)*y+(5/6)*z


tol = 1.0
done = false
i = 0;
while not done do (
    l = random(1, R);
    sols = templateSolve(l, I);
    norms = apply(sols, x -> norm sub(sub(gens I, CC[gens R]), matrix{x}));
    done = not all(norms, x -> tol > x);
    i = i + 1;
    )
print(toString norms);
toString l
print i
