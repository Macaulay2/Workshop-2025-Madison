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
    "MaxIter",
    -- Option symbol on pipeline methods: `AdjustParams => false` runs Greedy
    -- up to the buildHSymbolic particular solution (α := 0) and skips the
    -- adjustParams optimization step — the fast "MatrixHi mode" of the
    -- Greedy strategy. Default true.
    "AdjustParams",
    -- Option symbols used by the basis-search infrastructure in MartyushevClean.
    "RunGreedy",
    "RunGreedyAll"
}

EliminationTemplate = new Type of HashTable
ShiftSet = new Type of List
MonomialPartition = new Type of List

-- True when `a` is a monomial (single term with unit coefficient), so `a * b_i`
-- is a single monomial for every basis monomial `b_i`. Direction A uses this
-- to skip the graph-ring extension for Default / Larsson / Greedy, saving |B|
-- rows and matching the MatrixHi / paper template structure.
isMonomialAction = a -> a == leadMonomial a

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

getH0 = method(Options => {MonomialOrder => null, Strategy => null, AdjustParams => true})
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
    -- then lifts to the mainline H0 shape (nF x nB). For QQ input with
    -- AdjustParams => true, pin the whole computation (including adjustParams'
    -- RREF on α's) to ZZ/p to avoid rational coefficient blowup; the result
    -- is lifted back to R via per-coefficient ZZ→QQ. Shift monomial support
    -- — all the downstream (getTemplateHelper, getTemplate) consumes — is
    -- preserved by the lift, which is what matters for template construction.
    -- Downstream reads only monomials(H0^{i}) per row, ignoring coefficients.
    if (o.Strategy === "Greedy") then (
        pinToModP := (FF === QQ) and o.AdjustParams;
        Rw := if pinToModP then (
            (ZZ/32749)[gens R, MonomialOrder => (options R).MonomialOrder]
        ) else R;
        aw := if pinToModP then sub(a, Rw) else a;
        Jw := if pinToModP then sub(J, Rw) else J;
        Bw := if pinToModP then sub(B, Rw) else B;
        BlistGr := flatten entries Bw;
        (gp, resMons, BlistChk) := buildGapPolys(aw, matrix{BlistGr}, Jw);
        nFGr := numgens Jw;
        nBGr := #BlistGr;
        if #gp == 0 then return map(R^nFGr, R^nBGr, 0);
        FlistGr := flatten entries gens Jw;
        RB := (toList resMons) | BlistGr;
        (Hsym, Rext, alphaVars, perRow) := buildHSymbolic(FlistGr, gp);
        Hfinal := if o.AdjustParams then (
            adjustParams(FlistGr, Hsym, Rext, (alphaVars, RB))
        ) else if #alphaVars > 0 then (
            -- Particular α:=0 solution, projected back to Rw.
            finalSubst := map(Rw, Rext,
                apply(numgens Rw, i -> Rw_i) | apply(#alphaVars, k -> 0_Rw));
            matrix apply(numRows Hsym, i ->
                apply(numColumns Hsym, j -> finalSubst(Hsym_(i,j))))
        ) else matrix Hsym;
        H0w := liftGreedyHToH0(Hfinal, aw, BlistGr, nFGr);
        if not pinToModP then return H0w;
        -- Lift H0w from Rw (ZZ/p coefficients) back to R (QQ coefficients).
        -- Coefficients map via ZZ/p → ZZ (representative in [0, p)) → QQ.
        -- Monomials map via variable-to-variable.
        liftEntry := p -> (
            if p == 0_Rw then 0_R
            else (
                (mons, coeffs) := coefficients p;
                monsR := sub(mons, R);
                coeffsR := matrix apply(numRows coeffs, i -> apply(numColumns coeffs, j ->
                    sub(lift(coeffs_(i, j), ZZ), QQ)));
                first first entries (monsR * coeffsR)
            )
        );
        return matrix apply(numRows H0w, i -> apply(numColumns H0w, j -> liftEntry (H0w_(i,j))));
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
    else (error "Strategy not supported. getH0 accepts null (Default) / \"Larsson\" / \"Greedy\". For the particular-α:=0 mode (formerly Strategy => \"MatrixHi\") pass Strategy => \"Greedy\", AdjustParams => false.")
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

getTemplate = method(Options => {MonomialOrder => null, Strategy => null, AdjustParams => true})
getTemplate(EliminationTemplate) := o -> E -> (
    if E.cache#?"monomialPartition" and (o.Strategy === null or (E.cache#?"lastPartitionStrategy" and E.cache#"lastPartitionStrategy" === o.Strategy)) then (
        (E.cache#"shifts", E.cache#"monomialPartition")
    ) else (
        J := ideal E;
        R := ring J;
        a := sub(actionVariable E, R);

        -- Direction A: monomial action goes through [excess | residual | basis]
        -- with extractActionFromTemplate. No graph-ring lift, no (s - a) * b_k
        -- shifts, template matches Martyushev paper sizes.
        --
        -- Increment 3: polynomial action lifts to R_s = R[s] with the extra
        -- generator s - a, so that s IS a ring variable on R_s/Is (monomial
        -- action on Is in R_s). The same [excess | residual | basis] pipeline
        -- then handles both cases uniformly.
        if isMonomialAction a then (
            B := lift(basis(R/J), R);
            (shOrig, mpHelper) := getTemplateHelper(a, B, J, o);

            BlistMono := flatten entries B;
            BsetMono := set BlistMono;
            resMonsMono := select(apply(BlistMono, b -> a * b), m -> not BsetMono#?m);
            allMonsMono := union(set \ flatten \ entries \ monomials \ shiftPolynomials(shOrig, J));
            monsBMono := BsetMono;
            monsRMono := set resMonsMono;
            monsEMono := allMonsMono - monsBMono - monsRMono;
            mpMono := new MonomialPartition from rsort \ toList \ {monsEMono, monsRMono, monsBMono};

            E.cache#basis = B;
            E.cache#"shifts" = shOrig;
            E.cache#"monomialPartition" = mpMono;
            E.cache#"lastPartitionStrategy" = o.Strategy;
            E.cache#"isMonomialAction" = true;
            -- basisList / residualMonomials MUST match the template's B- and R-column
            -- ordering (mpMono#2 and mpMono#1 respectively), otherwise
            -- extractActionFromTemplate reads the wrong columns.
            E.cache#"basisList" = mpMono#2;
            E.cache#"residualMonomials" = mpMono#1;
            return (shOrig, mpMono);
        );

        -- Polynomial action: lift to R_s = R[s] / <s - a>, then drive the same
        -- [excess | residual | basis] pipeline on (s, Is) in R_s. After the
        -- lift s is a ring variable, so isMonomialAction(s) is true and the
        -- Greedy pipeline (buildGapPolys etc.) works on R_s even when the
        -- original action was polynomial on R.
        --
        -- Follow origin's shift construction: compute shifts from getH0 on the
        -- ORIGINAL (R, J) and lift them into R_s, then append |B| shifts for
        -- the (s - a) generator (one per basis monomial). Building shifts on
        -- the LIFTED ideal directly yields a much larger template whose RREF
        -- is numerically fragile during solution recovery.
        B := lift(basis(R/J), R);
        (shOrig, mpHelperOrig) := getTemplateHelper(a, B, J, o);

        K := coefficientRing R;
        ringVars := flatten entries vars R;
        MO := if not instance(o.MonomialOrder, Nothing) then o.MonomialOrder else (options R).MonomialOrder;
        Rs := K[prepend("s", ringVars), MonomialOrder => {Eliminate 1, MO}];

        toRs := map(Rs, R, apply(numgens R, i -> Rs_(i+1)));
        aS := toRs(a);
        actVarS := Rs_0;
        Is := ideal(toRs(gens J) | matrix{{actVarS - aS}});
        Bs := toRs(B);
        sortedBs := rsort flatten entries Bs;

        shLifted := new ShiftSet from (
            apply(shOrig, sh -> toRs(sh)) | {matrix {sortedBs}}
        );

        allMonsLifted := union(set \ flatten \ entries \ monomials \ shiftPolynomials(shLifted, Is));
        monsBLifted := set sortedBs;
        monsRLifted := set apply(sortedBs, b -> actVarS * b);
        monsELifted := allMonsLifted - monsBLifted - monsRLifted;
        mpLifted := new MonomialPartition from rsort \ toList \ {monsELifted, monsRLifted, monsBLifted};

        E.cache#basis = Bs;
        E.cache#"graphIdeal" = Is;
        E.cache#"shifts" = shLifted;
        E.cache#"monomialPartition" = mpLifted;
        E.cache#"lastPartitionStrategy" = o.Strategy;
        E.cache#"isMonomialAction" = true;
        E.cache#"liftedToRs" = true;
        E.cache#"liftedActionVar" = actVarS;
        E.cache#"basisList" = mpLifted#2;
        E.cache#"residualMonomials" = mpLifted#1;
        (shLifted, mpLifted)
    )
)

copyTemplate = method(Options => {})
copyTemplate(EliminationTemplate, Ideal) := o -> (E, J) -> (
    Rnew := ring J;
    FFnew := coefficientRing Rnew;
    aNew := sub(actionVariable E, Rnew);
    Enew := eliminationTemplate(aNew, J);

    -- Direction A monomial-direct (non-lifted): basis / shifts live in R,
    -- transplant via sub(..., Rnew). Increment 3 lifted templates have
    -- isMonomialAction = true too but their cache is in R_s, so they route
    -- through the graph-ideal branch below (which handles R_s transplant).
    if E.cache#?"isMonomialAction" and E.cache#"isMonomialAction"
        and not (E.cache#?"liftedToRs" and E.cache#"liftedToRs") then (
        Enew.cache#"isMonomialAction" = true;
        if E.cache#?basis then Enew.cache#basis = sub(E.cache#basis, Rnew);
        if E.cache#?"shifts" then Enew.cache#"shifts" = new ShiftSet from apply(E.cache#"shifts", sh -> sub(sh, Rnew));
        if E.cache#?"monomialPartition" then Enew.cache#"monomialPartition" =
            new MonomialPartition from apply(E.cache#"monomialPartition",
                mp -> apply(mp, m -> sub(m, Rnew)));
        if E.cache#?"basisList" then Enew.cache#"basisList" =
            apply(E.cache#"basisList", b -> sub(b, Rnew));
        if E.cache#?"residualMonomials" then Enew.cache#"residualMonomials" =
            apply(E.cache#"residualMonomials", m -> sub(m, Rnew));
        if E.cache#?"lastPartitionStrategy" then Enew.cache#"lastPartitionStrategy" = E.cache#"lastPartitionStrategy";
        if E.cache#?"lastMatrixStrategy" then Enew.cache#"lastMatrixStrategy" = E.cache#"lastMatrixStrategy";
        if E.cache#?"lastActionStrategy" then Enew.cache#"lastActionStrategy" = E.cache#"lastActionStrategy";
        if E.cache#?"shifts" and E.cache#?"monomialPartition" then (
            shiftsNew := Enew.cache#"shifts";
            mpNew := Enew.cache#"monomialPartition";
            allMonsMatNew := mpNew#0 | mpNew#1 | mpNew#2;
            Enew.cache#"templateMatrix" = sub(
                transpose fold(
                    apply(shiftPolynomials(shiftsNew, J),
                        m -> last coefficients(m, Monomials => allMonsMatNew)),
                    (u, v) -> u | v),
                coefficientRing Rnew);
        );
        return Enew;
    );

    if E.cache#?"graphIdeal" then (
	Rs := ring E.cache#"graphIdeal";
        Rsnew := FFnew[gens ring E.cache#"graphIdeal", MonomialOrder => (options Rs).MonomialOrder];

	if E.cache#?basis then Enew.cache#basis = sub(if E.cache#?basis then E.cache#basis else basis E, Rsnew);
        if E.cache#?"shifts" then Enew.cache#"shifts" = apply(E.cache#"shifts", sh -> sub(sh, Rsnew));
        if E.cache#?"monomialPartition" then Enew.cache#"monomialPartition" = apply(E.cache#"monomialPartition", mp -> apply(mp, m -> sub(m, Rsnew)));
        if E.cache#?"lastPartitionStrategy" then Enew.cache#"lastPartitionStrategy" = E.cache#"lastPartitionStrategy";
        if E.cache#?"lastMatrixStrategy" then Enew.cache#"lastMatrixStrategy" = E.cache#"lastMatrixStrategy";
	if E.cache#?"lastActionStrategy" then Enew.cache#"lastActionStrategy" = E.cache#"lastActionStrategy";

        -- Increment 3: propagate lifted-polynomial-action flags + basisList / residualMonomials markers.
        if E.cache#?"isMonomialAction" then Enew.cache#"isMonomialAction" = E.cache#"isMonomialAction";
        if E.cache#?"liftedToRs" then Enew.cache#"liftedToRs" = E.cache#"liftedToRs";
        if E.cache#?"liftedActionVar" then Enew.cache#"liftedActionVar" = sub(E.cache#"liftedActionVar", Rsnew);

	toRsnew := map(Rsnew, Rnew, apply(numgens Rnew, i -> Rsnew_(i+1)));
        JsGens := toRsnew(gens J);
        aS := toRsnew(aNew);
        actVar := Rsnew_0;

        Enew.cache#"graphIdeal" = ideal(JsGens | matrix{{actVar - aS}});

        if Enew.cache#?"monomialPartition" then (
            mpNewCopy := Enew.cache#"monomialPartition";
            Enew.cache#"basisList" = mpNewCopy#2;
            Enew.cache#"residualMonomials" = mpNewCopy#1;
        );

	if E.cache#?"templateMatrix" then (
            isLifted := E.cache#?"liftedToRs" and E.cache#"liftedToRs";
            shNew := Enew.cache#"shifts";
            mpNew := Enew.cache#"monomialPartition";
            IsNew := Enew.cache#"graphIdeal";
            Enew.cache#"templateMatrix" = if isLifted then (
                -- Increment 3 [E|R|B] layout: include residual block.
                allMonsMatL := mpNew#0 | mpNew#1 | mpNew#2;
                sub(
                    transpose fold(
                        apply(shiftPolynomials(shNew, IsNew),
                            m -> last coefficients(m, Monomials => allMonsMatL)),
                        (u, v) -> u | v),
                    coefficientRing Rsnew)
            ) else (
                -- Origin's [E|B] layout (unreachable in current code path, kept
                -- as the fallback target for a future LU-on-RREF-error guard).
                getTemplateMatrix(shNew, mpNew, IsNew)
            );
        );
    );
    Enew
)

getTemplateMatrix = method(Options => {MonomialOrder => null, Strategy => null, AdjustParams => true})
getTemplateMatrix(RingElement, Matrix, Ideal) := o -> (a, B, J) -> (
    getTemplateMatrix(eliminationTemplate(a, J), o)
)
getTemplateMatrix(ShiftSet, MonomialPartition, Ideal) := o -> (shifts, monomialPartition, J) -> (
    allMons := monomialPartition#0 | monomialPartition#2;
    sub(transpose fold(apply(shiftPolynomials(shifts, J), m -> last coefficients(m, Monomials => allMons)), (a,b) -> a|b), coefficientRing ring J)
)
getTemplateMatrix(EliminationTemplate) := o -> E -> (
    -- Unified path post-Increment-3: every strategy flows through
    -- getTemplate -> [E | R | B] fold. Strategy choice only affects the H0
    -- computation inside getH0 (Default: GB; Larsson: H0 mod syz F; Greedy:
    -- buildGapPolys + buildHSymbolic + adjustParams or alpha:=0). The
    -- column layout is always [E | R | B] for both monomial-direct and
    -- polynomial-lifted templates.
    if E.cache#?"templateMatrix" and (o.Strategy === null or (E.cache#?"lastMatrixStrategy" and E.cache#"lastMatrixStrategy" === o.Strategy)) then (
        E.cache#"templateMatrix"
    ) else (
        (shifts, monomialPartition) := getTemplate(E, o);
        -- The shifts / partition were built against graphIdeal = J or
        -- J_s = J + <s - a> (lifted case); use that ideal so generator
        -- counts match shiftPolynomials' assertion.
        idealForFold := if E.cache#?"graphIdeal" then E.cache#"graphIdeal" else ideal E;
        allMonsMat := monomialPartition#0 | monomialPartition#1 | monomialPartition#2;
        ret := sub(
            transpose fold(
                apply(shiftPolynomials(shifts, idealForFold),
                    m -> last coefficients(m, Monomials => allMonsMat)),
                (u, v) -> u | v),
            coefficientRing ring idealForFold);
        E.cache#"templateMatrix" = ret;
        E.cache#"lastMatrixStrategy" = o.Strategy;
        ret
    )
)

getActionMatrix = method(Options => {MonomialOrder => null, Strategy => null, AdjustParams => true})
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
    -- Unified action extraction post-Increment-3: RREF + pivot on the
    -- [E | R | B] template for every strategy and every action type. The
    -- action variable is the user's `a` for monomial-direct templates, or
    -- the cached `s` in R_s for polynomial-action templates (lifted to R_s
    -- = R[s] / <s - a>). Fallback to origin's LU split-solve on [E | B] if
    -- RREF raises a catchable error (QQ coefficient blowup on adversarial
    -- input).
    if E.cache#?"actionMatrix" and (o.Strategy === null or (E.cache#?"lastActionStrategy" and E.cache#"lastActionStrategy" === o.Strategy)) then (
        E.cache#"actionMatrix"
    ) else (
        getTemplate(E, o);  -- populates partition / basisList / residualMonomials
        templateMatrix := getTemplateMatrix(E, o);

        actVarForExtract := if E.cache#?"liftedToRs" and E.cache#"liftedToRs"
                            then E.cache#"liftedActionVar"
                            else actionVariable E;
        ret := try (
            extractActionFromTemplate(templateMatrix, E.cache#"residualMonomials",
                E.cache#"basisList", actVarForExtract)
        ) else (
            stderr << "-- [getActionMatrix] RREF failed on [E|R|B]; "
                   << "falling back to LU split-solve on [E|B]" << endl;
            mp1 := E.cache#"monomialPartition";
            shiftsNow := E.cache#"shifts";
            idealNow := if E.cache#?"graphIdeal" then E.cache#"graphIdeal" else ideal E;
            mpFallback := new MonomialPartition from {mp1#0 | mp1#1, {}, mp1#2};
            matFallback := getTemplateMatrix(shiftsNow, mpFallback, idealNow, o);
            getActionMatrix(actVarForExtract, mpFallback, matFallback, o)
        );
        E.cache#"actionMatrix" = ret;
        E.cache#"lastActionStrategy" = o.Strategy;
        ret
    )
)

-*

*-
getEigenMatrix = method(Options => {MonomialOrder => null, Strategy => null, AdjustParams => true})
getEigenMatrix(EliminationTemplate) := o -> (E) -> (
    Ma := getActionMatrix(E, o);
    -- Monomial-direct templates let downstream recoverSolutionsMatrixHi do
    -- its own normalization by the "1" slot, so return raw eigenvectors
    -- against the cached basis list. Polynomial-lifted templates (and any
    -- path that falls through to recoverSolutions) need eigenvectors
    -- pre-normalized so monomial values can be read directly.
    if (E.cache#?"isMonomialAction" and E.cache#"isMonomialAction")
        and not (E.cache#?"liftedToRs" and E.cache#"liftedToRs") then (
        (svals, P) := eigenvectors sub(Ma, CC);
        (matrix{E.cache#"basisList"}, P)
    ) else (
        (svals2, P2) := eigenvectors sub(Ma, CC);
        cleanEvecs := clean_(1e-10) (P2 * inverse diagonalMatrix(P2^{numColumns P2 - 1}));
        (transpose rsort basis E, cleanEvecs)
    )
)
getEigenMatrix(Ideal) := o -> (I) -> getEigenMatrix(random(1, ring I), I, o)
getEigenMatrix(RingElement, Ideal) := o -> (a, J) -> (
    E := eliminationTemplate(a, J);
    getEigenMatrix(E, o)
)

templateSolve = method(Options => {MonomialOrder => null, Strategy => null, AdjustParams => true})
templateSolve(EliminationTemplate) := o -> (E) -> (
    -- Unified solve post-Increment-3. Compute the action matrix Ma (RREF
    -- + pivot on [E | R | B]), eigendecompose over CC. If every ring
    -- variable appears in the basis, read coordinates directly from
    -- basis-indexed eigenvector slots (recoverSolutionsMatrixHi). Otherwise
    -- fall back to recoverSolutions, which uses RREF on the template to
    -- express excess / residual monomials in terms of B.
    Ma := getActionMatrix(E, o);
    R := ring ideal E;

    -- Monomial-direct: basis is in R, all vars present by construction
    -- (action is a ring variable, basis includes 1 and extends).
    if (E.cache#?"isMonomialAction" and E.cache#"isMonomialAction")
        and not (E.cache#?"liftedToRs" and E.cache#"liftedToRs") then (
        return recoverSolutionsMatrixHi(Ma, E.cache#"basisList", R);
    );

    -- Polynomial-lifted case: basis is in R_s, map back via s -> 0 to
    -- get a basis in R. If every ring variable is a basis monomial,
    -- recoverSolutionsMatrixHi works. Otherwise recoverSolutions reads
    -- missing vars from the template's E / R blocks.
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

    -- Use RREF on the full [E|R|B] template to express each excess monomial
    -- in terms of B. This handles over-determined templates (numTop > |E|+|R|)
    -- naturally, avoids QQ least-squares (not implemented), and matches the
    -- approach extractActionFromTemplate uses for residual columns.
    FF := coefficientRing Rnew;
    Rref := reducedRowEchelonForm templateMat;
    -- eReduction_(colIdx) = row in Rref whose pivot is at colIdx, or null.
    pivotRowOfCol := new MutableHashTable;
    for r from 0 to numRows Rref - 1 do (
        pivCol := position(0..numColumns Rref - 1, c -> Rref_(r,c) != 0_FF);
        if pivCol =!= null and not pivotRowOfCol#?pivCol then
            pivotRowOfCol#pivCol = r;
    );
    -- readMonValFromRref(col, basisVals) returns the CC value of the monomial
    -- in column `col` given the CC values of B monomials: row at `col` says
    --   1 * m_col + sum_b Rref[row, numE+numR+b] * b = 0
    -- so m_col = -sum_b Rref[row, numE+numR+b] * b_val.
    readMonValFromRref := (col, bVals) -> (
        if not pivotRowOfCol#?col then return null;
        row := pivotRowOfCol#col;
        val := 0_CC;
        for k from 0 to numB - 1 do
            val = val - sub(Rref_(row, numE + numR + k), CC) * bVals#k;
        val
    );
    
    solutions := {};
    varsList := flatten entries vars Rnew;

    -- Basis slot holding the monomial 1_Rnew. Eigenvectors normalize by this
    -- slot (done upstream by getEigenMatrix), so when this slot is tiny the
    -- eigenvector is spurious / numerically degenerate and must be skipped
    -- rather than blown up by the normalization.
    oneSlot := position(basisMonsRnew, b -> b == 1_Rnew);

    for rootIndex from 0 to numColumns M - 1 do (
        if oneSlot =!= null and abs M_(oneSlot, rootIndex) < 1e-10 then continue;

        monomialValues := new MutableHashTable;
        for i from 0 to #basisMonsRnew - 1 do (
            monomialValues#(basisMonsRnew#i) = M_(i, rootIndex);
        );

        -- Basis-monomial values as a vector aligned with monsB.
        bVals := apply(numB, j -> monomialValues#(toRnew(monsB#j)));

        root := {};
        for v in varsList do (
            if monomialValues#?v then (
                -- v is a basis monomial; read directly.
                root = append(root, monomialValues#v);
            )
            else (
                -- v lives in E or R block; reduce via RREF pivot row.
                vRs := toRs(v);
                posInE := position(monsE, m -> m == vRs);
                posInR := if posInE === null then position(monsR, m -> m == vRs) else null;
                colIdx := if posInE =!= null then posInE
                          else if posInR =!= null then numE + posInR
                          else null;
                if colIdx =!= null then (
                    val := readMonValFromRref(colIdx, bVals);
                    if val === null then (
                        -- RREF has no pivot at this column; shouldn't happen
                        -- for a well-posed template but fall back to the
                        -- polynomial normal form against the ideal.
                        r := v % J;
                        coeffs := last coefficients(r, Monomials => basisMonsRnew);
                        val = 0_CC;
                        for j from 0 to numB - 1 do
                            if monomialValues#?(basisMonsRnew#j) then
                                val = val + sub(coeffs_(j,0), CC) * monomialValues#(basisMonsRnew#j);
                    );
                    root = append(root, val);
                )
                else (
                    -- Failsafe: v is neither in B nor in E nor in R.
                    r := v % J;
                    coeffs := last coefficients(r, Monomials => basisMonsRnew);
                    val := 0_CC;
                    for j from 0 to numB - 1 do
                        if monomialValues#?(basisMonsRnew#j) then
                            val = val + sub(coeffs_(j,0), CC) * monomialValues#(basisMonsRnew#j);
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
-- Smoke test of the basic pipeline (getTemplateMatrix / getActionMatrix).
-- Template-size checks for this system are in the Greedy + AdjustParams=>false
-- TEST below.
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

TEST /// -- Greedy + AdjustParams=>false (alpha := 0 particular solution):
-- paper-reference sizes. Under Direction A all strategies already match
-- paper-reference sizes on monomial-action problems; this test pins the
-- alpha:=0 sizes specifically.
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  ET = eliminationTemplate(y, I);
  M = getTemplateMatrix(ET, Strategy => "Greedy", AdjustParams => false);
  assert(numRows M == 10 and numColumns M == 20)

  J = ideal(x^3+y^3+z^3-4,x^2-y-z-1,x-y^2+z-3)
  ET2 = eliminationTemplate(x, J);
  N = getTemplateMatrix(ET2, Strategy => "Greedy", AdjustParams => false);
  -- Greedy + AdjustParams=>false on this problem uses a min-degree
  -- buildHSymbolic H, size 20x33. Default's GB H0 happens to have smaller
  -- monomial support here (16x28 under Direction A); which mode wins is
  -- problem-dependent.
  assert(numRows N == 20 and numColumns N == 33)
///

TEST /// -- Direction A: on 5pt essential (monomial action, 0 free alphas),
-- Greedy flows through the unified getH0 but skips the (s-a) graph lift,
-- so its template matches paper-reference 10x20 (not the old 20x20).
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  ET = eliminationTemplate(y, I);
  Mg = getTemplateMatrix(ET, Strategy => "Greedy");
  assert(numRows Mg == 10 and numColumns Mg == 20)
///

TEST /// -- getActionMatrix on Greedy templates (alpha:=0 and adjustParams
-- modes): 10x10 with deg I distinct eigenvalues on 5pt essential.
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  d = degree I;
  E1 = eliminationTemplate(y, I);
  Ma1 = getActionMatrix(E1, Strategy => "Greedy", AdjustParams => false);
  assert(numRows Ma1 == d and numColumns Ma1 == d);
  assert(#eigenvalues sub(Ma1, CC) == d);
  E2 = eliminationTemplate(y, I);
  Ma2 = getActionMatrix(E2, Strategy => "Greedy");
  assert(numRows Ma2 == d and numColumns Ma2 == d);
  assert(#eigenvalues sub(Ma2, CC) == d);
///

TEST /// -- templateSolve via Greedy + AdjustParams=>false on 5pt essential
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  d = degree I;
  E1 = eliminationTemplate(y, I);
  sols1 = templateSolve(E1, Strategy => "Greedy", AdjustParams => false);
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

TEST /// -- Greedy + AdjustParams=>false with action variable x on 5pt
-- essential. Direction A handles the monomial action via the unified
-- [E|R|B] pipeline + recoverSolutionsMatrixHi; the legacy graph-ring
-- templateSolve path used to fail on 3-var systems with non-random actions.
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  d = degree I;
  E = eliminationTemplate(x, I);
  sols = templateSolve(E, Strategy => "Greedy", AdjustParams => false);
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

TEST /// -- Cross-validation under Direction A: on a 0-free-alpha problem
-- (5pt essential), every strategy skips the (s-a) graph lift on monomial
-- actions and uses the shared [excess|residual|basis] layout. Greedy's
-- buildHSymbolic H and Default's GB H0 happen to produce the same shifts
-- on this problem, so Mg == Md at 10x20. Greedy + AdjustParams=>false
-- (alpha:=0) reaches the same size 10x20 too; the action matrices of the
-- three share the char polynomial.
  R = QQ[x,y,z]
  Es = apply(4, i -> random(QQ^3, QQ^3));
  Ee = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
  I = ideal(Ee*transpose Ee * Ee - (1/2) * trace(Ee * transpose Ee) * Ee) + ideal(det Ee);
  E1 = eliminationTemplate(y, I);
  Mh = getTemplateMatrix(E1, Strategy => "Greedy", AdjustParams => false);
  E2 = eliminationTemplate(y, I);
  Mg = getTemplateMatrix(E2, Strategy => "Greedy");
  E3 = eliminationTemplate(y, I);
  Md = getTemplateMatrix(E3);
  Ah = getActionMatrix(E1, Strategy => "Greedy", AdjustParams => false);
  Ag = getActionMatrix(E2, Strategy => "Greedy");
  Ad = getActionMatrix(E3);
  -- Every strategy uses Direction A's downstream on monomial action
  -- and produces identical templates on this 0-free-alpha problem.
  assert(Mg == Md);
  assert(numRows Mh == numRows Mg);
  assert(Ag == Ad);
  -- Char polynomial of action is invariant under basis permutation.
  Rt = QQ[t];
  chiH = det(t * id_(Rt^(numRows Ah)) - sub(Ah, Rt));
  chiG = det(t * id_(Rt^(numRows Ag)) - sub(Ag, Rt));
  assert(chiH == chiG);
///

TEST /// -- Greedy on 6 Demazure cubics (no det): free alpha > 0,
-- adjustParams commits alphas to cancel excessive monomials, producing a
-- smaller H0 (fewer monomials per row) than Default's Groebner H0. Both
-- flow through the shared Direction A pipeline: Greedy's template should
-- have row count <= Default's (fewer shifts) and both solve to residual
-- < 1e-6.
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
  solsH = templateSolve(E3, Strategy => "Greedy", AdjustParams => false);
  solsD = templateSolve(E1);
  solsG = templateSolve(E2, Strategy => "Greedy");
  -- Greedy invariant: Greedy <= Default in rows on problems with free alphas.
  assert(numRows Mg <= numRows Md);
  assert(#solsD == d and #solsG == d and #solsH == d);
  Rc = CC[gens R];
  gensC = sub(gens I, Rc);
  assert(all(solsD, s -> 1e-6 > norm sub(gensC, matrix{s})));
  assert(all(solsG, s -> 1e-6 > norm sub(gensC, matrix{s})));
  assert(all(solsH, s -> 1e-6 > norm sub(gensC, matrix{s})));
///

TEST /// -- Greedy + AdjustParams=>false on the 3-variable docs system.
  R = QQ[x,y,z]
  J = ideal(x^3+y^3+z^3-4, x^2-y-z-1, x-y^2+z-3)
  d = degree J;
  E1 = eliminationTemplate(x, J);
  solsH = templateSolve(E1, Strategy => "Greedy", AdjustParams => false);
  assert(#solsH == d);
  assert(all(solsH, s -> 1e-6 > norm sub(sub(gens J, CC[gens R]), matrix{s})));
///

TEST /// -- Paper §3 small example on a 2-var system (deg I = 12): Greedy +
-- AdjustParams=>false produces a 7×19 template and solves to residual < 1e-6.
  R = QQ[x,y]
  I = ideal(x^4+y^2+x*y-3, x^2*y+y^3-2)
  d = degree I;
  E = eliminationTemplate(x, I);
  M = getTemplateMatrix(E, Strategy => "Greedy", AdjustParams => false);
  assert(numColumns M == 19);
  assert(numRows M < 19);  -- strictly smaller than Default's 19×19
  sols = templateSolve(E, Strategy => "Greedy", AdjustParams => false);
  assert(#sols == d);
  assert(all(sols, s -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{s})));
///

TEST /// -- #1 F+λ 8pt from Martyushev CVPR 2022 Table 1 over ZZ/32749.
-- Paper target: std 11×19, nstd 7×15. Greedy + AdjustParams=>false achieves
-- 11×20 (matches std up to one column). No tier-3 since we're in ZZ/p.
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
  M = getTemplateMatrix(E, Strategy => "Greedy", AdjustParams => false);
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

TEST /// -- copyTemplate on a Greedy-built template. Greedy populates the
-- same monomial partition / shifts cache as Default, so copyTemplate
-- transplants the structure into a perturbed ideal and the Greedy solve
-- on the copy still hits tier-3 residual.
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
  -- Default / Larsson / Greedy all flow through the shared Direction A
  -- [E|R|B] pipeline; only the H0 computation inside getH0 differs per
  -- strategy.
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
-- Strategy semantics (post-Increment-3 unified pipeline):
--   Every strategy flows through Direction A's [E|R|B] + RREF+pivot
--   extraction. Choice of Strategy only changes how H0 is computed inside
--   getH0:
--     null (Default)                    : H0 from Groebner change-of-basis
--     "Larsson"                         : H0 mod syz(F) (CVPR 2017)
--     "Greedy", AdjustParams => false   : alpha := 0 particular solution
--                                         of buildHSymbolic ("MatrixHi mode")
--     "Greedy"                          : buildHSymbolic + adjustParams
--                                         (Martyushev CVPR 2022 §4);
--                                         ZZ/p-pinned inner RREF on QQ input
--
-- Polynomial actions lift to R_s = R[s] / <s - a> and run the same
-- pipeline on (s, J_s) in R_s.
--
-- templateSolve picks recoverSolutionsMatrixHi when every ring variable is
-- a basis monomial (monomial-direct, or lifted problems whose basis happens
-- to contain all original vars) and falls back to recoverSolutions (RREF +
-- pivot) otherwise.
--------------------------------------------------------------------------

-- Demo 1: 5-point essential matrix (Demazure trace identity + det)
-- deg I = 10. Under Direction A, every strategy reaches paper std size
-- 10x20 on this monomial-action problem.
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
getTemplateMatrix ET                          -- Default   : 10 x 20
getTemplateMatrix(ET, Strategy => "Larsson")  -- Larsson   : 10 x 20
getTemplateMatrix(ET, Strategy => "Greedy", AdjustParams => false) -- alpha:=0 : 10 x 20  <- paper
getTemplateMatrix(ET, Strategy => "Greedy")   -- Greedy    : 10 x 20  (0 free alphas)

-- End-to-end solve + tier-3 residual check via Greedy.
ET2 = eliminationTemplate(y, I)
sols = templateSolve(ET2, Strategy => "Greedy")
assert(#sols == degree I)
assert(all(sols, s -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{s})))


-- Demo 2: 6 Demazure cubics (no det) -- adjustParams has leverage here.
-- 66 free alphas; adjustParams commits some of them to cut shifts in H0.
-- Every strategy shares Direction A's [E|R|B] pipeline. Default 27x34,
-- Larsson / Greedy 24x34 (cut via adjustParams or syzygy reduction),
-- Greedy + AdjustParams=>false (alpha:=0) 25x35.
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
Md = getTemplateMatrix(E1)                          -- Default   : 27 x 34
E2 = eliminationTemplate(y, I)
Mg = getTemplateMatrix(E2, Strategy => "Greedy")    -- Greedy    : 24 x 34  <- shift reduction
assert(numRows Mg <= numRows Md)
E3 = eliminationTemplate(y, I)
Mh = getTemplateMatrix(E3, Strategy => "Greedy", AdjustParams => false)  -- alpha:=0 : 25 x 35

-- Both solve correctly to residual < 1e-6.
solsG = templateSolve(E2, Strategy => "Greedy")
assert(#solsG == degree I)
assert(all(solsG, s -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{s})))


-- Demo 3: paper §3 example (2-var system, deg I = 12).
-- Direction A + Greedy + AdjustParams=>false reaches 7x19 (paper target).
restart
path = prepend("./", path)
needsPackage "EliminationTemplates"
R = QQ[x,y]
I = ideal(x^4 + y^2 + x*y - 3, x^2*y + y^3 - 2)
degree I  -- 12

E = eliminationTemplate(x, I)
getTemplateMatrix E                               -- Default   :  7 x 19
getTemplateMatrix(E, Strategy => "Greedy", AdjustParams => false)      -- alpha:=0  :  7 x 19
sols = templateSolve(E, Strategy => "Greedy", AdjustParams => false)
assert(#sols == degree I)
assert(all(sols, s -> 1e-6 > norm sub(sub(gens I, CC[gens R]), matrix{s})))

-- Cross-check: eigenvalues of the action matrix == x-coordinates of solutions.
Ma = getActionMatrix(E, Strategy => "Greedy", AdjustParams => false)
evals = eigenvalues sub(Ma, CC)
assert(#evals == degree I)


-- Demo 4: Paper Table 1 #1 (F+lambda 8pt) over ZZ/p; matches Martyushev std.
-- Shown over ZZ/32749 because Greedy over QQ is slower (coefficient growth).
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
getTemplateMatrix(E, Strategy => "Greedy", AdjustParams => false)  -- 11 x 20 (paper std: 11 x 19)




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
