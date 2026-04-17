-- bench_5pt_all.m2
-- Benchmark all EliminationTemplates strategies on the 5-point essential
-- matrix problem (Demazure cubics + det, deg I = 10).
--
-- For each strategy reports:
--   - template size (rows x cols)
--   - build time (cpuTime, seconds)
--   - #eigenvalues of action matrix          [tier-2]
--   - max |F_j(x*)| over recovered solutions [tier-3]
--
-- Strategies exercised:
--   Larsson   (H0 mod syzygy module, CVPR 2017)
--   MatrixHi  (per-monomial alphas, alphas := 0)
--   Greedy    (MatrixHi + adjustParams, Martyushev CVPR 2022)
--
-- Invoke from anywhere:
--     M2 /path/to/workshop/tests/bench_5pt_all.m2

pkgDir = currentFileDirectory | "../"
if not fileExists(pkgDir | "EliminationTemplates.m2") then (
    << "ERROR: EliminationTemplates.m2 not found at " << pkgDir << endl;
    exit 1;
);
path = prepend(pkgDir, path)
needsPackage("EliminationTemplates", Reload => true)

-- ========== Problem setup ==========
setRandomSeed 42
R = QQ[x,y,z]
Es = apply(4, i -> random(QQ^3, QQ^3))
E5 = x*Es#0 + y*Es#1 + z*Es#2 + Es#3
I5 = ideal(2*E5*transpose(E5)*E5 - trace(E5*transpose(E5))*E5, det E5)
aVar = y
dI = degree I5

<< "========================================================================" << endl
<< "  EliminationTemplates benchmark: 5-point essential matrix" << endl
<< "  action variable = " << aVar << "    deg(R/I) = " << dI << endl
<< "  tier-2 PASS iff #eigenvalues == " << dI << endl
<< "  tier-3 PASS iff max |F_j(x*)| < 1e-6" << endl
<< "========================================================================" << endl << endl

-- Substitute solutions into the original generators, return max residual.
evalResid = (sols, I) -> (
    if #sols == 0 then return 1e99;
    Rc := CC[gens ring I];
    gensC := sub(gens I, Rc);
    mx := 0.0;
    scan(sols, s -> (
        sM := matrix{apply(s, v -> sub(v, CC))};
        n := norm sub(gensC, sM);
        if n > mx then mx = n;
    ));
    mx
)

runStrategy = strategyName -> (
    << "---- " << strategyName << " ----" << endl;
    ET := eliminationTemplate(aVar, I5);
    t0 := cpuTime();
    M := getTemplateMatrix(ET, Strategy => strategyName);
    Ma := getActionMatrix(ET, Strategy => strategyName);
    sols := templateSolve(ET, Strategy => strategyName);
    elapsed := cpuTime() - t0;
    nev := #eigenvalues sub(Ma, CC);
    res := evalResid(sols, I5);
    sizeStr := toString numRows M | "x" | toString numColumns M;
    t2 := (nev == dI);
    t3 := (res < 1e-6);
    << "   template size: " << sizeStr << endl;
    << "   build+solve:   " << elapsed << "s" << endl;
    << "   #eigenvalues:  " << nev << " (expected " << dI << ")" << endl;
    << "   max |F(x*)|:   " << toString res << endl;
    << "   tier-2:        " << (if t2 then "PASS" else "FAIL") << endl;
    << "   tier-3:        " << (if t3 then "PASS" else "FAIL") << endl << endl;
    {strategyName, sizeStr, toString elapsed, toString nev, toString res,
     if t2 then "PASS" else "FAIL", if t3 then "PASS" else "FAIL"}
)

results = apply({"Larsson", "MatrixHi", "Greedy"}, runStrategy)

-- ========== Summary table ==========
<< "========================================================================" << endl
<< "Summary (5pt essential, deg I = " << dI << ", action var " << aVar << ")" << endl
<< "========================================================================" << endl
hdr = {"Strategy", "Size", "Time(s)", "#ev", "max|F|", "T2", "T3"};
widths = {16, 12, 10, 5, 14, 5, 5};
printRow = r -> (
    scan(#r, i -> << pad(r#i, widths#i) << "  ");
    << endl
);
printRow hdr;
printRow apply(widths, w -> concatenate(w:"-"));
scan(results, printRow);

<< endl << "=== done ===" << endl
exit 0
