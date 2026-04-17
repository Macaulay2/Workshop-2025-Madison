-- bench_greedy_effect.m2
-- Demonstrates that Strategy => "Greedy" (= MatrixHi + adjustParams,
-- Martyushev CVPR 2022) produces a strictly smaller template than
-- Strategy => "MatrixHi" (alphas := 0) on a problem with free α > 0.
--
-- Problem: 6 Demazure cubics of a random 3-linear essential matrix E(x,y,z)
-- (i.e. the 5-point essential problem without the det constraint). Over QQ
-- with seed 42, this ideal is 0-dim with degree 10 — same solutions as the
-- full 5pt essential, just with a redundant generator set. The MatrixHi
-- particular solution leaves 66 free α's; adjustParams commits some of them
-- to cancel excessive monomials.
--
-- Invoke from anywhere:
--     M2 /path/to/workshop/tests/bench_greedy_effect.m2

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
Em = x*Es#0 + y*Es#1 + z*Es#2 + Es#3
-- 6 Demazure cubics only (no det constraint). Still 0-dim with degree 10.
I5 = ideal(2*Em*transpose(Em)*Em - trace(Em*transpose(Em))*Em)
aVar = y
dI = degree I5

<< "========================================================================" << endl
<< "  Greedy effect: 6 Demazure cubics (5pt essential without det)" << endl
<< "  action variable = " << aVar << "    deg(R/I) = " << dI << endl
<< "========================================================================" << endl << endl

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
    t2 := (nev == dI);
    t3 := (res < 1e-6);
    << "   template size: " << numRows M << "x" << numColumns M << endl;
    << "   build+solve:   " << elapsed << "s" << endl;
    << "   #eigenvalues:  " << nev << " (expected " << dI << ")" << endl;
    << "   max |F(x*)|:   " << toString res << endl;
    << "   tier-2:        " << (if t2 then "PASS" else "FAIL") << endl;
    << "   tier-3:        " << (if t3 then "PASS" else "FAIL") << endl << endl;
    (numRows M, numColumns M, elapsed, nev, res, t2, t3)
)

(rH, cH, tH, evH, resH, t2H, t3H) = runStrategy "MatrixHi";
(rG, cG, tG, evG, resG, t2G, t3G) = runStrategy "Greedy";

<< "========================================================================" << endl
<< "Delta (MatrixHi - Greedy): " << (rH - rG) << " rows, " << (cH - cG) << " cols" << endl
<< "========================================================================" << endl

if (rH - rG) > 0 or (cH - cG) > 0 then (
    << "adjustParams REDUCED the template — this is the greedy's value." << endl;
) else (
    << "adjustParams did not reduce the template on this problem." << endl;
);

exit 0
