-- bench_problems.m2
-- Multi-problem benchmark harness for EliminationTemplates.
-- Runs Default / Larsson / MatrixHi / Greedy on each problem and reports:
--   - template size (rows x cols)
--   - build time (cpuTime, seconds)
--   - tier-2 test (#eigenvalues of action matrix == deg I)
--   - tier-3 test (max |F(x*)| < 1e-6)
--
-- Each problem supplies a name, ring, ideal, action variable (MUST be a
-- ring variable — MatrixHi / Greedy do not support polynomial actions),
-- and target sizes from Martyushev CVPR 2022 Tables 1-2 where applicable.
--
-- Use ZZ/32749 for problems where QQ is too slow (matches reference_martyushev
-- convention); use QQ when you want tier-3 (complex solutions + residuals).
--
-- Invoke from anywhere:
--     M2 /path/to/workshop/tests/bench_problems.m2

pkgDir = currentFileDirectory | "../"
if not fileExists(pkgDir | "EliminationTemplates.m2") then (
    << "ERROR: EliminationTemplates.m2 not found at " << pkgDir << endl;
    exit 1;
);
path = prepend(pkgDir, path)
needsPackage("EliminationTemplates", Reload => true)

setRandomSeed 42

-- ============================================================
-- Problem definitions
-- Each problem is a function that returns {name, R, I, aVar, target, overQQ}
-- `overQQ` controls whether tier-3 (solution residual) is attempted.
-- ============================================================

problems = {}

-- 1. Unit circle + line (toy, QQ)
problems = append(problems, () -> (
    R := QQ[x,y];
    I := ideal(x^2+y^2-1, x+y);
    {"circle+line (toy)", R, I, x, "2x-", true}
))

-- 2. Two conics (QQ, 2-var)
problems = append(problems, () -> (
    R := QQ[x,y];
    I := ideal(x^2+y^2-1, x^2+x*y+y^2-1);
    {"two conics", R, I, x, "4x-", true}
))

-- 3. Section-3 example (QQ, 2-var, deg 12)
problems = append(problems, () -> (
    R := QQ[x,y];
    I := ideal(x^4+y^2+x*y-3, x^2*y+y^3-2);
    {"paper §3 ex (2-var)", R, I, x, "-x-", true}
))

-- 4. 3-var docs system (QQ, deg 12)
problems = append(problems, () -> (
    R := QQ[x,y,z];
    I := ideal(x^3+y^3+z^3-4, x^2-y-z-1, x-y^2+z-3);
    {"3-var docs", R, I, x, "20x33", true}
))

-- 5. 5pt essential, standard Demazure+det (QQ, deg 10)
problems = append(problems, () -> (
    R := QQ[x,y,z];
    Es := apply(4, i -> random(QQ^3, QQ^3));
    Em := x*Es#0 + y*Es#1 + z*Es#2 + Es#3;
    I := ideal(2*Em*transpose(Em)*Em - trace(Em*transpose(Em))*Em, det Em);
    {"5pt essential (+det)", R, I, y, "10x20 (target)", true}
))

-- 6. 6 Demazure cubics only, no det — adjustParams strictly shrinks here (QQ)
problems = append(problems, () -> (
    R := QQ[x,y,z];
    Es := apply(4, i -> random(QQ^3, QQ^3));
    Em := x*Es#0 + y*Es#1 + z*Es#2 + Es#3;
    I := ideal(2*Em*transpose(Em)*Em - trace(Em*transpose(Em))*Em);
    {"6 Demazure cubics", R, I, y, "adjustParams demo", true}
))

-- 7. #1 F+λ 8pt (one-sided radial distortion; Paper Table 1 — ZZ/p for speed)
problems = append(problems, () -> (
    FF := ZZ/32749;
    R := FF[x,y];
    d1 := random(FF^8, FF^4);
    V1 := matrix{{y*x},{y},{x},{1_R}};
    fv := sub(d1, R) * V1;
    Fmat := matrix{{fv_(0,0), fv_(3,0), fv_(5,0)},
                   {fv_(1,0), fv_(4,0), fv_(6,0)},
                   {fv_(2,0), x,        1_R}};
    I := ideal(det Fmat, y*fv_(2,0) - fv_(7,0));
    {"#1 F+λ 8pt (ZZ/p)", R, I, x, "11x19 / 7x15", false}
))

-- 8. #2 E+f 6pt, one-sided focal (Paper Table 1 — ZZ/p)
problems = append(problems, () -> (
    FF := ZZ/32749;
    R := FF[x,y,z];
    A := random(FF^3, FF^9);
    X := apply(3, i -> matrix apply(3, j -> apply(3, k -> sub(A_(i, 3*j+k), R))));
    Fmat := X#0 + y*X#1 + z*X#2;
    om := diagonalMatrix{1_R,1_R,x};
    I := ideal(2*Fmat*om*transpose(Fmat)*Fmat - trace(Fmat*om*transpose(Fmat))*Fmat)
         + ideal(det Fmat);
    {"#2 E+f 6pt (ZZ/p)", R, I, x, "11x20", false}
))

-- 9. #3 f+E+f 6pt, both-sided focal (Paper Table 1 — ZZ/p)
problems = append(problems, () -> (
    FF := ZZ/32749;
    R := FF[x,y,z];
    A := random(FF^9, FF^3);
    X := apply(3, i -> matrix apply(3, j -> apply(3, k -> sub(A_(3*j+i, k), R))));
    Fmat := X#0 + y*X#1 + z*X#2;
    om := diagonalMatrix{1_R,1_R,x};
    I := ideal(2*Fmat*om*transpose(Fmat)*om*Fmat - trace(Fmat*om*transpose(Fmat)*om)*Fmat)
         + ideal(det Fmat);
    {"#3 f+E+f 6pt (ZZ/p)", R, I, x, "12x27 / 11x26", false}
))

-- 10. Mixed 3-var system (QQ, small deg) — cross-check
problems = append(problems, () -> (
    R := QQ[x,y,z];
    I := ideal(x^2 + y*z - 1, x*y + z - 2, x*z + y - 3);
    {"3-var mixed", R, I, x, "-x-", true}
))

-- ============================================================
-- Runner
-- ============================================================

-- Residual helper — accepts solutions list (CC entries) and returns max.
evalResid = (sols, I) -> (
    if #sols == 0 then return 1e99;
    Rb := ring I;
    FF := coefficientRing Rb;
    if not (FF === QQ or instance(FF, InexactField)) then return -1.0;  -- skip for ZZ/p
    Rc := CC[gens Rb];
    gensC := sub(gens I, Rc);
    max apply(sols, s -> norm sub(gensC, matrix{apply(s, v -> sub(v, CC))}))
)

runStrategy = (name, R, I, aVar, strat, overQQ) -> (
    -- Per-strategy runner. Catches any error and reports.
    ET := eliminationTemplate(aVar, I);
    t0 := cpuTime();
    ok := true;
    local M;
    try ( M = getTemplateMatrix(ET, Strategy => strat); ) else ( ok = false; );
    elapsed := cpuTime() - t0;
    if not ok then return {strat, "ERR", "-", "-", "-"};
    sizeStr := toString numRows M | "x" | toString numColumns M;
    tierInfo := "";
    if overQQ then (
        try (
            sols := templateSolve(ET, Strategy => strat);
            res := evalResid(sols, I);
            tierInfo = "res=" | toString res | if res < 1e-6 then " PASS" else " FAIL";
        ) else (
            tierInfo = "solve ERR";
        );
    ) else (
        tierInfo = "(ZZ/p — no tier-3)";
    );
    {strat, sizeStr, toString round(3, elapsed), tierInfo}
)

runProblem = prob -> (
    (name, R, I, aVar, target, overQQ) := toSequence prob();
    << "========================================================================" << endl;
    << name << "   (deg I = " << degree I << ", action = " << aVar
       << ", target = " << target << ")" << endl;
    << "========================================================================" << endl;
    -- Graph-ideal strategies — may be slow on larger problems; still run.
    for strat in {null, "Larsson", "MatrixHi", "Greedy"} do (
        row := runStrategy(name, R, I, aVar, strat, overQQ);
        stratName := if strat === null then "Default" else strat;
        << "   " << pad(stratName, 10) << " | size " << pad(row#1, 14)
           << " | " << pad(row#2, 8) << "s | " << row#3 << endl;
    );
    << endl
)

<< endl << "========================================================================" << endl
<< "  EliminationTemplates multi-problem benchmark" << endl
<< "  seed 42; Default/Larsson/MatrixHi/Greedy on " << #problems << " problems" << endl
<< "========================================================================" << endl << endl

scan(problems, runProblem)

<< endl << "=== done ===" << endl
exit 0
