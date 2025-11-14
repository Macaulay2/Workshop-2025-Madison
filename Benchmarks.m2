needsPackage "EliminationTemplates"

TestCase = new Type of HashTable;
testCase = method(Options => {"action polynomial" => null});
testCase (String, Ideal, ZZ, ZZ) := o -> (name, I, n, m) -> (
    R := ring I;
    a := if instance(o#"action polynomial", Nothing) then random(1, R) else (
        assert(instance(o#"action polynomial", R));
        o#"action polynomial"
	);
    new TestCase from {
        "name" => name,
        ideal => I,
        "dims" => n | " x " | m,
	    "action polynomial" => a
    }
)

runBenchmarks = method();
runBenchmarks = () -> (
    testCases := loadRelPoseBenchmarks({});
    --loadPnPBenchmarks(testCases);
    results := {};

    for testCase in testCases do (
        --print testCase#"name";
        I := testCase#ideal;
	    a := testCase#"action polynomial";
        template := eliminationTemplate(a, I);
        templateTest := timing getTemplateMatrix(template);
        templateTime = templateTest#0;
        M = templateTest#1;
        solveTime := (timing (templateSolve(template)))#0;

        results = results | {{testCase#"name", toString(numRows M) | " x " | toString(numColumns M), testCase#"dims", 
            toString(templateTime), toString(solveTime)}};
    );

    -- print the table
    nameLen := max(apply(results, r -> #toString(r#0)));
    myDimLen := #"Template Dim";
    litDimLen := #"Literature Dim";
    templateTimeLen := #"Template Time";
    solveTimeLen := #"Solve Time";
    header := "| " | pad("Problem", nameLen) | " | " | "Template Dim" | " | " | "Literature Dim" | " | ";
    header = header | "Template Time" | " | " | "Solve Time" | " |";
    separator := concatenate((nameLen + myDimLen + litDimLen + templateTimeLen + solveTimeLen + 16):"-");

    print separator;
    print header;
    print separator;

    for result in results do (
        row := "| " | pad(result#0, nameLen) | " | " | pad(result#1, myDimLen) | " | " | pad(result#2, litDimLen) | " | ";
        row = row | pad(result#3, templateTimeLen) | " | " | pad(result#4, solveTimeLen) | " |";
        print row;
    );
    print separator;
)

-- Relative pose benchmarks from:
-- https://openaccess.thecvf.com/content_cvpr_2017/papers/Kukelova_A_Clever_Elimination_CVPR_2017_paper.pdf
loadRelPoseBenchmarks = method();
loadRelPoseBenchmarks (List) := (testCases) -> (
    R := QQ[x,y,z];

    -- 5 pt relative pose
    Es := apply(4, i -> random(QQ^3, QQ^3));
    E := x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
    I := ideal(E * transpose E * E - (1/2) * trace(E * transpose E) * E);
    testCases = testCases | {testCase("Rel. Pose 5pt (random linear form)", I, 10, 20), 
        testCase("Rel. Pose 5pt (x variable)", I, 10, 20, "action polynomial" => x)};

    -- f+E+f 6pt relative pose
    R = QQ[w,x,y];
    Fs := apply(3, i -> random(QQ^3, QQ^3));
    F := x * Fs#0 + y * Fs#1 + Fs#2;
    Q := diagonalMatrix({1, 1, w});
    I = ideal(F * Q * transpose F * Q * F - (1/2) * trace(F * Q * transpose F * Q) * F) + ideal(det F);
    testCases = testCases | {testCase("Rel. pose + const. focal 6pt", I, 31, 46)};

    -- E+f 6pt relative pose
    I = ideal(F * Q * transpose F * F - (1/2) * trace(F * Q * transpose F) * F) + ideal(det F);
    testCases = testCases | {testCase("Rel. pose + one focal 6pt", I, 21, 30)};
    
    -- E+f+k 7pt relative pose
    R = QQ[w,x,y,lambda];
    mons := {x^2, y^2, lambda^2, x*y, x*lambda, y*lambda};
    coeffs := apply(6, i -> random(QQ));
    h := sum(0..#mons-1, i -> coeffs#i * mons#i);  -- random quadratic function
    Fs = apply(4, i -> random(QQ^3, QQ^3));
    F = x * Fs#0 + y * Fs#1 + lambda * Fs#2 + Fs#3;
    Q = sub(Q, R);
    I = ideal(F * Q * transpose F * Q * F - (1/2) * trace(F * Q * transpose F * Q) * F) + ideal(det F) + ideal(lambda * y - h);
    testCases | {testCase("Rel. pose 7pt one-sided focal + rad. dist.", I, 185, 204)}
)

-- Perspective n-point benchmarks from:
-- https://www.bmva-archive.org.uk/bmvc/2015/papers/paper078/paper078.pdf
loadPnPBenchmarks = method();
loadPnPBenchmarks (List) := (testCases) -> (
    R := QQ[x_1..x_9];
    M := random(QQ^9, QQ^9);
    Fs := apply(9, i->random(QQ^3, QQ^3));
    f = vector(flatten entries F);

    -- Optimal PnP (Cayley) 124 × 164 
    -- Optimal PnP (quaternion) 630 × 710
    -- Optimal PnP (rot. matrix) 1936 × 1976
)