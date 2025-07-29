needsPackage "EliminationTemplates"

TestCase = new Type of HashTable;
testCase = method();
testCase (String, Ideal, ZZ, ZZ) := (name, I, n, m) -> (
    new TestCase from {
        "name" => name,
        ideal => I,
        "dims" => n | " x " | m
    }
)

runBenchmarks = method();
runBenchmarks = () -> (
    testCases = loadBenchmarks({});
    testCases = loadKukelovaBenchmarks(testCases);
    results = {};

    for testCase in testCases do (
        I = testCase#ideal;
        template = eliminationTemplate(random(1, ring I), I);
        M = getTemplateMatrix(template);
        results = results | {{testCase#"name", toString(numRows M) | " x " | toString(numColumns M), testCase#"dims"}};
    );

    -- print the table
    nameLen = max(apply(results, r -> #toString(r#0)));
    myDimLen = #"Template Dim";
    litDimLen = #"Literature Dim";
    header = "| " | pad("Problem", nameLen) | " | " | "Template Dim" | " | " | "Literature Dim" | " |";
    separator = concatenate((nameLen + myDimLen + litDimLen + 10):"-");

    print separator;
    print header;
    print separator;

    for result in results do (
        row = "| " | pad(result#0, nameLen) | " | " | pad(result#1, myDimLen) | " | " | pad(result#2, litDimLen) | " |";
        print row;
    );
    print separator;
)

-- load the 3 benchmarks from:
-- https://openaccess.thecvf.com/content_cvpr_2017/papers/Kukelova_A_Clever_Elimination_CVPR_2017_paper.pdf
loadKukelovaBenchmarks = method();
loadKukelovaBenchmarks (List) := (testCases) -> (
    R = QQ[w,x,y];
    Fs = apply(3, i -> random(QQ^3, QQ^3));
    F = x * Fs#0 + y * Fs#1 + Fs#2;
    Q = diagonalMatrix({1, 1, w});

    -- f+E+f relative pose
    I = ideal(F * Q * transpose F * Q * F - (1/2) * trace(F * Q * transpose F * Q) * F) + ideal(det F);
    testCases = testCases | {testCase("Rel. pose + const. focal 6pt", I, 31, 46)};

    -- E+f 6pt relative pose
    I = ideal(F * Q * transpose F * F - (1/2) * trace(F * Q * transpose F) * F) + ideal(det F);
    testCases = testCases | {testCase("Rel. pose + one focal 6pt", I, 21, 30)};
    
    -- E+f+k 7pt relative pose
    R = QQ[w,x,y,lambda];
    mons = {x^2, y^2, lambda^2, x*y, x*lambda, y*lambda};
    coeffs = apply(6, i -> random(QQ));
    h = sum(0..#mons-1, i -> coeffs#i * mons#i);  -- random quadratic function
    Fs = apply(4, i -> random(QQ^3, QQ^3));
    F = x * Fs#0 + y * Fs#1 + lambda * Fs#2 + Fs#3;
    Q = sub(Q, R);
    I = ideal(F * Q * transpose F * Q * F - (1/2) * trace(F * Q * transpose F * Q) * F) + ideal(det F) + ideal(lambda * y - h);
    testCases | {testCase("Rel. pose 7pt one-sided focal + rad. dist.", I, 185, 204)}
)

loadBenchmarks = method();
loadBenchmarks (List) := (testCases) -> (
    R = QQ[x,y,z];
    Es = apply(4, i -> random(QQ^3, QQ^3));
    E = x * Es#0 + y * Es#1 + z * Es#2 + Es#3;
    I = ideal(E * transpose E * E - (1/2) * trace(E * transpose E) * E);
    testCases | {testCase("Rel. Pose 5pt", I, 10, 20)}

-*
    Tests from https://openaccess.thecvf.com/content_cvpr_2017/papers/Larsson_Efficient_Solvers_for_CVPR_2017_paper.pdf
    -- Rel. pose + rad. dist. 8pt Kukelova et al. [33] (*) 32 × 48
    -- Rel. pose + rad. dist. 6pt Kukelova et al. [33] (*) 238 × 290 
    -- Rel. pose + 2 rad. dist. 9pt Kukelova et al. [33] (*) 179 × 203
    Rel. pose 8pt one-sided rad. dist. Kuang et al. [30] 12 × 24 11 × 20 11 × 20
    TDOA offset rank 2, 7,4 pts Kuang et al. [28] 20 × 15 20 × 15 20 × 15
    P3.5P + focal Wu [52] 20 × 43 24 × 45 20 × 44
    Rel. pose 6pt ones-sided rad. dist. Kuang et al. [30] 48 × 70 34 × 60 34 × 60
    TDOA offset rank 2, 5,6 pts Kuang et al. [28] 105 × 83 105 × 83 40 × 42
    Rolling shutter pose Saurer et al. [44] (*) 48 × 56 50 × 55 47 × 55
    Generalized P4P + scale Ventura et al. [51] (*) 48 × 56 50 × 55 47 × 55
    Stitching + const. focal + rad. dist. 3pt Naroditsky et al. [39] 54 × 77 96 × 108 48 × 66
    TDOA offset rank 3, 9,5 pts Kuang et al. [28] 70 × 31 70 × 31 70 × 31
    TDOA offset rank 3, 7,6 pts Kuang et al. [28] 255 × 157 255 × 157 75 × 57
    Generalized rel. pose 6pt Stewenius ´ et al. [48] 60 × 120‡ 135 × 164 99 × 163
    Optimal PnP Hesch et al. [21] 120 × 120 93 × 116 88 × 115
    Triangulation from satellite im. Zheng et al. [53] (*) 93 × 120 93 × 116 88 × 115
    Optimal PnP (Cayley) Nakano [38] (*) 124 × 164 186 × 161 118 × 158
    P4P + focal + rad. dist. Bujnak et al. [5] (*) 136 × 152 140 × 144 140 × 156
    Weak PnP Larsson et al. [35] 234 × 276 568 × 498† 189 × 232 †
    Weak PnP (2x2 sym) Larsson et al. [35] 104 × 90 83 × 90† 49 × 59†
    Rolling shutter R6P Albl et al. [2] (*) 196 × 216 222 × 230 204 × 224
    Optimal pose w dir 4pt Svarm¨ et al. [49] 280 × 252 371 × 351 203 × 239
    Rel. pose w dir. 3pt Saurer et al. [45] (*) 411 × 489 287 × 324 210 × 255
    Rel. pose w dir. 3pt (using sym.) - - 94 × 111† 40 × 57†
    Abs. pose quivers Kuang et al. [27] 372 × 386 420 × 406 217 × 253
    Rel. pose w angle 4pt Li et al. [36] (*) 270 × 290 280 × 304 266 × 329
    Refractive P5P Haner et al. [16] 280 × 399 410 × 480 240 × 324
    TDOA offset rank 3, 6,8 pts Kuang et al. [28] 1359 × 754 1359 × 754 356 × 345
    Optimal PnP Zheng et al. [54] (*) 575 × 656 812 × 704 521 × 601
    Optimal PnP (using sym.) Zheng et al. [54] (*) 348 × 376 484 × 408† 302 × 342†
    Optimal pose w dir 3pt Svarm¨ et al. [49] 1, 260 × 1, 278 918 × 726 544 × 592
    Optimal PnP (quaternion) Nakano [38] (*) 630 × 710 958 × 693 604 × 684
    Refractive P6P + focal Haner et al. [16] 648 × 917 2, 196 × 1, 913† 636 × 851†
    Rel. pose + const. focal + rad. dist. 7pt Jiang et al. [23] 886 × 1, 011 1, 393 × 1, 237 581 × 862
    Dual-Receiver TDOA 5pt Burgess et al. [7] 2, 625 × 2, 352 850 × 1, 167 455 × 768
    Optimal PnP (rot. matrix) Nakano [38] (*) 1, 936 × 1, 976 1, 698 × 1, 153 1, 102 × 1, 135
    L2 3 view triangulation (Relaxed) Kukelova et al. [34] (*) 274 × 305 399 × 384 239 × 290
    L2 3 view triangulation Kukelova et al. [34] (*) 1, 866 × 1, 975 2, 647 × 2, 584 1, 759 × 2, 013
    *-
)