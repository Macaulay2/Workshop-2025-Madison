-- Theta-adjustment driver for the Greedy strategy.
--
-- This file isolates the "find theta assignments" logic into one place so
-- getH0 can call a single entrypoint.
--
-- Design intent (mapped to greedyAG.mw):
-- 1) Start from H = H0 + Theta * H1, where H1 is a syzygy direction matrix.
-- 2) Solve theta assignments greedily to maximize zero columns in the affine
--    template-support state.
-- 3) Instantiate H with the discovered assignments and return the adjusted H0.
--
-- In this M2 package, theta search is implemented by affine state
-- utilities in GreedyHelpers.m2:
--   buildRowWiseAffineState
--   runRowWiseShortlistGreedy
--   instantiateHRowWiseFromAssignments

defaultGreedyThetaOptions = () -> new MutableHashTable from {
    "shortlistMax" => 200,
    "shortlistMin" => 20,
    "rowRestarts" => 4,
    "restartSeed" => 12345,
    "restartSeedStep" => 7919,
    "lookaheadTopK" => 6,
    "lookaheadPairPool" => 10,
    "lookaheadOnZeroOnly" => true,
    "trialWithRollback" => true,
    "fallbackLarsson" => true,
    "debug" => false,
    "progressEvery" => 10
};

mergeGreedyThetaOptions = (userOptions) -> (
    base := defaultGreedyThetaOptions();
    if instance(userOptions, HashTable) then (
        scan(keys userOptions, k -> base#k = userOptions#k);
    );
    base
);

debugTheta = (opts, msg) -> (
    if opts#?"debug" and opts#"debug" then print("Greedy debug: " | msg);
);

timedTheta = (opts, label, thunk) -> (
    if opts#?"debug" and opts#"debug" then (
        print("Greedy debug: start " | label);
        t := timing thunk();
        print("Greedy debug: end " | label | " (" | toString(t#0) | " s)");
        t#1
    ) else thunk()
);

copyOpts = (opts) -> (
    out := new MutableHashTable from {};
    scan(keys opts, k -> out#k = opts#k);
    out
);

monomialSupportCardinality = (H) -> (
    #set flatten entries monomials H
);

runRowWiseWithRestarts = (affineState, opts) -> (
    restarts := if opts#?"rowRestarts" then max(1, opts#"rowRestarts") else 1;
    baseSeed := if opts#?"restartSeed" then opts#"restartSeed" else 12345;
    seedStep := if opts#?"restartSeedStep" then opts#"restartSeedStep" else 7919;
    bestRun := null;

    for r from 0 to restarts-1 do (
        runOpts := copyOpts(opts);
        runOpts#"restartId" = r;
        runOpts#"restartSeed" = baseSeed + r * seedStep;
        rowRun := timedTheta(opts, "row-wise shortlist greedy restart " | toString(r + 1) | "/" | toString(restarts), () -> runRowWiseShortlistGreedy(affineState, runOpts));
        debugTheta(opts, "restart " | toString(r + 1) | " result rounds=" | toString(rowRun#"rounds") | ", zeroCols=" | toString(rowRun#"zeroColumns"));

        if bestRun === null then (
            bestRun = rowRun;
        ) else if (rowRun#"zeroColumns" > bestRun#"zeroColumns") or ((rowRun#"zeroColumns" == bestRun#"zeroColumns") and (rowRun#"rounds" < bestRun#"rounds")) then (
            bestRun = rowRun;
        );
    );
    bestRun
);

-- Entry point called by getH0(..., Strategy => "Greedy").
--
-- Inputs
--   a : action polynomial
--   B : quotient-ring basis used by getH0/getTemplate
--   J : ideal
--   H0: baseline matrix returned by algebraic construction
--   H1: syzygy-direction matrix (transpose sub(syz(gens J), ring J))
--   userOptions: hash-table with tuning flags (optional)
--
-- Output
--   Adjusted H0 in the same orientation as getH0 expects.
greedyAGThetaAdjustH0 = (a, B, J, H0, H1, userOptions) -> (
    opts := mergeGreedyThetaOptions(userOptions);
    baseR := ring J;

    -- H0e is row-wise orientation used by affine helper routines.
    H0e := timedTheta(opts, "transpose H0", () -> transpose H0);

    -- Stage A: Build affine state from H0e and H1.
    -- This compiles every relevant equation into an affine form in theta vars.
    affineState := timedTheta(opts, "build affine state", () -> buildRowWiseAffineStateWithOptions(H0e, H1, baseR, opts));
    debugTheta(opts, "affine state summary numWCols=" | toString(affineState#"numWCols") | ", zeroCols=" | toString(affineState#"totalZeroCols"));

    -- Stage B: Row-wise greedy pass with restart seeds/tie-break variation.
    rowRun := timedTheta(opts, "row-wise restart sweep", () -> runRowWiseWithRestarts(affineState, opts));
    debugTheta(opts, "row-wise result rounds=" | toString(rowRun#"rounds") | ", zeroCols=" | toString(rowRun#"zeroColumns"));
    rowGain := rowRun#"zeroColumns" - affineState#"totalZeroCols";

    -- If no gain is found, keep baseline H0 exactly.
    if rowGain <= 0 then (
        if opts#?"fallbackLarsson" and opts#"fallbackLarsson" then (
            larssonH0 := timedTheta(opts, "Larsson fallback", () -> H0 % image(syz(gens(J))));
            baseScore := monomialSupportCardinality(H0);
            larssonScore := monomialSupportCardinality(larssonH0);
            debugTheta(opts, "row-wise gain is zero; fallback scores baseline=" | toString(baseScore) | ", Larsson=" | toString(larssonScore));
            if larssonScore < baseScore then return larssonH0;
        );
        debugTheta(opts, "row-wise gain is zero; returning baseline H0");
        return H0
    );

    -- Stage C: Instantiate adjusted H with row-wise assignments.
    HbestE := timedTheta(opts, "instantiate H from assignments", () -> instantiateHRowWiseFromAssignments(H0e, H1, rowRun#"assignments", baseR));
    transpose HbestE
);
