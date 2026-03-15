-- Helper routines for Greedy strategy in EliminationTemplates
--
-- This file keeps a single affine assignment engine and follows the
-- excessive-monomial prioritization idea used in greedyAG.mw (adjustParams).

thetaKey = (i, j) -> toString(i) | "|" | toString(j);

isZeroScalar = (v) -> v == 0_(ring v);

debugEnabled = (options) -> options#?"debug" and options#"debug";

debugPrint = (options, msg) -> (
    if debugEnabled(options) then print("Greedy debug: " | msg);
);

timedDebug = (options, label, thunk) -> (
    if debugEnabled(options) then (
        print("Greedy debug: start " | label);
        t := timing thunk();
        print("Greedy debug: end " | label | " (" | toString(t#0) | " s)");
        t#1
    ) else thunk()
);

copyMutableHash = (H) -> (
    C := new MutableHashTable from {};
    scan(keys H, k -> C#k = H#k);
    C
);

beginAffineTransaction = () -> new MutableHashTable from {
    "log" => new MutableHashTable from {},
    "size" => 0
};

recordAffineTxn = (txn, op) -> (
    if txn =!= null then (
        k := txn#"size";
        (txn#"log")#k = op;
        txn#"size" = k + 1;
    );
);

rollbackAffineTransaction = (txn) -> (
    if txn === null then return null;
    log := txn#"log";
    size := txn#"size";
    if size > 0 then (
        scan(reverse(0..size-1), idx -> (
            op := log#idx;
            tag := op#0;
            if tag == "assignAdd" then (
                remove(op#1, op#2);
            ) else if tag == "eqA" then (
                (op#1)#"a" = op#2;
            ) else if tag == "eqB" then (
                (op#1)#(op#2) = op#3;
            ) else if tag == "zeroCount" then (
                (op#1)#(op#2) = op#3;
            ) else if tag == "isZeroCol" then (
                (op#1)#(op#2) = op#3;
            ) else if tag == "totalZeroCols" then (
                (op#1)#"totalZeroCols" = op#2;
            );
        ));
    );
    null
);

polyToScalarMap = (f, baseR) -> (
    K := coefficientRing baseR;
    M := new MutableHashTable from {};
    if f =!= 0_baseR then (
        (monsMat, coeffsMat) := coefficients f;
        mons := flatten entries monsMat;
        coeffs := flatten entries coeffsMat;
        scan(0..#mons-1, q -> (
            c := coeffs#q;
            if c =!= 0_K then (
                if M#?(mons#q) then M#(mons#q) = M#(mons#q) + c else M#(mons#q) = c
            );
        ));
    );
    M
);

buildRowWiseAffineStateWithOptions = (H0e, H1, baseR, options) -> (
    d := numRows H0e;
    nGen := numColumns H0e;
    l := numRows H1;
    K := coefficientRing baseR;
    monsSeen := new MutableHashTable from {};
    emitDebug := options =!= null and options#?"debug" and options#"debug";

    if emitDebug then print("Greedy debug: buildRowWiseAffineState dims H0e=" | toString(d) | "x" | toString(nGen) | ", H1=" | toString(l) | "x" | toString(nGen));

    mapsTimed := timing (
        h0Maps := apply(d, i -> apply(nGen, k -> (
            map0 := polyToScalarMap(H0e_(i,k), baseR);
            scan(keys map0, m -> monsSeen#m = true);
            map0
        )));
        h1Maps := apply(l, j -> apply(nGen, k -> (
            map1 := polyToScalarMap(H1_(j,k), baseR);
            scan(keys map1, m -> monsSeen#m = true);
            map1
        )));
    );
    if emitDebug then print("Greedy debug: buildRowWiseAffineState maps built in " | toString(mapsTimed#0) | " s");

    mons := toList set(keys monsSeen);
    if #mons == 0 then (
        return new MutableHashTable from {
            "cols" => new MutableHashTable from {},
            "incidence" => new MutableHashTable from {},
            "assignments" => new MutableHashTable from {},
            "zeroCount" => new MutableHashTable from {},
            "isZeroCol" => new MutableHashTable from {},
            "totalZeroCols" => 0,
            "d" => d,
            "l" => l,
            "nGen" => nGen,
            "numWCols" => 0,
            "mons" => {},
            "baseR" => baseR
        }
    );

    monIndex := new MutableHashTable from {};
    scan(0..#mons-1, q -> monIndex#(mons#q) = q);

    numWCols := nGen * #mons;
    cols := new MutableHashTable from {};
    incidence := new MutableHashTable from {};

    ensureEq := (c, i) -> (
        rowMap := if cols#?c then cols#c else (
            t := new MutableHashTable from {};
            cols#c = t;
            t
        );
        if rowMap#?i then rowMap#i else (
            eq := new MutableHashTable from {
                "a" => 0_K,
                "b" => new MutableHashTable from {}
            };
            rowMap#i = eq;
            eq
        )
    );

    assembleTimed := timing (
        scan(0..d-1, i -> (
            scan(0..nGen-1, k -> (
                map0 := (h0Maps#i)#k;
                scan(keys map0, m -> (
                    c := k * #mons + monIndex#m;
                    eq := ensureEq(c, i);
                    eq#"a" = eq#"a" + map0#m;
                ));

                scan(0..l-1, j -> (
                    map1 := (h1Maps#j)#k;
                    scan(keys map1, m -> (
                        coeff := map1#m;
                        if coeff =!= 0_K then (
                            c := k * #mons + monIndex#m;
                            eq := ensureEq(c, i);
                            bMap := eq#"b";
                            if bMap#?j then bMap#j = bMap#j + coeff else bMap#j = coeff;
                            tKey := thetaKey(i, j);
                            if incidence#?tKey then incidence#tKey#c = true else (
                                inc := new MutableHashTable from {};
                                inc#c = true;
                                incidence#tKey = inc;
                            );
                        );
                    ));
                ));
            ));
        ));
    );
    if emitDebug then print("Greedy debug: buildRowWiseAffineState equations assembled in " | toString(assembleTimed#0) | " s");

    -- Objective tracks realized support of instantiated H:
    -- unassigned theta entries are treated as 0 in instantiateHRowWiseFromAssignments.
    eqIsZero := (eq) -> isZeroScalar(eq#"a");

    zeroCount := new MutableHashTable from {};
    isZeroCol := new MutableHashTable from {};
    totalZero := 0;
    zeroInitTimed := timing (
        scan(0..numWCols-1, c -> (
            zc := d;
            if cols#?c then (
                rowMap := cols#c;
                scan(keys rowMap, i -> (
                    if not eqIsZero(rowMap#i) then zc = zc - 1;
                ));
            );
            zeroCount#c = zc;
            isZeroCol#c = (zc == d);
            if zc == d then totalZero = totalZero + 1;
        ));
    );
    if emitDebug then print("Greedy debug: buildRowWiseAffineState zero-count initialized in " | toString(zeroInitTimed#0) | " s");
    if emitDebug then print("Greedy debug: buildRowWiseAffineState summary mons=" | toString(#mons) | ", numWCols=" | toString(numWCols) | ", zeroCols=" | toString(totalZero));

    new MutableHashTable from {
        "cols" => cols,
        "incidence" => incidence,
        "assignments" => new MutableHashTable from {},
        "zeroCount" => zeroCount,
        "isZeroCol" => isZeroCol,
        "totalZeroCols" => totalZero,
        "d" => d,
        "l" => l,
        "nGen" => nGen,
        "numWCols" => numWCols,
        "mons" => mons,
        "baseR" => baseR
    }
);

buildRowWiseAffineState = (H0e, H1, baseR) -> buildRowWiseAffineStateWithOptions(H0e, H1, baseR, null);

copyAffineState = (state) -> (
    colsOld := state#"cols";
    colsNew := new MutableHashTable from {};
    scan(keys colsOld, c -> (
        rowOld := colsOld#c;
        rowNew := new MutableHashTable from {};
        scan(keys rowOld, i -> (
            eqOld := rowOld#i;
            eqNew := new MutableHashTable from {
                "a" => eqOld#"a",
                "b" => copyMutableHash(eqOld#"b")
            };
            rowNew#i = eqNew;
        ));
        colsNew#c = rowNew;
    ));

    new MutableHashTable from {
        "cols" => colsNew,
        "incidence" => state#"incidence",
        "assignments" => copyMutableHash(state#"assignments"),
        "zeroCount" => copyMutableHash(state#"zeroCount"),
        "isZeroCol" => copyMutableHash(state#"isZeroCol"),
        "totalZeroCols" => state#"totalZeroCols",
        "d" => state#"d",
        "l" => state#"l",
        "nGen" => state#"nGen",
        "numWCols" => state#"numWCols",
        "mons" => state#"mons",
        "baseR" => state#"baseR"
    }
);

applyThetaAssignmentWithTxn = (state, i, j, val, txn) -> (
    cols := state#"cols";
    incidence := state#"incidence";
    assigned := state#"assignments";
    zeroCount := state#"zeroCount";
    isZeroCol := state#"isZeroCol";
    d := state#"d";

    eqIsZero := (eq) -> isZeroScalar(eq#"a");

    tKey := thetaKey(i, j);
    if assigned#?tKey then (
        if assigned#tKey == val then return true;
        return false;
    );
    recordAffineTxn(txn, {"assignAdd", assigned, tKey});
    assigned#tKey = val;

    if incidence#?tKey then (
        scan(keys(incidence#tKey), c -> (
            if cols#?c then (
                rowMap := cols#c;
                if rowMap#?i then (
                    eq := rowMap#i;
                    wasZero := eqIsZero(eq);
                    bMap := eq#"b";
                    if bMap#?j then (
                        coeff := bMap#j;
                        if coeff =!= 0_(ring coeff) then (
                            recordAffineTxn(txn, {"eqA", eq, eq#"a"});
                            recordAffineTxn(txn, {"eqB", bMap, j, coeff});
                            eq#"a" = eq#"a" + val * coeff;
                            bMap#j = 0_(ring coeff);
                        );
                    );
                    nowZero := eqIsZero(eq);
                    if wasZero and not nowZero then (
                        recordAffineTxn(txn, {"zeroCount", zeroCount, c, zeroCount#c});
                        zeroCount#c = zeroCount#c - 1;
                    ) else if (not wasZero) and nowZero then (
                        recordAffineTxn(txn, {"zeroCount", zeroCount, c, zeroCount#c});
                        zeroCount#c = zeroCount#c + 1;
                    );
                    wasColZero := isZeroCol#c;
                    nowColZero := (zeroCount#c == d);
                    if wasColZero =!= nowColZero then (
                        recordAffineTxn(txn, {"isZeroCol", isZeroCol, c, wasColZero});
                        isZeroCol#c = nowColZero;
                        recordAffineTxn(txn, {"totalZeroCols", state, state#"totalZeroCols"});
                        if wasColZero and not nowColZero then state#"totalZeroCols" = state#"totalZeroCols" - 1;
                        if (not wasColZero) and nowColZero then state#"totalZeroCols" = state#"totalZeroCols" + 1;
                    );
                );
            );
        ));
    );
    true
);

applyThetaAssignment = (state, i, j, val) -> applyThetaAssignmentWithTxn(state, i, j, val, null);

enforceZeroForColumnsAffineWithTxn = (state, colsToZero, txn) -> (
    cols := state#"cols";
    d := state#"d";
    isZeroCol := state#"isZeroCol";

    activeVars := (eq) -> select(keys(eq#"b"), j -> not isZeroScalar((eq#"b")#j));

    scan(colsToZero, c -> (
        scan(0..d-1, i -> (
            if cols#?c and (cols#c)#?i then (
                eq := (cols#c)#i;
                if not (isZeroScalar(eq#"a") and #activeVars(eq) == 0) then (
                    vars := activeVars(eq);
                    if #vars == 0 then return false;
                    pivot := first vars;
                    coeff := (eq#"b")#pivot;
                    if isZeroScalar(coeff) then return false;
                    val := - (eq#"a") / coeff;
                    if not applyThetaAssignmentWithTxn(state, i, pivot, val, txn) then return false;
                );
            );
        ));
        if not isZeroCol#c then return false;
    ));
    true
);

enforceZeroForColumnsAffine = (state, colsToZero) -> enforceZeroForColumnsAffineWithTxn(state, colsToZero, null);

columnFreeThetaCount = (state, c) -> (
    if not ((state#"cols")#?c) then 0 else (
        rowMap := (state#"cols")#c;
        total := 0;
        scan(keys rowMap, i -> (
            eq := rowMap#i;
            total = total + #(select(keys(eq#"b"), j -> not isZeroScalar((eq#"b")#j)));
        ));
        total
    )
);

columnUnresolvedCount = (state, c) -> (
    d := state#"d";
    zc := (state#"zeroCount")#c;
    d - zc
);

deterministicMix = (seed, x) -> (
    m := 2147483647;
    y := ((seed + 104729) * 1103515245 + (x + 1) * 12345 + 1013904223) % m;
    if y < 0 then y + m else y
);

columnScoreStats = (state, c, thetaImpactCache) -> (
    if not ((state#"cols")#?c) then (
        {0, 0, 0, 0}
    ) else (
        rowMap := (state#"cols")#c;
        freeTheta := 0;
        unresolved := columnUnresolvedCount(state, c);
        thetaSeen := new MutableHashTable from {};
        scan(keys rowMap, i -> (
            eq := rowMap#i;
            active := select(keys(eq#"b"), j -> not isZeroScalar((eq#"b")#j));
            freeTheta = freeTheta + #active;
            scan(active, j -> thetaSeen#(thetaKey(i, j)) = true);
        ));
        downstream := 0;
        scan(keys thetaSeen, tKey -> (
            if thetaImpactCache#?tKey then (
                downstream = downstream + thetaImpactCache#tKey;
            ) else (
                impact := if (state#"incidence")#?tKey then #select(keys((state#"incidence")#tKey), cc -> not ((state#"isZeroCol")#cc)) else 0;
                thetaImpactCache#tKey = impact;
                downstream = downstream + impact;
            );
        ));
        {freeTheta, unresolved, downstream, #keys thetaSeen}
    )
);

candidateShortlist = (state, options, roundNo) -> (
    shortMax := if options#?"shortlistMax" then options#"shortlistMax" else 200;
    shortMin := if options#?"shortlistMin" then options#"shortlistMin" else 20;
    restartSeed := if options#?"restartSeed" then options#"restartSeed" else 0;
    activeCols := select(0..((state#"numWCols")-1), c -> not ((state#"isZeroCol")#c));
    activeCount := #activeCols;
    if activeCount == 0 then return {};

    shortSize := min(shortMax, max(shortMin, (activeCount + 9) // 10));
    shortSize = min(shortSize, activeCount);

    thetaImpactCache := new MutableHashTable from {};
    scores := toList apply(activeCols, c -> (
        stats := columnScoreStats(state, c, thetaImpactCache);
        tie := deterministicMix(restartSeed + roundNo * 7919, c);
        {stats#0, stats#1, -(stats#2), -(stats#3), tie, c}
    ));
    sorted := sort scores;
    if shortSize == 0 then {} else apply(0..shortSize-1, idx -> (sorted#idx)#5)
);

evaluateRowTrial = (state, colsToZero, baseZero, options) -> (
    useRollback := if options#?"trialWithRollback" then options#"trialWithRollback" else true;
    if useRollback then (
        txn := beginAffineTransaction();
        enforceTimed := timing enforceZeroForColumnsAffineWithTxn(state, colsToZero, txn);
        ok := enforceTimed#1;
        gain := if ok then state#"totalZeroCols" - baseZero else -1;
        rollbackTimed := timing rollbackAffineTransaction(txn);
        new HashTable from {
            "ok" => ok,
            "gain" => gain,
            "copyTime" => 0,
            "enforceTime" => enforceTimed#0,
            "rollbackTime" => rollbackTimed#0,
            "trialTime" => enforceTimed#0 + rollbackTimed#0
        }
    ) else (
        trialTimed := timing copyAffineState(state);
        trial := trialTimed#1;
        enforceTimed2 := timing enforceZeroForColumnsAffine(trial, colsToZero);
        ok2 := enforceTimed2#1;
        gain2 := if ok2 then trial#"totalZeroCols" - baseZero else -1;
        new HashTable from {
            "ok" => ok2,
            "gain" => gain2,
            "copyTime" => trialTimed#0,
            "enforceTime" => enforceTimed2#0,
            "rollbackTime" => 0,
            "trialTime" => trialTimed#0 + enforceTimed2#0
        }
    )
);

runRowWiseShortlistGreedy = (affineState, options) -> (
    state := copyAffineState(affineState);
    restartId := if options#?"restartId" then options#"restartId" else 0;
    debugPrint(options, "row-wise init: restart=" | toString(restartId) | ", numWCols=" | toString(state#"numWCols") | ", zeroCols=" | toString(state#"totalZeroCols"));
    improved := true;
    rounds := 0;
    progressEvery := if options#?"progressEvery" then options#"progressEvery" else 10;
    lookaheadTopK := if options#?"lookaheadTopK" then options#"lookaheadTopK" else 0;
    lookaheadPairPool := if options#?"lookaheadPairPool" then options#"lookaheadPairPool" else lookaheadTopK;
    lookaheadOnZeroOnly := if options#?"lookaheadOnZeroOnly" then options#"lookaheadOnZeroOnly" else true;
    while improved do (
        rounds = rounds + 1;
        improved = false;
        baseZero := state#"totalZeroCols";
        bestGain := 0;
        bestCols := {};
        shortlist := candidateShortlist(state, options, rounds);
        debugPrint(options, "row-wise round " | toString(rounds) | ": shortlist=" | toString(#shortlist) | ", baseZero=" | toString(baseZero));
        feasibleCount := 0;
        copyTime := 0;
        enforceTime := 0;
        rollbackTime := 0;
        trialTime := 0;
        lookaheadFeasible := 0;
        lookaheadTrialTime := 0;
        if #shortlist > 0 then (
            for idx from 0 to #shortlist-1 do (
                c := shortlist#idx;
                trial := evaluateRowTrial(state, {c}, baseZero, options);
                copyTime = copyTime + trial#"copyTime";
                enforceTime = enforceTime + trial#"enforceTime";
                rollbackTime = rollbackTime + trial#"rollbackTime";
                trialTime = trialTime + trial#"trialTime";
                if trial#"ok" then (
                    feasibleCount = feasibleCount + 1;
                    gain := trial#"gain";
                    if gain > bestGain then (
                        bestGain = gain;
                        bestCols = {c};
                    );
                );
                if debugEnabled(options) and (((idx + 1) % progressEvery == 0) or (idx + 1 == #shortlist)) then (
                    print("Greedy debug: row-wise round " | toString(rounds) | " progress " | toString(idx + 1) | "/" | toString(#shortlist) | ", feasible=" | toString(feasibleCount) | ", bestGain=" | toString(bestGain));
                );
            );
        );

        if (#shortlist > 1) and (lookaheadTopK > 1) and ((not lookaheadOnZeroOnly) or (bestGain <= 0)) then (
            headCount := min(lookaheadTopK, #shortlist);
            pairPool := min(max(headCount, lookaheadPairPool), #shortlist);
            for h from 0 to headCount-1 do (
                c1 := shortlist#h;
                for t from 0 to pairPool-1 do (
                    c2 := shortlist#t;
                    if c2 =!= c1 then (
                        trial2 := evaluateRowTrial(state, {c1, c2}, baseZero, options);
                        copyTime = copyTime + trial2#"copyTime";
                        enforceTime = enforceTime + trial2#"enforceTime";
                        rollbackTime = rollbackTime + trial2#"rollbackTime";
                        trialTime = trialTime + trial2#"trialTime";
                        lookaheadTrialTime = lookaheadTrialTime + trial2#"trialTime";
                        if trial2#"ok" then (
                            lookaheadFeasible = lookaheadFeasible + 1;
                            gain2 := trial2#"gain";
                            if gain2 > bestGain then (
                                bestGain = gain2;
                                bestCols = {c1, c2};
                            );
                        );
                    );
                );
            );
        );

        debugPrint(options, "row-wise round " | toString(rounds) | " summary: feasible=" | toString(feasibleCount) | "/" | toString(#shortlist) | ", lookaheadFeasible=" | toString(lookaheadFeasible) | ", copyTime=" | toString(copyTime) | " s, enforceTime=" | toString(enforceTime) | " s, rollbackTime=" | toString(rollbackTime) | " s, trialTime=" | toString(trialTime) | " s, lookaheadTrialTime=" | toString(lookaheadTrialTime) | " s, bestGain=" | toString(bestGain));
        if bestGain > 0 and #bestCols > 0 then (
            commitTimed := timing enforceZeroForColumnsAffine(state, bestCols);
            if not commitTimed#1 then error "row-wise commit failed after successful trial";
            improved = true;
            debugPrint(options, "row-wise round " | toString(rounds) | " accepted cols=" | toString(bestCols) | "; zeroCols=" | toString(state#"totalZeroCols"));
        ) else (
            debugPrint(options, "row-wise round " | toString(rounds) | " no improvement");
        );
    );
    debugPrint(options, "row-wise done: rounds=" | toString(rounds) | ", zeroCols=" | toString(state#"totalZeroCols"));
    new HashTable from {
        "assignments" => state#"assignments",
        "zeroColumns" => state#"totalZeroCols",
        "rounds" => rounds,
        "state" => state
    }
);

instantiateHRowWiseFromAssignments = (H0e, H1, assignments, baseR) -> (
    d := numRows H0e;
    nGen := numColumns H0e;
    l := numRows H1;
    K := coefficientRing baseR;
    valFor := (i, j) -> (
        key := thetaKey(i, j);
        if assignments#?key then assignments#key else 0_K
    );
    matrix apply(d, i -> apply(nGen, k -> (
        s := H0e_(i,k);
        scan(0..l-1, j -> (
            coeff := valFor(i, j);
            if coeff =!= 0_K then s = s + sub(coeff, baseR) * H1_(j,k);
        ));
        s
    )))
);

buildAffineColumnInfo = (state) -> (
    mons := state#"mons";
    nGen := state#"nGen";
    flatten for k from 0 to nGen-1 list (
        apply(#mons, j -> new HashTable from {
            "generatorIndex" => k,
            "monomial" => mons#j
        })
    )
);

computeExcessiveMonomials = (a, B, J, columnInfo, baseR) -> (
    gensJ := toList J_*;
    allMons := set flatten apply(columnInfo, info -> (
        genIdx := info#"generatorIndex";
        m := info#"monomial";
        flatten entries monomials(m * sub(gensJ#genIdx, baseR))
    ));
    monsB := set flatten entries(lift(B, baseR));
    monsR := set flatten entries(a * lift(B, baseR)) - set flatten entries(lift(B, baseR));
    toList(allMons - union(monsR, monsB))
);

columnsForExcessiveMonomial = (e, columnInfo, J, baseR) -> (
    gensJ := toList J_*;
    select(0..#columnInfo - 1, c -> (
        info := columnInfo#c;
        genIdx := info#"generatorIndex";
        m := info#"monomial";
        coefficient(e, m * sub(gensJ#genIdx, baseR)) =!= 0_baseR
    ))
);

buildExcessiveToColumnsIndex = (excessiveMons, columnInfo, J, baseR, options) -> (
    mapE := new MutableHashTable from {};
    progressEvery := if options#?"indexProgressEvery" then options#"indexProgressEvery" else 100;
    if debugEnabled(options) then print("Greedy debug: index build start for " | toString(#excessiveMons) | " excessive monomials");
    if #excessiveMons > 0 then (
        for idx from 0 to #excessiveMons-1 do (
            e := excessiveMons#idx;
            mapE#e = columnsForExcessiveMonomial(e, columnInfo, J, baseR);
            if debugEnabled(options) and (((idx + 1) % progressEvery == 0) or (idx + 1 == #excessiveMons)) then (
                print("Greedy debug: index build progress " | toString(idx + 1) | "/" | toString(#excessiveMons));
            );
        )
    );
    mapE
);

shortlistExcessiveMons = (state, excessiveMons, getCols, options, roundNo) -> (
    shortMax := if options#?"shortlistMax" then options#"shortlistMax" else 50;
    shortMin := if options#?"shortlistMin" then options#"shortlistMin" else 5;
    probeFactor := if options#?"excessiveProbeFactor" then options#"excessiveProbeFactor" else 3;

    activeCount := #excessiveMons;
    if activeCount == 0 then return {};

    shortSize := min(shortMax, max(shortMin, (activeCount + 9) // 10));
    shortSize = min(shortSize, activeCount);
    probeSize := min(activeCount, max(shortSize, shortSize * probeFactor));

    probe := if activeCount <= probeSize then excessiveMons else (
        start := ((roundNo - 1) * probeSize) % activeCount;
        apply(0..probeSize-1, t -> excessiveMons#((start + t) % activeCount))
    );

    candidates := select(probe, e -> (
        colsE := getCols(e);
        any(colsE, c -> not ((state#"isZeroCol")#c))
    ));
    if #candidates == 0 then return {};

    -- Inspired by greedyAG adjustParams ordering: prioritize monomials with broader impact.
    scores := apply(candidates, e -> (
        colsE := getCols(e);
        activeImpact := #select(colsE, c -> not ((state#"isZeroCol")#c));
        {-activeImpact, #colsE, toString(e), e}
    ));
    sorted := sort scores;
    if shortSize == 0 then {} else apply(0..shortSize-1, idx -> (sorted#idx)#3)
);

runColumnWiseShortlistGreedyAffine = (affineState, excessiveMons, columnInfo, J, baseR, options) -> (
    state := copyAffineState(affineState);
    debugPrint(options, "column-wise init: excessiveMons=" | toString(#excessiveMons) | ", zeroCols=" | toString(state#"totalZeroCols"));
    precomputeAll := options#?"precomputeExcessiveIndex" and options#"precomputeExcessiveIndex";
    e2cols := if precomputeAll then timedDebug(options, "build excessive->columns index", () -> buildExcessiveToColumnsIndex(excessiveMons, columnInfo, J, baseR, options)) else new MutableHashTable from {};
    cacheHits := 0;
    cacheMisses := 0;
    getCols := e -> (
        if e2cols#?e then (
            cacheHits = cacheHits + 1;
            e2cols#e
        ) else (
            cacheMisses = cacheMisses + 1;
            cols := columnsForExcessiveMonomial(e, columnInfo, J, baseR);
            e2cols#e = cols;
            cols
        )
    );

    improved := true;
    rounds := 0;
    progressEvery := if options#?"progressEvery" then options#"progressEvery" else 10;
    while improved do (
        rounds = rounds + 1;
        improved = false;
        baseZero := state#"totalZeroCols";
        bestGain := 0;
        bestTrial := null;

        shortlist := shortlistExcessiveMons(state, excessiveMons, getCols, options, rounds);
        debugPrint(options, "column-wise round " | toString(rounds) | ": shortlist=" | toString(#shortlist) | ", baseZero=" | toString(baseZero));
        feasibleCount := 0;
        copyTime := 0;
        enforceTime := 0;
        trialTime := 0;
        if #shortlist > 0 then (
            for idx from 0 to #shortlist-1 do (
                e := shortlist#idx;
                colsE := getCols(e);
                if #colsE > 0 then (
                    trialTimed := timing copyAffineState(state);
                    trial := trialTimed#1;
                    copyTime = copyTime + trialTimed#0;
                    enforceTimed := timing enforceZeroForColumnsAffine(trial, colsE);
                    ok := enforceTimed#1;
                    enforceTime = enforceTime + enforceTimed#0;
                    trialTime = trialTime + trialTimed#0 + enforceTimed#0;
                    if ok then (
                        feasibleCount = feasibleCount + 1;
                        gain := trial#"totalZeroCols" - baseZero;
                        if gain > bestGain then (
                            bestGain = gain;
                            bestTrial = trial;
                        );
                    );
                );
                if debugEnabled(options) and (((idx + 1) % progressEvery == 0) or (idx + 1 == #shortlist)) then (
                    print("Greedy debug: column-wise round " | toString(rounds) | " progress " | toString(idx + 1) | "/" | toString(#shortlist) | ", feasible=" | toString(feasibleCount) | ", bestGain=" | toString(bestGain));
                );
            );
        );
        debugPrint(options, "column-wise round " | toString(rounds) | " summary: feasible=" | toString(feasibleCount) | "/" | toString(#shortlist) | ", copyTime=" | toString(copyTime) | " s, enforceTime=" | toString(enforceTime) | " s, trialTime=" | toString(trialTime) | " s, bestGain=" | toString(bestGain) | ", cacheSize=" | toString(#keys e2cols) | ", cacheHits=" | toString(cacheHits) | ", cacheMisses=" | toString(cacheMisses));

        if bestGain > 0 then (
            state = bestTrial;
            improved = true;
            debugPrint(options, "column-wise round " | toString(rounds) | " accepted; zeroCols=" | toString(state#"totalZeroCols"));
        ) else (
            debugPrint(options, "column-wise round " | toString(rounds) | " no improvement");
        );
    );
    debugPrint(options, "column-wise done: rounds=" | toString(rounds) | ", zeroCols=" | toString(state#"totalZeroCols"));

    new HashTable from {
        "assignments" => state#"assignments",
        "zeroColumns" => state#"totalZeroCols",
        "rounds" => rounds,
        "state" => state
    }
);
