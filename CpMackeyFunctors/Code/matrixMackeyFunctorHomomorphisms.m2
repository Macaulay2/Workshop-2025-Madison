
-- needsPackage Complexes???
isAbelianCategory CpMackeyFunctor := M -> true
isMorphism MackeyFunctorHomomorphism := f -> true

-- input: a matrix of mackey homs {fij: Mj -> Ni}

matrixMackeyFunctorHomomorphism = method()
matrixMackeyFunctorHomomorphism(List) := MackeyFunctorHomomorphism => Fij -> (
    -- check dimensions of list
    -- TODO

    m := length Fij--first dimension
    n := length Fij_0--second dimension

    -- check sources and targets of fij
    -- TODO

    -- lists of sources and targets

    listOfSources := for i to m-1 list (source Fij_0_i);
    listOfTargets := for j to n-1 list (source Fij_j_0);
    listOfUnderlyingMaps := for i to m-1 list (for j to n-1 list(getUnderlyingMap(Fij_i_j)));
    listOfFixedMaps := for i to m-1 list (for j to n-1 list(getFixedMap(Fij_i_j)));


    sourceMatrixHomomorphism := directSum listOfSources;
    targetMatrixHomomorphism := directSum listOfTargets;


)

-- SheafMap |  SheafMap := SheafMap => SheafMap.concatCols = maps -> map(
--     target maps#0, directSum apply(maps, source), concatCols apply(maps, matrix))
MackeyFunctorHomomorphism | MackeyFunctorHomomorphism := MackeyFunctorHomomorphism => MackeyFunctorHomomorphism.concatCols = maps -> (
    --blah
)

-- SheafMap || SheafMap := SheafMap => SheafMap.concatRows = maps -> map(
--     directSum apply(maps, target), source maps#0, concatRows autotruncate maps)
MackeyFunctorHomomorphism || MackeyFunctorHomomorphism := MackeyFunctorHomomorphism => MackeyFunctorHomomorphism.concatRows = maps -> (
    --blah
)

MackeyFunctorHomomorphism.concatBlocks = maps -> MackeyFunctorHomomorphism.concatRows apply(maps, MackeyFunctorHomomorphism.concatCols)
MackeyFunctorHomomorphism.matrix = opts -> MackeyFunctorHomomorphism.concatBlocks