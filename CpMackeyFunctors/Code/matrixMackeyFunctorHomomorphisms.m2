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
