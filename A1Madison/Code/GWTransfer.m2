transferGW = method()
transferGW (GrothendieckWittClassOverAlgebra) := GrothendieckWittClass => (alpha) -> (
    -- Catch all the errors

    kk := coefficientRing getRing alpha;

    makeDiagonalForm(kk, toSequence apply(getDiagonalEntriesOverAlgebra(alpha), i -> AlgebraicTrace(getRing alpha, i)))
)