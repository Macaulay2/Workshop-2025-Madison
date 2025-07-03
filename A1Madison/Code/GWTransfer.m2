transferGW = method()
transferGW (GrothendieckWittClass, Ring, Ring, Ideal) := GrothendieckWittClass => (alpha, kk, C, I) -> (
    -- Catch all the errors

    makeDiagonalForm(kk, apply(getDiagonalEntries(alpha), i -> AlgebraicTrace(C, I, i, kk)))
)