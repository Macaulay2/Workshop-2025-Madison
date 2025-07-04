transferGW = method()
transferGW (GrothendieckWittClass) := GrothendieckWittClass => (alpha) -> (

    kk := coefficientRing getAlgebra alpha;

    makeDiagonalForm(kk, toSequence apply(getDiagonalEntries(alpha), i -> algebraicTrace(getAlgebra alpha, i)))
)