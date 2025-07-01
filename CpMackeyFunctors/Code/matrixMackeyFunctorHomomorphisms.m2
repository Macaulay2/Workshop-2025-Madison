-- needsPackage Complexes???
isAbelianCategory CpMackeyFunctor := M -> true
isMorphism MackeyFunctorHomomorphism := f -> true

MackeyFunctorHomomorphism | MackeyFunctorHomomorphism := MackeyFunctorHomomorphism => MackeyFunctorHomomorphism.concatCols = maps -> (
    if not all(maps, f -> target f === target maps#0) or not all(maps, f -> getPrimeOrder source f === getPrimeOrder source maps#0) then
        error "MackeyFunctorHomomorphism.concatCols: all maps must have the same target and prime order";
    if #maps === 0 then
        error "MackeyFunctorHomomorphism.concatCols: no maps provided";
    map(target maps#0, directSum apply(maps, source), concatCols apply(maps, getUnderlyingMap), concatCols apply(maps, getFixedMap))
)

MackeyFunctorHomomorphism || MackeyFunctorHomomorphism := MackeyFunctorHomomorphism => MackeyFunctorHomomorphism.concatRows = maps -> (
    if not all(maps, f -> source f === source maps#0) or not all(maps, f -> getPrimeOrder source f === getPrimeOrder source maps#0) then
        error "MackeyFunctorHomomorphism.concatRows: all maps must have the same source and prime order";
    if #maps === 0 then
        error "MackeyFunctorHomomorphism.concatRows: no maps provided";
    map(directSum apply(maps, target), source maps#0, concatRows apply(maps, getUnderlyingMap), concatRows apply(maps, getFixedMap))
)

MackeyFunctorHomomorphism.concatBlocks = maps -> MackeyFunctorHomomorphism.concatRows apply(maps, MackeyFunctorHomomorphism.concatCols)
MackeyFunctorHomomorphism.matrix = opts -> MackeyFunctorHomomorphism.concatBlocks