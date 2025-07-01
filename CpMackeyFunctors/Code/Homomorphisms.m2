protect symbol Domain
protect symbol Codomain

MackeyFunctorHomomorphism = new Type of HashTable
MackeyFunctorHomomorphism.synonym = "Mackey Functor homomorphism"

isWellDefined MackeyFunctorHomomorphism := Boolean => F ->(
    -- Verify the keys are correct
    if not (F#?Domain and F#?Codomain and F#?UnderlyingMap and F#?FixedMap) then return false;
    
    -- Verify domain and codomain are indeed Mackey functors
    if not (class F.Domain === CpMackeyFunctor and class F.Codomain === CpMackeyFunctor) then return false;

    -- Verify that the primes are the same
    if not (F.Domain.PrimeOrder == F.Codomain.PrimeOrder) then return false;

    -- Verify F.UnderlyingMap and F.FixedMap have the right domain and codomain
    if not (source F.UnderlyingMap == getUnderlyingModule F.Domain and target F.UnderlyingMap == getUnderlyingModule F.Codomain) then return false;
    if not (source F.FixedMap == getFixedModule F.Domain and target F.FixedMap == getFixedModule F.Codomain) then return false;

    -- Check commutes with restriction
    if not (F.UnderlyingMap * F.Domain.Res == F.Codomain.Res * F.FixedMap) then (print " -- the given morphism does not commute with restriction"; return false);

    -- Check commutes with transfer
    if not (F.FixedMap * F.Domain.Tr == F.Codomain.Tr * F.UnderlyingMap) then (print " -- the given morphism does not commute with transfer"; return false);

    -- Check commutes with conjugation
    if not (F.UnderlyingMap * F.Domain.Conj == F.Codomain.Conj * F.UnderlyingMap) then (print " -- the given morphism does not commute with conjugation"; return false);

    true
)

-- Arguments:
-- 1. Target
-- 2. Source
-- 3. Underlying map
-- 4. Fixed-point map
map(CpMackeyFunctor, CpMackeyFunctor, Matrix, Matrix) := MackeyFunctorHomomorphism => opts -> (N,M,u,f) -> (
    F := new MackeyFunctorHomomorphism from {
        symbol Codomain => N,
        symbol Domain => M,
        symbol UnderlyingMap => map(getUnderlyingModule(N),getUnderlyingModule(M),u),
        symbol FixedMap =>  map(getFixedModule(N),getFixedModule(M),f),
        symbol cache => new CacheTable
        };
    if isWellDefined F then (
        return F
    )
    else (
        error "Mackey Functor homomorphism is not well-defined";
	)
)

source(MackeyFunctorHomomorphism) := CpMackeyFunctor => F -> (
    return F.Domain
)

target(MackeyFunctorHomomorphism) := CpMackeyFunctor => F -> (
    return F.Codomain
)

CpMackeyFunctor#id = X -> map(X, X, id_(X.Underlying), id_(X.Fixed))

-- This is the linearization map A -> RU
complexLinearizationMap = method()
complexLinearizationMap(ZZ) := MackeyFunctorHomomorphism => p -> (
    map(makeComplexRepresentationMackeyFunctor p, makeBurnsideMackeyFunctor p, matrix {{1}}, matrix {{1,1}} || matrix (for i to p-2 list {0,1}))
)

-- This is the linearization map A -> RO
realLinearizationMap = method()
realLinearizationMap(ZZ) := MackeyFunctorHomomorphism => p -> (
    RO := makeRealRepresentationMackeyFunctor p;
    map(RO, makeBurnsideMackeyFunctor p, matrix {{1}}, matrix {{1,1}} || matrix (for i to (rank (getFixedModule RO) - 2) list {0,1}))
)

-- Given a Mackey functor M and vector x in fixed module, produce map A -> M
makeUniversalMapFixed = method()
makeUniversalMapFixed(CpMackeyFunctor,Vector) := MackeyFunctorHomomorphism => (M,x) -> (
    X := matrix x;
    if (target X) != M.Fixed then (
        error "element is not in fixed module";
    )
    else (
        p := M.PrimeOrder;
        A := makeFixedFreeMackeyFunctor(p);
        U := M.Res * X;
        F := X | (M.Tr * M.Res * X);
        return map(M, A, U, F);
    )
)

-- Given a Mackey functor M and vector x in underlying module, produce map B -> M
makeUniversalMapUnderlying = method()
makeUniversalMapUnderlying(CpMackeyFunctor,Vector) := MackeyFunctorHomomorphism => (M,x) -> (
    X := matrix x;

    if (target X) != M.Underlying then (
        error "element is not in underlying module";
    )
    else (
        p := M.PrimeOrder;
        B := makeUnderlyingFreeMackeyFunctor(p);
        U := matrix {for i to p-1 list ((M.Conj)^i) * X};
        F := M.Tr * X;
        return map(M, B, U, F);
    )
)


MackeyFunctorHomomorphism * MackeyFunctorHomomorphism := MackeyFunctorHomomorphism => (G,F) ->(
    -- todo: uncomment after Sasha pushes == method
    -- if not (F.Codomain == G.Domain) then error "Mackey functor maps are not composable";
    map(G.Codomain, F.Domain, G.UnderlyingMap * F.UnderlyingMap, G.FixedMap * F.FixedMap)
)