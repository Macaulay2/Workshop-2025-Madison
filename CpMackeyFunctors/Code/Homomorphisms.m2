protect symbol Domain
protect symbol Codomain
protect symbol UnderlyingMap
protect symbol FixedMap

MackeyFunctorHomomorphism = new Type of HashTable
MackeyFunctorHomomorphism.synonym = "Mackey Functor homomorphism"

isWellDefinedCpMackeyFunctorHomomorphism = method()
isWellDefinedCpMackeyFunctorHomomorphism MackeyFunctorHomomorphism := Boolean => F ->(
    -- Verify the keys are correct
    if not (F#?Domain and F#?Codomain and F#?UnderlyingMap and F#?FixedMap) then return false;
    
    -- Verify domain and codomain are indeed Mackey functors
    if not (class F.Domain == CpMackeyFunctor and class F.Codomain == CpMackeyFunctor) then return false;

    -- Verify F.UnderlyingMap and F.FixedMap have the right domain and codomain
    if not (source F.UnderlyingMap == getUnderlyingModule F.Domain and target F.UnderlyingMap == getUnderlyingModule F.Codomain) then return false;
    if not (source F.FixedMap == getFixedModule F.Domain and target F.FixedMap == getFixedModule F.Codomain) then return false;

    -- Check commutes with restriction

    -- Check commutes with transfer

    -- Check commutes with conjugation
)

makeMackeyFunctorHomomorphism = method()
makeMackeyFunctorHomomorphism(CpMackeyFunctor, CpMackeyFunctor, Matrix, Matrix) := MackeyFunctorHomomorphism => (N,M,f,g) -> (
    F := new MackeyFunctorHomomorphism from {
        symbol Codomain => N,
        symbol Domain => M,
        symbol UnderlyingMap => f,
        symbol FixedMap => g,
        symbol cache => new CacheTable
        };
    if isWellDefinedCpMackeyFunctorHomomorphism F then (
        return F
    )
    else (
        error "Mackey Functor homomorphism is not well-defined";
	)
)