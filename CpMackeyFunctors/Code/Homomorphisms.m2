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

makeMackeyFunctorHomomorphism = method()
makeMackeyFunctorHomomorphism(CpMackeyFunctor, CpMackeyFunctor, Matrix, Matrix) := MackeyFunctorHomomorphism => (N,M,u,f) -> (
    F := new MackeyFunctorHomomorphism from {
        symbol Codomain => N,
        symbol Domain => M,
        symbol UnderlyingMap => u,
        symbol FixedMap => f,
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
