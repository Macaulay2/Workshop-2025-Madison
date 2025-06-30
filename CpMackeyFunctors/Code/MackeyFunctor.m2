protect symbol PrimeOrder
protect symbol Underlying
protect symbol Fixed
protect symbol Res
protect symbol Tr
protect symbol Conj

CpMackeyFunctor = new Type of HashTable
CpMackeyFunctor.synonym = "Cp Mackey Functor"

-- Check if a Mackey functor is a well-defined
isWellDefinedCpMackeyFunctor = method()
isWellDefinedCpMackeyFunctor CpMackeyFunctor := Boolean => M -> (
    ------------------------
    -- General type-checking
    ------------------------
    
    -- Ensure all the keys in the hash table defining a Mackey functor are indeed defined
    if not (M#?PrimeOrder and M#?Res and M#?Tr and M#?Conj and M#?Underlying and M#?Fixed) then (print "-- Hashtable does not have correct keys"; return false);

    -- Check that the input p is actually a prime number
    if not (class M.PrimeOrder === ZZ and isPrime(M.PrimeOrder)) then (print "-- p is not prime"; return false);
    
    -- Check fixed and underlying modules are Z-modules
    if not (isModule M.Fixed and isModule M.Underlying) then return false;
    if not (ring M.Fixed === ZZ and ring M.Underlying === ZZ) then (print "-- objects are not abelian groups"; return false);

    -- Check source and target of restriction are correct
    if not source(M.Res) == M.Fixed then (print " -- the source of res should be the fixed module"; return false);
    if not target(M.Res) == M.Underlying then (print " -- the target of res should be the underlying module"; return false);

    -- Check source and target of transfer are correct
    if not source(M.Tr) == M.Underlying then (print " -- the source of tr should be the underlying module"; return false);
    if not target(M.Tr) == M.Fixed then (print " -- the target of tr should be the fixed module"; return false);
    
    -- Check source and target of conjugation are correct
    if not source(M.Conj) == M.Underlying then (print " -- the source of conj should be the underlying module"; return false);
    if not target(M.Conj) == M.Underlying then (print " -- the target of conj should be the underlying module"; return false);

    ---------
    -- Axioms
    ---------

    -- Axiom 1: Conj is an automorphism of order dividing p
    if not isIsomorphism(M.Conj) then return false;
    if not (M.Conj)^(M.PrimeOrder) == id_(M.Underlying) then (print "-- Conj is not an automorphism of order dividing p"; return false);

    -- Axiom 2: res and tr are homomorphisms (item 3 in overleaf)
    if not isWellDefined M.Tr then (print " -- tr is not a homomorphism"; return false);
    if not isWellDefined M.Res then (print " -- res is not a homomorphism"; return false);

    -- Axiom 3: c*res = res and tr*c = tr
    if M.Conj * M.Res != M.Res then (print " -- c * res is not equal to res"; return false);
    if M.Tr * M.Conj != M.Tr then (print " -- tr * c is not equal to tr"; return false);

    -- Axiom 4: res * tr = sum of all conjugates (item 5 in overleaf)
    if not (M.Res * M.Tr == sum for i to M.PrimeOrder-1 list M.Conj^i) then (print " -- res * tr is not equal to the sum of all conjugates"; return false);

    true
)

makeCpMackeyFunctor = method()
-- Ordering for the input is:
-- 1. A prime number p
-- 2. Restriction matrix
-- 3. Transfer matrix
-- 4. Conjugation matrix
makeCpMackeyFunctor(ZZ,Matrix,Matrix,Matrix) := CpMackeyFunctor => (p,R,T,C) ->(
    M := new CpMackeyFunctor from {
        symbol PrimeOrder => p,
        symbol Underlying => source T,          -- extract the underlying module from the transfer homomorphism
        symbol Fixed => target T,               -- extract the fixed module from the transfer homomorphism
        symbol Res => R,
        symbol Tr => T,
        symbol Conj => C,
        symbol cache => new CacheTable
        };
    if isWellDefinedCpMackeyFunctor M then (
        return M
    )
    else (
        error "Mackey Functor is not well-defined";
	)
)

-- todo: verify that this still works when the first input is an Array


-------------------------
-- Recovering cached data
-------------------------

getPrimeOrder = method()
getPrimeOrder CpMackeyFunctor := ZZ => M ->(
    M.PrimeOrder
)

getUnderlyingModule = method()
getUnderlyingModule CpMackeyFunctor := Module => M ->(
    M.Underlying
)

getFixedModule = method()
getFixedModule CpMackeyFunctor := Module => M ->(
    M.Fixed
)

getRestriction = method()
getRestriction CpMackeyFunctor := Matrix => M ->(
    M.Res
)

getTransfer = method()
getTransfer CpMackeyFunctor := Matrix => M ->(
    M.Tr
)

getConjugation = method()
getConjugation CpMackeyFunctor := Matrix => M ->(
    M.Conj
)
