CpMackeyFunctor = new Type of HashTable
CpMackeyFunctor.synonym = "Cp Mackey Functor"

isWellDefinedCpMackeyFunctor = method()
isWellDefinedCpMackeyFunctor CpMackeyFunctor := Boolean => M ->  (
    if not (M#?Res and M#?Tr and M#?Conj and M#?Underlying and M#?Fixed) then return false;
    
    -- Check fixed and underlying modules are Z-modules
    if not (isModule M.Fixed and isModule M.Underlying) then return false;
    if not (ring M.Fixed === ZZ and ring M.Underlying === ZZ) then return false;

    -- Check source and target of restriction are correct
    if not source(M.Res) == M.Fixed then return false;
    if not target(M.Res) == M.Underlying then return false;

    -- Check source and target of transfer are correct
    if not source(M.Tr) == M.Underlying then return false;
    if not target(M.Tr) == M.Fixed then return false;
    
    -- Check source and target of conjugation are correct
    if not source(M.Conj) == M.Underlying then return false;
    if not target(M.Conj) == M.Underlying then return false;

    return true
    -- TODO: check all axioms
)


makeCpMackeyFunctor = method()
-- Ordering for the input is:
-- 1. List of modules
--    [Underlying, Fixed]
-- 2. Restriction
-- 3. Transfer
-- 4. Conjugation
makeCpMackeyFunctor(List,Matrix,Matrix,Matrix) := CpMackeyFunctor => (L,R,T,C) ->(
    U := L#0;
    F := L#1;
    M := new CpMackeyFunctor from {
        symbol Underlying => U,
        symbol Fixed => F,
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

-- todo: test this still works when the first input is an Array