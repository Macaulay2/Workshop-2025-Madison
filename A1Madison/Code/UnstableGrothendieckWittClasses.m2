-- Input: A matrix
-- Output: Boolean that gives whether the matrix defines a nondegenerate symmetric bilinear form over a field of characteristic not 2

isWellDefinedGWu = method()
-- First version of this function treats the case where a is a Number, (eg. an element of CC_53, RR_53, QQ, or ZZ)
isWellDefinedGWu (Matrix, Number) := Boolean => (M, a) -> (

    -- Return false if the matrix isn't defined over a field
    if not isField ring M then return false;

    -- If matrix is defined over the complex numbers, allow scalar to be one of complex, real, rational, or integral. 
    if instance(ring(M), ComplexField) and not (instance(ring(a), ComplexField) or instance(ring(a), RealField) or ring(a)===QQ or (a)===ZZ) then return false;

    -- If matrix is defined over the real numbers, allow scalar to be one of real, rational, or integral. 
    if instance(ring(M), RealField) and not (instance(ring(a), RealField) or ring(a)===QQ or ring(a)===ZZ) then return false;

    -- If matrix is defined over the rationals, allow scalar to be one of rational, or integral. 
    if ring(M)===QQ and not (ring(a)===QQ or ring(a)===ZZ) then return false;

    -- If matrix is defined over a finite field, allow scalar to be integral or an element of that Galois field. 
    if instance(ring(M), GaloisField) and not ring(a)===ZZ then return false;

    -- Then check that M is a well-defined element of GW(k)
    isWellDefinedGW(M)
    )
isWellDefinedGWu (Matrix, RingElement) := Boolean => (M, a) -> (

    -- Return false if the matrix isn't defined over a field
    if not isField ring M then return false;

    -- If matrix is defined over the complex numbers, allow scalar to be one of complex, real, rational, or integral. 
    if instance(ring(M), ComplexField) or instance(ring(M), RealField) or ring(M)===QQ return false;

    -- If matrix is defined over a finite field, allow scalar to be integral or an element of that Galois field. 
    if instance(ring(M), GaloisField) and not (ring(a)===ZZ or (instance(ring(a), GaloisField) and ring(M).order == ring(a).order)) then return false;

    -- If matrix is defined over an arbitrary field, allow scalar to be either integral or an element of that field. 
    if not (ring(a)===ring(M) or ring(a)===ZZ) then return false;

    -- Then check that M is a well-defined element of GW(k)
    isWellDefinedGW(M)
    )

-- We define GrothendieckWittClass to be a new type, meant to represent the isomorphism class 
-- of a nondegenerate symmetric bilinear form over a field of characteristic not 2

UnstableGrothendieckWittClass = new Type of HashTable
UnstableGrothendieckWittClass.synonym = "Unstable Grothendieck-Witt Class"

-- Input: A matrix M representing a nondegenerate symmetric bilinear form over a field of characteristic not 2
-- Output: The GrothendieckWittClass representing the symmetric bilinear form determined by M

makeGWuClass = method()
makeGWuClass (Matrix, Number) := UnstableGrothendieckWittClass => (M, a) -> (
   if isWellDefinedGWu (M, a) then (
        new UnstableGrothendieckWittClass from {
            symbol matrix => M,
            symbol cache => new CacheTable,
            symbol scalar => sub(a, ring(M))
            }
        )
    else (
        error "makeGWuClass called on a pair that does not produce a well-defined element of the unstable Grothendieck-Witt group.";
	)
    )
makeGWuClass (Matrix, RingElement) := UnstableGrothendieckWittClass => (M, a) -> (
   if isWellDefinedGWu (M, a) then (
        new UnstableGrothendieckWittClass from {
            symbol matrix => M,
            symbol cache => new CacheTable,
            symbol scalar => sub(a, ring(M))
            }
        )
    else (
        error "makeGWuClass called on a pair that does not produce a well-defined element of the unstable Grothendieck-Witt group.";
	)
    )
makeGWuClass (Matrix) := UnstableGrothendieckWittClass => (M) -> (
   if isWellDefinedGWu (M, det(M)) then (
        new UnstableGrothendieckWittClass from {
            symbol matrix => M,
            symbol cache => new CacheTable,
            symbol scalar => sub(det(M), ring(M))
            }
        )
    else (
        error "makeGWuClass called on a pair that does not produce a well-defined element of the unstable Grothendieck-Witt group.";
	)
    )
makeGWuClass (Matrix) := UnstableGrothendieckWittClass => (M) -> (
   if isWellDefinedGWu (M, det(M)) then (
        new UnstableGrothendieckWittClass from {
            symbol matrix => M,
            symbol cache => new CacheTable,
            symbol scalar => sub(det(M), ring(M))
            }
        )
    else (
        error "makeGWuClass called on a pair that does not produce a well-defined element of the unstable Grothendieck-Witt group.";
	)
    )