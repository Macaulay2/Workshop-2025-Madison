-- Input: A matrix and a scalar
-- Output: Boolean that gives whether the matrix defines a well-defined class of the unstable Grothendieck-Witt group. 

isWellDefinedGWu = method()
-- First version of this function treats the case where a is a Number, (eg. an element of CC_53, RR_53, QQ, or ZZ)
isWellDefinedGWu (Matrix, Number) := Boolean => (M, a) -> (

    -- Return false if the matrix isn't defined over a field
    if not isField ring M then return false;

    -- If matrix is defined over the complex numbers, allow scalar to be one of complex, real, rational, or integral. 
    if instance(ring M, ComplexField) and not (instance(ring a, ComplexField) or instance(ring a, RealField) or ring a===QQ or (a)===ZZ) then return false;

    -- If matrix is defined over the real numbers, allow scalar to be one of real, rational, or integral. 
    if instance(ring M, RealField) and not (instance(ring a, RealField) or ring a===QQ or ring a===ZZ) then return false;

    -- If matrix is defined over the rationals, allow scalar to be one of rational, or integral. 
    if ring M===QQ and not (ring a===QQ or ring a===ZZ) then return false;

    -- If matrix is defined over a finite field, allow scalar then the only scalars allowed are integral. The case of the scalar being over the same Galois field is treated in the next variant. 
    if instance(ring M, GaloisField) and not ring a===ZZ then return false;

    -- Then check that M is a well-defined element of GW(k)
    isWellDefinedGW(M)
    )

-- Second version of this function treats the case where a is a RingElement (eg. an element of a Galois field)
isWellDefinedGWu (Matrix, RingElement) := Boolean => (M, a) -> (

    -- Return false if the matrix isn't defined over a field
    if not isField ring M then return false;

    -- If matrix is defined over the complex numbers, allow scalar to be one of complex, real, rational, or integral. 
    if instance(ring M, ComplexField) or instance(ring M, RealField) or ring M===QQ then return false;

    -- If matrix is defined over a finite field, allow scalar to be an element of that Galois field. The case of a being an integer is treated in the previous variant. 
    if instance(ring M, GaloisField) and not (ring a===ZZ or (instance(ring a, GaloisField) and (ring M).order == (ring a).order)) then return false;

    -- If matrix is defined over an arbitrary field, allow scalar to be either integral or an element of that field. 
    if not (ring a===ring M or ring a===ZZ) then return false;

    -- Then check that M is a well-defined element of GW(k)
    isWellDefinedGW(M)
    )

-- We define UnstableGrothendieckWittClass to be a new type, meant to represent the isomorphism class 
-- of a nondegenerate symmetric bilinear form over a field of characteristic not 2 together with the data of a scalar.

UnstableGrothendieckWittClass = new Type of HashTable
UnstableGrothendieckWittClass.synonym = "Unstable Grothendieck-Witt Class"

-- Input: Either a matrix M or a matrix-scalar pair (M,a) representing a well-defined element of the unstable Grothendieck-Witt group. 
-- Output: The GrothendieckWittClass representing the symmetric bilinear form determined by M

makeGWuClass = method()
-- First version of this function treats the case of an input (M,a) where a is a Number, (eg. an element of CC_53, RR_53, QQ, or ZZ)
makeGWuClass (Matrix, Number) := UnstableGrothendieckWittClass => (M, a) -> (
   if isWellDefinedGWu (M, a) then (
        new UnstableGrothendieckWittClass from {
            symbol matrix => M,
            symbol cache => new CacheTable,
            symbol scalar => sub(a, ring M)
            }
        )
    else (
        error "makeGWuClass called on a pair that does not produce a well-defined element of the unstable Grothendieck-Witt group.";
	)
    )
-- Second version of this function treats the case of an input (M,a) where a is a Number (eg. an element of a Galois field)
makeGWuClass (Matrix, RingElement) := UnstableGrothendieckWittClass => (M, a) -> (
   if isWellDefinedGWu (M, a) then (
        new UnstableGrothendieckWittClass from {
            symbol matrix => M,
            symbol cache => new CacheTable,
            symbol scalar => sub(a, ring M)
            }
        )
    else (
        error "makeGWuClass called on a pair that does not produce a well-defined element of the unstable Grothendieck-Witt group.";
	)
    )
-- Third version of this function treats the case of an input M, where a is assumed to be the determinant of M. 
makeGWuClass (Matrix) := UnstableGrothendieckWittClass => (M) -> (
   if isWellDefinedGWu (M, det(M)) then (
        new UnstableGrothendieckWittClass from {
            symbol matrix => M,
            symbol cache => new CacheTable,
            symbol scalar => sub(det(M), ring M)
            }
        )
    else (
        error "makeGWuClass called on a matrix that does not produce a well-defined Grothendieck-Witt class.";
	)
    )


-- Input: An UnstableGrothendieckWittClass
-- Output: A string for printing the underlying scalar

getScalar = method()
getScalar UnstableGrothendieckWittClass := alpha -> (
    alpha.scalar
)

-- Input: An unstable Grothendieck-Witt class beta
-- Output: The base field of beta

getBaseField UnstableGrothendieckWittClass := Ring => beta -> (
    ring getMatrix beta
    )

-- Input: A GrothendieckWittClass representing a symmetric bilinear form determined by a matrix M
-- Output: The matrix M

getMatrix UnstableGrothendieckWittClass := Matrix => alpha -> (
    alpha.matrix
    )

-- Input: Two Grothendieck-Witt classes beta and gamma over the same field
-- Output: The direct sum of beta and gamma

addGWu = method()
addGWu (UnstableGrothendieckWittClass,UnstableGrothendieckWittClass) := UnstableGrothendieckWittClass => (beta,gamma) -> (
    Kb := getBaseField beta;
    Kg := getBaseField gamma;
    
    -- Galois field case
    if instance(Kb, GaloisField) and instance(Kg, GaloisField) then (
	-- Return an error if the underlying fields of the two classes are different
	if not Kb.order == Kg.order then
	    error "these classes have different underlying fields";
	return makeGWuClass(getMatrix beta ++ sub(getMatrix gamma, Kb), getScalar beta * sub(getScalar gamma, Kb));
	);
    
    -- Remaining cases
    if not Kb === Kg then
	error "these classes have different underlying fields";
    makeGWuClass(getMatrix beta ++ getMatrix gamma, getScalar beta * getScalar gamma)
    )

addGWuDivisorial = method()
addGWuDivisorial (List, List) := UnstableGrothendieckWittClass => (classList, rootList) -> (
    n := #classList;

    -- Return an error if list of roots is of different size than list of classes
    if n != #rootList then
        error "need same number of classes and roots";

    -- Return an error if lists are empty
    if n == 0 then
        error "the empty sum is the additive identity of the unstable Grothendieck-Witt group over the field of interes; please construct this as makeGWuClass(matrix(k,{}),1)";

    -- Return an error if the base fields are different for the list of GWu classes
    baseFieldList := apply(classList, getBaseField);
    matrixList := apply(classList, getMatrix);
    scalarList := apply(classList, getScalar);
    multiplicityList := apply(classList, i -> rank(getMatrix(i)));
    isGaloisField := apply(baseFieldList, i -> instance(i, GaloisField));
    if (not instance(baseFieldList#0, GaloisField) and not same baseFieldList) or (isGaloisField#0 and (not same isGaloisField or not same apply(baseFieldList, i -> i.order))) then 
        error "the list of GWu classes should have the same base field";
    if not fieldsAreCompatible(baseFieldList, rootList) then
        error "the roots must be in the base field of the classes";
    newForm := directSum(matrixList);
    newScalar := product(scalarList);
    for i from 0 to n-1 do (
        for j from i+1 to n-1 do (
            newScalar = newScalar * (rootList#i-rootList#j)^(2*multiplicityList#i*multiplicityList#j);
        );
    );
    makeGWuClass(newForm,newScalar)
    )

fieldsAreCompatible = method()
fieldsAreCompatible (List, List) := Boolean => (baseFieldList, rootList) -> (
    n := #baseFieldList;
    for i from 0 to n-1 do (
        -- If matrix is defined over the complex numbers, allow root to be one of complex, real, rational, or integral. 
        if instance(baseFieldList#i, ComplexField) and not (instance(ring rootList#i, ComplexField) or instance(ring rootList#i, RealField) or ring rootList#i===QQ or ring rootList#i===ZZ) then return false;

        -- If matrix is defined over the real numbers, allow scalar to be one of real, rational, or integral. 
        if instance(baseFieldList#i, RealField) and not (instance(ring rootList#i, RealField) or ring rootList#i===QQ or ring rootList#i===ZZ) then return false;

        -- If matrix is defined over the rationals, allow scalar to be one of rational, or integral. 
        if baseFieldList#i===QQ and not (ring rootList#i===QQ or ring rootList#i===ZZ) then return false;

        -- If matrix is defined over a finite field, allow scalar then the only scalars allowed are integral. The case of the scalar being over the same Galois field is treated in the next variant. 
        if instance(baseFieldList#i, GaloisField) and not (ring rootList#i===ZZ or (instance(ring rootList#i, GaloisField) and (baseFieldList#i).order == (ring rootList#i).order)) then return false;
    );
    true
)
