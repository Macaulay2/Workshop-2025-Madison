-- Input: A matrix and a scalar
-- Output: Boolean that gives whether the matrix defines a well-defined class of the unstable Grothendieck-Witt group. 

isWellDefinedGWu = method()
-- First version of this function treats the case where a is a Number, (eg. an element of CC_53, RR_53, QQ, or ZZ)
isWellDefinedGWu (Matrix, Number) := Boolean => (M, a) -> (

    -- Return false if the matrix isn't defined over a field
    if not isField ring M then return false;

    -- Return false if a is zero in the field
    if a == 0 then return false;

    -- If matrix is defined over the complex numbers, allow scalar to be one of complex, real, rational, or integral. 
    if instance(ring M, ComplexField) and not (instance(ring a, ComplexField) or instance(ring a, RealField) or ring a === QQ or ring a === ZZ) then return false;

    -- If matrix is defined over the real numbers, allow scalar to be one of real, rational, or integral. 
    if instance(ring M, RealField) and not ((instance(ring a, RealField) or ring a === QQ or ring a === ZZ) and (sign(a) == sign(det(M)))) then return false;

    -- If matrix is defined over the rationals, allow scalar to be one of rational, or integral. 
    if ring M === QQ and not ((ring a === QQ or ring a === ZZ) and (getSquarefreePart det M == getSquarefreePart a)) then return false;

    -- If matrix is defined over a finite field, allow scalar then the only scalars allowed are integral. The case of the scalar being over the same Galois field is treated in the next variant. 
    if instance(ring M, GaloisField) and not (ring a === ZZ and isGFSquare(det M) == isGFSquare(sub(a, ring M))) then return false;

    -- Then check that M is a well-defined element of GW(k)
    isWellDefinedGW M
    )

-- Second version of this function treats the case where a is a RingElement (eg. an element of a Galois field)
isWellDefinedGWu (Matrix, RingElement) := Boolean => (M, a) -> (

    -- Return false if the matrix isn't defined over a field
    if not isField ring M then return false;

    -- Return false if a is zero in the field
    if a == 0 then return false;

    -- If matrix is defined over the complex numbers, allow scalar to be one of complex, real, rational, or integral. 
    if instance(ring M, ComplexField) or instance(ring M, RealField) or ring M === QQ then return false;

    -- If matrix is defined over a finite field, allow scalar to be an element of that Galois field. The case of a being an integer is treated in the previous variant. 
    if instance(ring M, GaloisField) and not (instance(ring a, GaloisField) and (ring M).order == (ring a).order and isGFSquare(det M) == isGFSquare(a)) then return false

    -- If matrix is defined over an arbitrary field, scalars being equal to the determinant of the matrix are allowed automatically. 
    else if ring a === ring M and det(M) =!= a then print "Warning, the function is not able to verify if the determinant of M and a agree up to squares."

    else if not instance(ring M, GaloisField) and ring a =!= ring M then return false;

    -- Then check that M is a well-defined element of GW(k)
    isWellDefinedGW M
    )

-- We define UnstableGrothendieckWittClass to be a new type, meant to represent the isomorphism class 
-- of a nondegenerate symmetric bilinear form over a field of characteristic not 2 together with the data of a scalar.

UnstableGrothendieckWittClass = new Type of HashTable
UnstableGrothendieckWittClass.synonym = "Unstable Grothendieck-Witt Class"

-- Input: An UnstableGrothendieckWittClass
-- Output: A net for printing the underlying dat

net UnstableGrothendieckWittClass := Net => alpha -> (
    net (getMatrix alpha, getScalar alpha)
    )

-- Input: A GrothendieckWittClass
-- Output: A string for printing the underlying matrix

texMath UnstableGrothendieckWittClass := String => alpha -> (
    texMath (getMatrix alpha, getScalar alpha)
    )

-- Input: Either a matrix M or a matrix-scalar pair (M,a) representing a well-defined element of the unstable Grothendieck-Witt group. 
-- Output: The GrothendieckWittClass representing the symmetric bilinear form determined by M

makeGWuClass = method()

-- First version of this function treats the case of an input (M,a) where a is a Number (eg. an element of CC_53, RR_53, QQ, or ZZ)
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
   if isWellDefinedGWu (M, det M) then (
        new UnstableGrothendieckWittClass from {
            symbol matrix => M,
            symbol cache => new CacheTable,
            symbol scalar => sub(det M, ring M)
            }
        )
    else (
        error "makeGWuClass called on a matrix that does not produce a well-defined element of the unstable Grothendieck-Witt group.";
	)
    )

-- Fourth version of this function treats the case of an input a GrothendieckWittClass alpha and a Number (eg. an element of CC_53, RR_53, QQ, or ZZ)
makeGWuClass (GrothendieckWittClass, Number) := UnstableGrothendieckWittClass => (alpha, a) -> (
   if isWellDefinedGWu (getMatrix alpha, a) then (
        new UnstableGrothendieckWittClass from {
            symbol matrix => getMatrix alpha,
            symbol cache => new CacheTable,
            symbol scalar => sub(a, getBaseField alpha)
            }
        )
    else (
        error "makeGWuClass called on a pair that does not produce a well-defined element of the unstable Grothendieck-Witt group.";
	)
    )

-- Fifth version of this function treats the case of an input a GrothendieckWittClass alpha and a RingElement (eg. an element of a Galois field)
makeGWuClass (GrothendieckWittClass, RingElement) := UnstableGrothendieckWittClass => (alpha, a) -> (
   if isWellDefinedGWu (getMatrix alpha, a) then (
        new UnstableGrothendieckWittClass from {
            symbol matrix => getMatrix alpha,
            symbol cache => new CacheTable,
            symbol scalar => sub(a, getBaseField alpha)
            }
        )
    else (
        error "makeGWuClass called on a pair that does not produce a well-defined element of the unstable Grothendieck-Witt group.";
	)
    )

-- Sixth version of this function treats the case of an input a GrothendieckWittClass alpha, where a is assumed to be the determinant of the Gram matrix of alpha.  
makeGWuClass (GrothendieckWittClass) := UnstableGrothendieckWittClass => (alpha) -> (
   if isWellDefinedGWu (getMatrix alpha, det getMatrix alpha) then (
        new UnstableGrothendieckWittClass from {
            symbol matrix => getMatrix alpha,
            symbol cache => new CacheTable,
            symbol scalar => sub(det getMatrix alpha, getBaseField alpha)
            }
        )
    else (
        error "makeGWuClass called on a form that does not produce a well-defined element of the unstable Grothendieck-Witt group.";
	)
    )

-- Input: An unstable Grothendieck-Witt class
-- Output: Its stable part

makeStableGWClass = method()
makeStableGWClass (UnstableGrothendieckWittClass) := GrothendieckWittClass => alpha -> (
    makeGWClass getMatrix alpha
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

-- Input: List of GWu(k) classes, list of elements of k
-- Output: The divisorial sum of the GWu(k) classes as a GWu(k) class

addGWuDivisorial = method()
addGWuDivisorial (List, List) := UnstableGrothendieckWittClass => (classList, rootList) -> (
    n := #classList;
    baseFieldList := apply(classList, getBaseField);
    matrixList := apply(classList, getMatrix);
    scalarList := apply(classList, getScalar);
    multiplicityList := apply(classList, i -> rank(getMatrix(i)));
    isGaloisField := apply(baseFieldList, i -> instance(i, GaloisField));

    -- Return an error if list of roots is of different size than list of classes
    if n != #rootList then
        error "need same number of classes and roots";

    -- Return an error if lists are empty
    if n == 0 then
        error "the empty sum is the additive identity of the unstable Grothendieck-Witt group over the field of interest; please construct this as makeGWuClass(matrix(k,{}),1)";

    -- Return an error if the base fields are different for the list of GWu classes
    if (not instance(baseFieldList#0, GaloisField) and not same baseFieldList) or (isGaloisField#0 and (not same isGaloisField or not same apply(baseFieldList, i -> i.order))) then 
        error "the list of GWu classes should have the same base field";
    
    -- Return an error if the roots are not in the correct field
    if not fieldsAreCompatible(baseFieldList, rootList) then
        error "the roots must be in the base field of the classes";

    -- Create the sum matrix and scalar    
    newForm := directSum(matrixList);
    newScalar := product(scalarList);
    for i from 0 to n-1 do (
        for j from i+1 to n-1 do ( -- We require j > i 
            newScalar = newScalar * (rootList#i - rootList#j)^(2 * multiplicityList#i * multiplicityList#j);
        );
    );
    makeGWuClass(newForm,newScalar)
    )

-- Input: List of unstable Grothendieck-Witt classes, list of numbers/ring elements
-- Output: Boolean of whether the elements of the second list are elements of the field corresponding to the GWu class

fieldsAreCompatible = method()
fieldsAreCompatible (List, List) := Boolean => (baseFieldList, rootList) -> (
    n := #baseFieldList;
    for i from 0 to n-1 do (
        -- If matrix is defined over the complex numbers, allow root to be one of complex, real, rational, or integral. 
        if instance(baseFieldList#i, ComplexField) and not (instance(ring rootList#i, ComplexField) or instance(ring rootList#i, RealField) or ring rootList#i === QQ or ring rootList#i === ZZ) then return false;

        -- If matrix is defined over the real numbers, allow scalar to be one of real, rational, or integral. 
        if instance(baseFieldList#i, RealField) and not (instance(ring rootList#i, RealField) or ring rootList#i === QQ or ring rootList#i === ZZ) then return false;

        -- If matrix is defined over the rationals, allow scalar to be one of rational, or integral. 
        if baseFieldList#i === QQ and not (ring rootList#i === QQ or ring rootList#i === ZZ) then return false;

        -- If matrix is defined over a finite field, allow scalar then the only scalars allowed are integral. The case of the scalar being over the same Galois field is treated in the next variant. 
        if instance(baseFieldList#i, GaloisField) and not (ring rootList#i === ZZ or (instance(ring rootList#i, GaloisField) and (baseFieldList#i).order == (ring rootList#i).order)) then return false;
    );
    true
)

-- Input: An unstable Grothendieck-Witt class beta over QQ, RR, CC, or a finite field of characteristic not 2
-- Output: A diagonalized form of beta, with squarefree entries on the diagonal
getDiagonalClass UnstableGrothendieckWittClass := UnstableGrothendieckWittClass => beta -> (

    -- Check if the diagonal class has already been computed; if so, recall it from the cache
    if beta.cache.?getDiagonalClass then return beta.cache.getDiagonalClass;

    getDiagonalClassOfBetaMatrix := diagonalizeAndSimplifyViaCongruence getMatrix beta;

    -- The computed diagonal class gets stored in the cache
    beta.cache.getDiagonalClass = makeGWuClass(getDiagonalClassOfBetaMatrix, getScalar beta);
    makeGWuClass(getDiagonalClassOfBetaMatrix, getScalar beta)
    )

-- Input: Two unstable Grothendieck-Witt classes over CC, RR, QQ, or a finite field of characteristic not 2
-- Output: Boolean that gives whether the classes are isomorphic

isIsomorphicForm (UnstableGrothendieckWittClass,UnstableGrothendieckWittClass) := Boolean => (alpha,beta) -> (
    k1 := ring getScalar alpha;
    k2 := ring getScalar beta;

    -- Ensure both base fields are supported
    if not (instance(k1, ComplexField) or instance(k1, RealField) or k1 === QQ or (instance(k1, GaloisField) and k1.char != 2)) then
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
    if not (instance(k2, ComplexField) or instance(k2, RealField) or k2 === QQ or (instance(k2, GaloisField) and k2.char != 2)) then
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
    
    -- Over CC, the scalars are automatically in the same square class, so just need to check the forms
    if ((instance(k1, ComplexField) and instance(k2, ComplexField)) or (instance(k1, RealField) and instance(k2, RealField)) or (k1 === QQ and k2 === QQ)) then (
        return (isIsomorphicForm(getMatrix alpha, getMatrix beta) and getScalar alpha == getScalar beta);
        )
    
    -- Over a finite field, the scalars are in the same square class if and only if they are either both squares or both not squares 
    else if (instance(k1, GaloisField) and instance(k2, GaloisField) and k1.char !=2 and k2.char != 2 and k1.order == k2.order) then (
        return (isIsomorphicForm(getMatrix alpha, getMatrix beta) and getScalar alpha == sub(getScalar beta, k1));
        )
    -- If we get here, then the base fields are not the same
    else
	    error "Base fields are not the same";
    )
