loadPackage "NumericalAlgebraicGeometry";

-- Input: A reduced pointed rational function q = f/g
-- Output: A pair (M,a) where M is a matrix and a is a scalar (the determinant of M)

getGlobalUnstableA1Degree = method(Options => {linearTolerance => 1e-6})

getGlobalUnstableA1Degree RingElement := UnstableGrothendieckWittClass => opts -> q -> (

    R := ring q;

    -- Extract numerator f from q
    f := numerator(sub(q, frac R));
       
    -- Extract numerator g from q and first normalize g
    g := denominator(sub(q, frac R));
    g = g/leadCoefficient(f);

    -- then normalize the leading coefficient of f
    f = f/leadCoefficient(f);
    
    -- Get the underlying ring and ensure it is a field
    kk := coefficientRing ring(f);
    if not isField kk then kk = toField kk;
    
    -- Check whether the rational function has isolated zeros
    if dim ideal(f) > 0 then 
        error "rational function does not have isolated zeros";
	
    -- Check whether the number of variables matches the number of polynomials
    S := ring f;
    u := (gens ring q)#0;

    if #(gens S) != 1 then
        error "the number of variables does not match the number of polynomials";    

    if (degree f)#0 <= (degree g)#0 then
        error "the rational function is not pointed";
    
    -- If the field is CC, ask the user to run the computation using the method getGlobalUnstableA1Degree(f,g)
    if instance(kk, ComplexField) then
        error "getGlobalUnstableA1Degree does not work over the complex numbers for a rational function f/g. Instead, use getGlobalUnstableA1Degree(f,g).";    

    -- If the field is RR, ask the user to run the computation over QQ instead and then base change to RR
    if instance(kk, RealField) then error "getGlobalUnstableA1Degree method does not work over the reals. Instead, define the polynomials over QQ to output an UnstableGrothendieckWittClass. Then extract the form, base change it to RR, and run getSumDecomposition().";    

    -- Initialize a polynomial ring in X and Y in which to compute the Bezoutian
    X := local X;
    Y := local Y;
    R' := kk[X,Y];
    -- R' := QQ[i][X,Y]; 
    
    fX := sub(f,{u => X});
    fY := sub(f,{u => Y});
    gX := sub(g,{u => X});
    gY := sub(g,{u => Y});

    D := lift((fX * gY - fY * gX)/(X-Y),R');
    
    m := degree(X,D);
    n := degree(Y,D);
        
    B := mutableMatrix id_(kk^(m+1));  
    
    for i from 0 to m do(
	for j from 0 to n do
    	B_(i,j) = coefficient(X^i*Y^j,D)
	);
    
    makeGWuClass matrix B
)

-- Input: A pair (f,g) of univariate polynomials for which f/g is a pointed rational function
-- Output: A pair (M,a) where M is a matrix and a is a scalar (the determinant of M)

getGlobalUnstableA1Degree (RingElement, RingElement) := UnstableGrothendieckWittClass => opts -> (f, g) -> (

   linTol := opts.linearTolerance;

   if not ((ring f === ring g) and length gens ring f == 1) then
        error "the two polynomials must be in the same univariate polynomial ring";

    R := ring f;

    -- Normalize the leading coefficient of f
    f = f/leadCoefficient(f);
    
    -- Get the underlying ring and ensure it is a field
    kk := coefficientRing ring(f);
    if not isField kk then kk = toField kk;
    
    -- Check whether the rational function has isolated zeros
    if dim ideal(f) > 0 then 
        error "rational function does not have isolated zeros";
	
    -- Check whether the number of variables matches the number of polynomials
    S := ring f;
    u := (gens ring f)#0;

    if #(gens S) != 1 then
        error "the number of variables does not match the number of polynomials";    
    
    -- Check if rational f/g  function is pointed
    if (degree f)#0 <= (degree g)#0 then
        error "the rational function is not pointed";
    
    if instance(kk, ComplexField) then (

	 -- Use NumericalAlgebraicGeometry to compute roots (numerical in CC)
	 r1 := roots f;
	 r2 := roots g;

	 -- Cancel common roots
	 removeCommonApprox := (L1,L2,tol) -> (
	     L1new := {};
	     L2new := L2;
	     for r in L1 list (
		 i := position(L2new, s -> abs(r-s) < tol);
		 if i === null then (
		     L1new = append(L1new, r)
		     ) else (
		     L2new = (take(L2new,i)) | (drop(L2new,i+1));
		     )
		 );
	     (L1new, L2new)
	     );

	 -- Cancel common roots numerically
	 (r1, r2) = removeCommonApprox(r1, r2, linTol);

	 -- Rebuild cleaned numerator and denominator
	 x := (gens ring f)#0;
	 fr := product(r1, r -> (x - r));
	 fr = fr/leadCoefficient(fr);
	 gr := (leadCoefficient g)*product(r2, r -> (x - r));
   	
	 makeGWuClass(
            id_(CC^((degree fr)#0)),
            promote((-1)^(((degree fr)#0^2 - (degree fr)#0)/2), kk) * getResultant(fr, gr)
	    )
    ) else if instance(kk, RealField) then (
        error "getGlobalUnstableA1Degree method does not work over the reals."
    ) else (
    getGlobalUnstableA1Degree(f/g)
	 )
    )

-- Input: A pair (q,r) where q is a rational function and r is a root of q
-- Output: An unstable Grothendieck-Witt class

getLocalUnstableA1Degree = method()
getLocalUnstableA1Degree (RingElement, Number) := (UnstableGrothendieckWittClass) => (q, r) -> (

    if not (instance(ring q, PolynomialRing) or instance(ring q, FractionField)) then
        error "input must be in polynomial ring or fraction field";
        
    kk := coefficientRing ring q;

    if not (kk === QQ or (instance(kk, GaloisField) and kk.char != 2)) then 
        error "only implemented over QQ and finite fields of characteristic not 2";

    -- If the base field is QQ, allow the root to be integer or rational
    if kk === QQ and not (ring r === QQ or ring r === ZZ) then error "root not from the base field of the polynomial";

    -- If the base field is a finite field, allow the root to be integer, rational, or from the same finite field
    if instance(kk, GaloisField) and not (ring r === QQ or ring r === ZZ or (instance(ring r, GaloisField) and kk.order == (ring r).order)) then error "root not from the base field of the polynomial";

    if numgens ring q != 1 then error "must input function of one variable";

    u := (gens ring q)#0;

    q = sub(q, frac ring q);

    -- Extract numerator f from q
    f := numerator(q);
       
    -- Extract denominator g from q
    g := denominator(q);
    
    -- Check whether the rational function has isolated zeros
    if dim ideal(f) > 0 then 
        error "rational function does not have isolated zeros";
	
    -- Check whether the number of variables matches the number of polynomials
    if not f(r) == 0 then
        error "the field element is not a zero of the function";

    if (degree f)#0 <= (degree g)#0 then
        error "the rational function is not pointed"; 

    m := getMultiplicity(f, r);

    F := (u - sub(r, frac ring q))^m * g/f;

    makeAntidiagonalUnstableForm(kk, F(r), m)
)

getLocalUnstableA1Degree (RingElement, RingElement) := (UnstableGrothendieckWittClass) => (q, r) -> (

    if not (instance(ring q, PolynomialRing) or instance(ring q, FractionField)) then
        error "input must be in polynomial ring or fraction field";
        
    kk := coefficientRing ring q;

    if not (kk === QQ or (instance(kk, GaloisField) and kk.char != 2)) then 
        error "only implemented over QQ and finite fields of characteristic not 2";

    -- If the base field is QQ, allow the root to be integer or rational
    if kk === QQ and not (ring r === QQ or ring r === ZZ) then error "root not from the base field of the polynomial";

    -- If the base field is a finite field, allow the root to be integer, rational, or from the same finite field
    if instance(kk, GaloisField) and not (ring r === QQ or ring r === ZZ or (instance(ring r, GaloisField) and kk.order == (ring r).order)) then error "root not from the base field of the polynomial";

    if numgens ring q != 1 then error "must input function of one variable";

    u := (gens ring q)#0;

    q = sub(q, frac ring q);

    -- Extract numerator f from q
    f := numerator(q);
       
    -- Extract denominator g from q
    g := denominator(q);
    
    -- Check whether the rational function has isolated zeros
    if dim ideal(f) > 0 then 
        error "rational function does not have isolated zeros";
	
    -- Check whether the number of variables matches the number of polynomials
    if not f(r) == 0 then
        error "the field element is not a zero of the function";
    
    if (degree f)#0 <= (degree g)#0 then
        error "the rational function is not pointed"; 

    m := getMultiplicity(f, r);

    F := (u - sub(r, frac ring q))^m * g/f;

    makeAntidiagonalUnstableForm(kk, F(r), m)
)

-- Variant that takes in numerator and denominator separately
getLocalUnstableA1Degree (RingElement, RingElement, Number) := (UnstableGrothendieckWittClass) => (f, g, r) -> (

    if not (instance(ring f, PolynomialRing) and instance(ring g, PolynomialRing) and ring f === ring g) then
        error "both input polynomials must be defined over the same univariate polynomial ring";
        
    kk := coefficientRing ring f;

    if not (instance(kk, ComplexField) or kk === QQ or (instance(kk, GaloisField) and kk.char != 2)) then 
        error "only implemented over CC, QQ, and finite fields of characteristic not 2";

    -- If the base field is CC, allow the root to be complex, real, rational, or an integer. 
    if instance(kk, ComplexField) and not (instance(ring r, ComplexField) or instance(ring r, RealField) or ring r === QQ or ring r === ZZ) then error "root not from the base field of the polynomials";

    -- If the base field is QQ, allow the root to be integer or rational
    if kk === QQ and not (ring r === QQ or ring r === ZZ) then error "root not from the base field of the polynomial";

    -- If the base field is a finite field, allow the root to be integer, rational, or from the same finite field
    if instance(kk, GaloisField) and not (ring r === QQ or ring r === ZZ or (instance(ring r, GaloisField) and kk.order == (ring r).order)) then error "root not from the base field of the polynomial";

    if numgens ring f != 1 then error "must input function of one variable";

    
    -- Check whether the rational function has isolated zeros
    if dim ideal(f) > 0 then 
        error "rational function does not have isolated zeros";
	
    -- Check whether the number of variables matches the number of polynomials
    if not f(r) == 0 then
        error "the field element is not a zero of the function";
    
    u := (gens ring f)#0;

    if (degree f)#0 <= (degree g)#0 then
        error "the rational function is not pointed"; 
    
    if instance(kk, ComplexField) then
        return getLocalUnstableA1DegreeCC(f, g, r);

    m := getMultiplicity(f, r);

    F := (u - sub(r, frac ring f))^m * g/f;

    makeAntidiagonalUnstableForm(kk, F(r), m)
)

-- Variant that takes in numerator and denominator separately
getLocalUnstableA1Degree (RingElement, RingElement, RingElement) := (UnstableGrothendieckWittClass) => (f, g, r) -> (

    if not (instance(ring f, PolynomialRing) and instance(ring g, PolynomialRing) and ring f === ring g) then
        error "both input polynomials must be defined over the same univariate polynomial ring";
        
    kk := coefficientRing ring f;

    if not (kk === QQ or (instance(kk, GaloisField) and kk.char != 2)) then 
        error "QQ, and finite fields of characteristic not 2";

    -- If the base field is CC, allow the root to be complex, real, rational, or an integer. 
    if instance(kk, ComplexField) and not (instance(ring r, ComplexField) or instance(ring r, RealField) or ring r === QQ or ring r === ZZ) then error "root not from the base field of the polynomials";

    -- If the base field is QQ, allow the root to be integer or rational
    if kk === QQ and not (ring r === QQ or ring r === ZZ) then error "root not from the base field of the polynomial";

    -- If the base field is a finite field, allow the root to be integer, rational, or from the same finite field
    if instance(kk, GaloisField) and not (ring r === QQ or ring r === ZZ or (instance(ring r, GaloisField) and kk.order == (ring r).order)) then error "root not from the base field of the polynomial";

    if numgens ring f != 1 then error "must input function of one variable";

    
    -- Check whether the rational function has isolated zeros
    if dim ideal(f) > 0 then 
        error "rational function does not have isolated zeros";
	
    -- Check whether the number of variables matches the number of polynomials
    if not f(r) == 0 then
        error "the field element is not a zero of the function";
    
    u := (gens ring f)#0;

    if (degree f)#0 <= (degree g)#0 then
        error "the rational function is not pointed"; 

    m := getMultiplicity(f, r);

    F := (u - sub(r, frac ring f))^m * g/f;

    makeAntidiagonalUnstableForm(kk, F(r), m)
)

-- Input: A rational function f/g, a root of f, and the multiplicity of that root
-- Output: An unstable Grothendieck-Witt class
getLocalUnstableA1DegreeCC = method()
getLocalUnstableA1DegreeCC(RingElement, RingElement, Number) := UnstableGrothendieckWittClass => (f, g, r) -> (
    Sf := solveSystem {f};
    Sg := solveSystem {g};

    rootsf := apply(Sf, i -> i.Coordinates);
    rootsg := apply(Sg, i -> i.Coordinates);

    rootsfNotr := select(rootsf, i -> not areEqual(i#0, sub(r, CC_53)));

    outputFormRank := number(rootsf, i -> areEqual(i#0, sub(r, CC_53)));

    LDdenom := sub(1, CC_53);

    -- compute the denominator of the local degree as the evaluation of the product of (x-ri) where ri range over the roots not equal to r
    for i from 0 to (#rootsfNotr) - 1 do (
        LDdenom = LDdenom * (r - rootsfNotr#i#0)
    );

    makeGWuClass(id_(CC_53^(outputFormRank)), (-1)^((outputFormRank^2 - outputFormRank)/2)*(g(r)/LDdenom)^outputFormRank)        
)


-- Input: a polynomial in one variable and a root
-- Output: multiplicity of the polynomial

getMultiplicity = method()
getMultiplicity(RingElement, Number) := ZZ => (f, r) -> (
    -- return an error if the polynomial isn't in one variable, or is a polynomial at all
    if not instance(ring f, PolynomialRing) or numgens ring f != 1 then
        error "need polynomial with one variable";

    -- return an error if the root isn't an element of the base field of the polynomial
    try r = sub(r, coefficientRing ring f) else error "entered root must be in base field of polynomial";

    var := (gens ring f)#0;
    multiplicity := 0;

    --for each root add one more to the multiplicity
    while f(r) == 0 do (
        multiplicity += 1;
        f = sub(f/(var - r), ring f);
    );
    multiplicity
)

getMultiplicity(RingElement, RingElement) := ZZ => (f, r) -> (
    -- return an error if the polynomial isn't in one variable, or is a polynomial at all
    if not instance(ring f, PolynomialRing) or numgens ring f != 1 then
        error "need polynomial with one variable";

    -- return an error if the root isn't an element of the base field of the polynomial
    try r = sub(r, coefficientRing ring f) else error "entered root must be in base field of polynomial";

    var := (gens ring f)#0;
    multiplicity := 0;

    --for each root add one more to the multiplicity
    while f(r) == 0 do (
        multiplicity += 1;
        f = sub(f/(var - r), ring f);
    );
    multiplicity
)
