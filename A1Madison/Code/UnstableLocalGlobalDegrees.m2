-------------------------------------------------------------
-- unstable A1-Brouwer degree methods (task 3.12 from Overleaf)
-------------------------------------------------------------

-- Input: A pointed rational function q = f/g

-- Output: A pair (M,a) where M is a matrix and a is a scalar (the determinant of M)

getGlobalUnstableA1Degree = method()
getGlobalUnstableA1Degree RingElement := (Matrix,Number) => q -> (

    R := ring q;

    -- Extract numerator f from q
    f := numerator(sub(q, frac R));
       
    -- Extract numerator g from q
    g := denominator(sub(q, frac R));
    
    -- Get the underlying ring and ensure it is a field
    kk := coefficientRing ring(f);
    if not isField kk then kk = toField kk;
    
    -- Check whether the rational function has isolated zeros
    if dim ideal(f) > 0 then 
        error "rational function does not have isolated zeros";
	
    -- Check whether the number of variables matches the number of polynomials
    S := ring f;
    if #(gens S) != 1 then
        error "the number of variables does not match the number of polynomials";     
    -- If the field is CC, output the Grothendieck-Witt class of an identity matrix of the appropriate rank
    if instance(kk, ComplexField) then (
    	rankAlgebra := getGlobalAlgebraRank {q};
    	return makeGWuClass id_(CC^rankAlgebra);
        );

    -- If the field is RR, ask the user to run the computation over QQ instead and then base change to RR
    if instance(kk, RealField) then error "getGlobalUnstableA1Degree method does not work over the reals. Instead, define the polynomials over QQ to output an unstableGrothendieckWittClass. Then extract the matrix, base change it to RR, and run getSumDecomposition().";    

    -- Initialize a polynomial ring in X and Y in which to compute the Bezoutian
    X := local X;
    Y := local Y;
    R' := kk[X,Y];
    
    fX := sub(f,{x => X});
    fY := sub(f,{x => Y});
    gX := sub(g,{x => X});
    gY := sub(g,{x => Y});
    
    D := lift((fX * gY - fY * gX)/(X-Y),R');
    
    m := degree (X,D);
    n := degree (Y,D);
        
    B := mutableMatrix id_(QQ^(m+1));  
    
    for i from 0 to m do(
	for j from 0 to n do
    	B_(i,j) := coefficient(X^i*Y^j,D)
	);
    
     makeGWuClass matrix B
     )

-- Input: A rational function f/g, a root of f, and the multiplicity of that root

-- Output: A pair (M,a) where M is a matrix and a is a scalar

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

    u := (gens ring q)#0

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

    m := getMultiplicity(f, r);

    F := (u - sub(r, frac ring q))^m * g/f;

    makeAntidiagonalUnstableForm(kk, F(r), m)
)

-- Input: a polynomial in one variable and a root
-- Output: multiplicity of the polynomial

getMultiplicity = method()
getMultiplicity(RingElement, Number) := ZZ => (f, r) ->(
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

getMultiplicity(RingElement, RingElement) := ZZ => (f, r) ->(
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
