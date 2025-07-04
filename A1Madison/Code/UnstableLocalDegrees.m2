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

    x := (gens ring q)#0;

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

    F := (x-sub(r, frac ring q))^m * g/f;

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

    x := (gens ring q)#0;

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

    F := (x-sub(r, frac ring q))^m * g/f;

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