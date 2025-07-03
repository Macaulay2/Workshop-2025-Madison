-- Input: a polynomial in one variable and a root
-- Output: multiplicity of the polynomial

getMultiplicity = method()
getMultiplicity(RingElement, Number) := ZZ => (f, r) ->(
    -- return an error if the polynomial isn't in one variable, or is a polynomial at all
    if not instance(ring f, PolynomialRing) or numgens ring f != 1 then
        error "need polynomial with one variable";

    -- return an error if the root isn't an element of the base field of the polynomial
    try d = sub(r, coefficientRing ring f) else error "entered root must be in base field of polynomial";

    var := (gens ring f)#0;
    multiplicity := 0;

    --for each root add one more to the multiplicity
    while f(r) == 0 do (
        multiplicity += 1;
        f = sub(f/(var - r), ring f);
    );
    multiplicity
)

getMultiplicity = method()
getMultiplicity(RingElement, RingElement) := ZZ => (f, r) ->(
    -- return an error if the polynomial isn't in one variable, or is a polynomial at all
    if not instance(ring f, PolynomialRing) or numgens ring f != 1 then
        error "need polynomial with one variable";

    -- return an error if the root isn't an element of the base field of the polynomial
    try d = sub(r, coefficientRing ring f) else error "entered root must be in base field of polynomial";

    var := (gens ring f)#0;
    multiplicity := 0;

    --for each root add one more to the multiplicity
    while f(r) == 0 do (
        multiplicity += 1;
        f = sub(f/(var - r), ring f);
    );
    multiplicity
)