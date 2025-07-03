-------------------------------------
-- unstable A1-Brouwer degree methods (task 3.12 from Overleaf)
-------------------------------------

-- Input: A rational function f/g

-- Output: A pair (M,a) where M is a matrix and a is a scalar

getLocalUnstableA1Degree = method()
getLocalUnstableA1Degree (RingElement, Number, ZZ) := (UnstableGrothendieckWittClass) => (q, r, m) -> (

    kk := coefficientRing ring q;

    x := (gens ring q)#0;

    R := kk[x];
    fracR := frac R;

    -- Extract numerator f from q
    f := numerator(sub(q, fracR));
       
    -- Extract numerator g from q
    g := denominator(sub(q, fracR));
    
    -- Get the underlying ring and ensure it is a field
    --kk := coefficientRing ring(f);
    if not isField kk then kk = toField kk;
    
    -- Check whether the rational function has isolated zeros
    if dim ideal(f) > 0 then 
        error "rational function does not have isolated zeros";
	
    -- Check whether the number of variables matches the number of polynomials
    
    S := ring f;
    if #(gens S) != 1 then
        error "the number of variables does not match the number of polynomials";

    if not f(r) == 0 then
        error "the field element is not a zero of the function";

    F := (sub(x,fracR)-sub(r,fracR))^m * g/f;

    makeAntidiagonalUnstableForm(kk, F(r), m)
)