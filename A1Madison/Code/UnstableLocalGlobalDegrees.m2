-------------------------------------
-- unstable A1-Brouwer degree methods (task 3.12 from Overleaf)
-------------------------------------

-- Input: A rational function f/g where gcd(f,g) = 1

-- Output: A pair (M,a) where M is a matrix and a is a scalar

getGlobalUnstableA1Degree = method()
getGlobalUnstableA1Degree RingElement := (Matrix,Number) => q -> (
    
    -- Extract numerator f from q
    f := numerator(q);
    
    -- Extract numerator g from q
    g := denominator(q);
    
    -- Check q = f/g is reduced
    if gcd(f,g) != 1 then
    	error "rational function is not reduced";
       
    -- Get the underlying ring and ensure it is a field
    kk := coefficientRing ring(q)
    q not isField kk then kk = toField kk;
    
    -- Check whether the rational function has isolated zeros
    if dim ideal(f) > 0 then 
        error "rational function does not have isolated zeros";
    
    -- Check whether the number of variables matches the number of polynomials
    
    S := ring f;
    if #(gens S) != 1 then
        error "the number of variables does not match the number of polynomials";
    
    -- If the field is CC, output the Grothendieck-Witt class of an identity matrix of the appropriate rank
    if instance(kk, ComplexField) then (
    	rankAlgebra := getGlobalAlgebraRank list{q};
    	return makeGWuClass id_(CC^rankAlgebra);
        );
	
    
)
