-------------------------------------
-- unstable A1-Brouwer degree methods (task 3.12 from Overleaf)
-------------------------------------

-- Input: A rational function f/g where gcd(f,g) = 1

-- Output: A pair (M,a) where M is a matrix and a is a scalar

getGlobalUnstableA1Degree = method()
getGlobalUnstableA1Degree RingElement := (Matrix,Number) => q -> (
    
    -- Extract numerator f from q
    f := numerator(q)
    
    -- Extract numerator g from q
    g := denominator(q)
    
    -- Check q = f/g is reduced
    if gcd(f,g) != 1 then
    	error "rational function is not reduced";
       
    -- Get the underlying ring and ensure it is a field
    kk := coefficientRing ring(f)
)
