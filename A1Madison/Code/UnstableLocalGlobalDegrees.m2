-------------------------------------------------------------
-- unstable A1-Brouwer degree methods (task 3.12 from Overleaf)
-------------------------------------------------------------

-- Input: A rational function q = f/g

-- Output: A pair (M,a) where M is a matrix and a is a scalar (the determinant of M)

getGlobalUnstableA1Degree = method()
getGlobalUnstableA1Degree RingElement := (Matrix,Number) => q -> (

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
    
    fX = sub(f,{x => X});
    fY = sub(f,{x => Y});
    gX = sub(g,{x => X});
    gY = sub(g,{x => Y});
    
    D = lift((fX * gY - fY * gX)/(X-Y),R');
    
    m := degree (X,D);
    n := degree (Y,D);
        
    B := mutableMatrix id_(QQ^(m+1));  
    
    for i from 0 to m do(
	for j from 0 to n do
    	B_(i,j) = coefficient(X^i*Y^j,D)
	);
    
     makeGWuClass matrix B
     )
