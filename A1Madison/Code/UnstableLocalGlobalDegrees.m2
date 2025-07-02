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
    	error "rational function is not reduced"
       
    -- Get the underlying ring and ensure it is a field
    kk := coefficientRing ring(q)
    q not isField kk then kk = toField kk
    
    -- Check whether the rational function has isolated zeros
    if dim ideal(f) > 0 then 
        error "rational function does not have isolated zeros"
    
    -- Check whether the number of variables matches the number of polynomials
    
    S := ring f --or q?
    if #(gens S) != 1 then
        error "the number of variables does not match the number of polynomials"
    
    -- If the field is CC, output the Grothendieck-Witt class of an identity matrix of the appropriate rank
    if instance(kk, ComplexField) then (
    	rankAlgebra := getGlobalAlgebraRank {q};
    	return makeGWuClass id_(CC^rankAlgebra);
        )

    -- If the field is RR, ask the user to run the computation over QQ instead and then base change to RR
    if instance(kk, RealField) then error "getGlobalUnstableA1Degree method does not work over the reals. Instead, define the polynomials over QQ to output an unstableGrothendieckWittClass. Then extract the matrix, base change it to RR, and run getSumDecomposition().";    
    -- Create internal rings/matrices
    -- Initialize a polynomial ring in X_i's and Y_i's in which to compute the Bezoutian
--    X := local X;
--    Y := local Y;
--    R := kk(monoid[X | Y]);

--    numD   := sub(f,X)*sub(g,Y) - sub(f,Y)*sub(g,X)
--    denomD := X - Y
--    qD     := numD/ denomD

    n := 1;

    -- Create internal rings/matrices
    
    -- Initialize a polynomial ring in X_i's and Y_i's in which to compute the Bezoutian
    X := local X;
    Y := local Y;
    R := kk(monoid[X_1..X_n | Y_1..Y_n]);

    -- Create an (n x n) matrix D to be populated by \Delta_{ij} from the Brazelton-McKean-Pauli paper
    D := "";
    try D = mutableMatrix id_((frac R)^n) else D = mutableMatrix id_(R^n);
    
    for i from 0 to n - 1 do (
	for j from 0 to n - 1 do (
	    -- Iterate through the entries of the matrix D and populate it as follows
            -- Create the list {Y_1,...,Y_(j-1), X_j,...,X_n}. Note that Macaulay2 is 0-indexed hence the difference in notation. 
	    targetList1 := apply(toList(Y_1..Y_j | X_(j+1)..X_n), i->i_R);
            -- Create the list {Y_1,...,Y_j, X_(j+1),...,X_n}. Note that Macaulay2 is 0-indexed hence the difference in notation. 
	    targetList2 := apply(toList (Y_1..Y_(j+1) | X_(j+2)..X_n), i->i_R);
            -- Suppose our endomorphisms are given in the variables x_1,...,x_n
    --        -- Map f_i(x_1,...,x_n) to f_i(Y_1,...,Y_(j-1), X_j,...,X_n) resp.
            -- Take the difference f_i(Y_1,...,Y_(j-1), X_j,...,X_n) - f_i(Y_1,...,Y_j, X_(j+1),...,X_n)
            numeratorD := (map(R, S, targetList1)) (({f}_i) * ({g}_i)) - (map(R, S, targetList2)) (({f}_i) * ({g}_i)); 
            -- Divide this by X_j - Y_j. Note that Macaulay2 is 0-indexed hence the difference in notation. 
	    D_(i,j) = numeratorD / ((X_(j+1))_R - (Y_(j+1))_R); 
	    ); 
        );
    	 
    -- Set up the local variables bezDet and bezDetR
    bezDet:="";
    bezDetR:="";
      
    -- The determinant of D is interpreted as an element of Frac(k[x_1..x_n]), so we can try to lift it to k[x_1..x_n]           
    if liftable(det D, R) then bezDetR = lift(det D, R);
    
    -- In some computations, applying lift(-,R) doesn't work, so we instead lift the numerator and
    -- then divide by a lift of the denominator (which will be a scalar) to the coefficient ring kk
    if not liftable(det D, R) then (
	bezDet = lift(numerator det D, R) / lift(denominator det D, coefficientRing R);
    	bezDetR = lift(bezDet, R);
	);

    -- Define formal variables X_i and Y_i that replace x_i
    RX := kk[X_1..X_n]; 
    RY := kk[Y_1..Y_n];

    -- mapxtoX replaces all instances of x_i with X_i; mapxtoY replaces all instances of y_i with Y_i
    mapxtoX := map(RX,S,toList(X_1..X_n)); 
    mapxtoY := map(RY,S,toList(Y_1..Y_n));

    -- Compute the standard basis of kk[X_1,...,X_n]/(f_1,...,f_n)
    standBasisX := basis (RX/(ideal (leadTerm (mapxtoX ideal {f})))); 
    standBasisY := basis (RY/(ideal (leadTerm (mapxtoY ideal {f})))); 
      
)
