\-------------------------------------
-- unstable A1-Brouwer degree methods (task 3.12 from Overleaf)
-------------------------------------

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
        )

    -- If the field is RR, ask the user to run the computation over QQ instead and then base change to RR
    if instance(kk, RealField) then error "getGlobalUnstableA1Degree method does not work over the reals. Instead, define the polynomials over QQ to output an unstableGrothendieckWittClass. Then extract the matrix, base change it to RR, and run getSumDecomposition().";    

    n := 1;

    -- Create internal rings/matrices
    
    -- Initialize a polynomial ring in X and Y in which to compute the Bezoutian
    X := local X;
    Y := local Y;
    R' := kk[X,Y];
    
fX = sub(f,{x => X})

fY = sub(f,{x => Y})

gX = sub(g,{x => X})

gY = sub(g,{x => Y})

D = lift((fX * gY - fY * gX)/(X-Y),R')

    -- Create an (n x n) matrix D to be populated by \Delta_{ij} from the Brazelton-McKean-Pauli paper
    D := "";
    try D = mutableMatrix id_((frac R')^n) else D = mutableMatrix id_(R^n);
    
    for i from 0 to n - 1 do (
	for j from 0 to n - 1 do (
	    -- Iterate through the entries of the matrix D and populate it as follows
            -- Create the list {Y_1,...,Y_(j-1), X_j,...,X_n}. Note that Macaulay2 is 0-indexed hence the difference in notation. 
	    targetList1 := apply(toList(Y_1..Y_j | X_(j+1)..X_n), i->i_R);
            -- Create the list {Y_1,...,Y_j, X_(j+1),...,X_n}. Note that Macaulay2 is 0-indexed hence the difference in notation. 
	    targetList2 := apply(toList (Y_1..Y_(j+1) | X_(j+2)..X_n), i->i_R);
            -- Suppose our endomorphisms are given in the variables x_1,...,x_n
            -- Map f_i(x_1,...,x_n) to f_i(Y_1,...,Y_(j-1), X_j,...,X_n) resp.
            -- Take the difference f(X)g(Y) - f(Y)g(X)
            numeratorD := (map(R, S, targetList1)) (({f}_i) * ({g}_i)) - (map(R, S, targetList2)) (({f}_i) * ({g}_i)); 
            -- Divide this by X - Y. Note that Macaulay2 is 0-indexed hence the difference in notation. 
	    D_(i,j) = numeratorD / ((X_(j+1))_R - (Y_(j+1))_R); 
	    ); 
        );
    	 
    -- Set up the local variables bezDet and bezDetR
    bezDet:="";
    bezDetR:="";
      
    -- The determinant of D is interpreted as an element of Frac(k[x_1..x_n]), so we can try to lift it to k[x_1..x_n]           
    if liftable(D, R') then bezDetR = lift(D, R');
    
    -- In some computations, applying lift(-,R) doesn't work, so we instead lift the numerator and
    -- then divide by a lift of the denominator (which will be a scalar) to the coefficient ring kk
    if not liftable(D, R') then (
	bezDet = lift(numerator D, R') / lift(denominator D, coefficientRing R');
    	bezDetR = lift(bezDet, R');
	);

-----------------------------------
-- test example, need to generalize
-----------------------------------

R := QQ[x]

f = x^2 + x - 2
g = 3*x + 5

S := QQ[X,Y]

fX = sub(f,{x => X})

fY = sub(f,{x => Y})

gX = sub(g,{x => X})

gY = sub(g,{x => Y})

bez = sub((fX * gY - fY * gX)/(X-Y),S)


-- test

    -- Create the Bezoutian matrix B for the symmetric bilinear form by reading off the coefficients. 
    -- B is an (m x m) matrix. The coefficient B_(i,j) is the coefficient of the (ith basis vector x jth basis vector) in the tensor product.
    -- phi0 maps the coefficient to kk
    (M,C) := coefficients bez
    B := mutableMatrix id_(QQ^2);
    for i from 0 to 1 do (
        for j from 0 to 1 do
           if i == 0 then B_(i,j) = coefficient(M_(i,3-j),bez) else B_(i,j) = coefficient(M_(1-i,1-j),bez);
        );
    
    (makeGWClass matrix B, det matrix B) -- check makeGWuClass
        
----------------
-- end here
----------------    
    -- Define formal variables X and Y that replace x
    S := QQ[X,s]
    gXs := sub(gX,S)
    RX := S/(s*gXs - 1) 
    
    T := QQ[Y,t]
    gYt := sub(gY,T)
    RY := T/(t*gYt-1)

    -- mapxtoX replaces all instances of x with X; mapxtoY replaces all instances of x with Y
    mapxtoX := map(RX,R,{X})
    mapxtoY := map(RY,R,{Y})

    -- Compute the standard basis of kk[X,s]/(s*g(X) - 1, s*f(X))
    
    -- f/g
    f' = (mapxtoX f)*s
    
    standBasisX := basis (RX/(ideal (leadTerm (mapxtoX ideal f)))) 
    standBasisY := basis (RY/(ideal (leadTerm (mapxtoY ideal f)))); 
 
    -- Define an ideal (f_1(X),...,f_n(X))
    id1 := ideal apply(toList(0..n-1), i-> mapxtoX(Endo_i)); 
    -- Define an ideal (f_1(Y),...,f_n(Y))
    id2 := ideal apply(toList(0..n-1), i-> mapxtoY(Endo_i)); 

    -- Take the sum of ideals (f_1(X),...,f_n(X)) + (f_1(Y),...,f_n(Y)) in the ring kk[X_1..Y_n]
    promotedEndo := sub(id1, R) + sub(id2, R); 

    -- Take the sum of ideals (f_1(X),...,f_n(X)) + (f_1(Y),...,f_n(Y)) in the ring kk[X_1..Y_n]
    promotedEndo := sub(id1, R) + sub(id2, R); 

    -- Here we're using that (R/I) \otimes_R (R/J) = R/(I+J) in order to express Q(f) \otimes Q(f),
    -- where X's are the variables in first part and Y's are variables in second part
    Rquot := R/promotedEndo; 

    -- Move the standard bases to the quotient ring
    sBXProm := sub(standBasisX, Rquot); 
    sBYProm := sub(standBasisY, Rquot); 
    
    -- Reduce the bezDetR determinant subject to the ideal in the X's and Y's
    bezDetRed := bezDetR % promotedEndo;

    -- Define a ring map that takes the coefficients to the field kk instead of considering it as an element of the quotient ring
    phi0 := map(kk,Rquot, (toList ((2*n):0))); 

    -- m is the dimension of the basis for the algebra
    m := numColumns sBXProm;

    -- Create the Bezoutian matrix B for the symmetric bilinear form by reading off the coefficients. 
    -- B is an (m x m) matrix. The coefficient B_(i,j) is the coefficient of the (ith basis vector x jth basis vector) in the tensor product.
    -- phi0 maps the coefficient to kk
    B := mutableMatrix id_(kk^m);

    for i from 0 to m - 1 do (
        for j from 0 to m - 1 do
            B_(i,j) = phi0(coefficient((sBXProm_(0,i)**sBYProm_(0,j))_(0,0), bezDetRed));
        );
    makeGWClass matrix B
      
)

