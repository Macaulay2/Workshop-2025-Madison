---------------------------------------------------------------
-- unstable A1-Brouwer degree methods (task 3.12 from Overleaf)
---------------------------------------------------------------

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
    
     makeGWClass matrix B
     )
    
    
    
    
   --   L1 := toList (0..m)
   --   L2 := toList (0..n)
   --  blist := apply (L1 ** L2,(i,j) -> (X^i)*(Y^j))


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
    
    m := degree (X,Bezzz)
    n := degree (Y,Bezzz)
   --   L1 := toList (0..m)
   --   L2 := toList (0..n)
   --  blist := apply (L1 ** L2,(i,j) -> (X^i)*(Y^j))
    
    B := mutableMatrix id_(QQ^(m+1));  
    -- number of possible monomial terms in the Bezoutian (even if they are zero)
    
    for i from 0 to m do(
	for j from 0 to n do
    	B_(i,j) = coefficient(X^i*Y^j,Bezzz)
	)
Bezzz = 3*X^2 + 3*X*Y + 3*Y^2 + X + Y 
    
Bez = X^4*Y^4 - 5*X^4*Y^2 - 16*X^3*Y^3 - 5*X^2*Y^4 + 7*X^4*Y + 39*X^3*Y^2 + 39*X^2*Y^3 + 7*X*Y^4 + X^4 - 29*X^3*Y - 84*X^2*Y^2 - 29*X*Y^3 + Y^4 - 14*X^3 + 63*X^2*Y + 63*X*Y^2 - 14*Y^3 + 11*X^2 - 63*X*Y + 11*Y^2 + 38*X + 38*Y - 68     
     
--    (M,C) := coefficients bez
--    B := mutableMatrix id_(QQ^b);
--    for i from 0 to 1 do (
--        for j from 0 to 1 do
--           if i == 0 then B_(i,j) = coefficient(M_(i,3-j),bez) else B_(i,j) = coefficient(M_(1-i,1-j),bez);
--        );    
--    (makeGWClass matrix B, det matrix B) -- check makeGWuClass
    
   
   
        
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

