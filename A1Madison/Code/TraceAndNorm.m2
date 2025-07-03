restart    
MultiplicationMatrix = (C,a) -> (
	B:=basis(C);
	r:=degree C;
	Q:=(a)*(transpose B)*B;
	toVector := q -> last coefficients(q,Monomials=>B);
	Matrep := q -> (M:=toVector(q*B_(0,0));i=1;while i<r do
	    (M=M|(toVector (q*B_(0,i))) ; i=i+1); M);
	lift(Matrep a, coefficientRing C))

 AlgebraicTrace = (C,a) -> (
	M=MultiplicationMatrix(C,a);
	trace M)
    
 AlgebraicNorm = (C,a) -> (
	M=MultiplicationMatrix(C,a);
	det M)
