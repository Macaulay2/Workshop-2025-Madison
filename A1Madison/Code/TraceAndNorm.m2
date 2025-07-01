restart

    AlgebraicTrace = (C,I,a,FF) -> (
	B:=basis(C/I);
	r:=degree I;
	Q:=(a)*(transpose B)*B;
	toVector := q -> last coefficients(q,Monomials=>B);
	fieldTrace := q -> (M:=toVector(q*B_(0,0));i=1;while i<r do
	    (M=M|(toVector (q*B_(0,i))) ; i=i+1); trace M);
	lift(fieldTrace a,FF))
    
 AlgebraicNorm = (C,I,a,FF) -> (
	B:=basis(C/I);
	r:=degree I;
	Q:=(a)*(transpose B)*B;
	toVector := q -> last coefficients(q,Monomials=>B);
	fieldNorm := q -> (M:=toVector(q*B_(0,0));i=1;while i<r do
	    (M=M|(toVector (q*B_(0,i))) ; i=i+1); det M);
	lift(fieldNorm a,FF))
       



    
MatrixM = (C,I,a,FF) -> (
	B:=basis(C/I);
	r:=degree I;
	Q:=(a)*(transpose B)*B;
	toVector := q -> last coefficients(q,Monomials=>B);
	Matrep := q -> (M:=toVector(q*B_(0,0));i=1;while i<r do
	    (M=M|(toVector (q*B_(0,i))) ; i=i+1); M);
	lift(Matrep a,FF))
