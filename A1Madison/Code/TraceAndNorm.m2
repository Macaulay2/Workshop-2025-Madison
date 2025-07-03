restart    

MultiplicationMatrix=method()


MultiplicationMatrix(Ring,Thing):= (C,a) -> (
	B:=basis(C);
	r:=degree C;
	Q:=(a)*(transpose B)*B;
	toVector := q -> last coefficients(q,Monomials=>B);
	Matrep := q -> (M:=toVector(q*B_(0,0));i=1;while i<r do
	    (M=M|(toVector (q*B_(0,i))) ; i=i+1); M);
	lift(Matrep a, coefficientRing C))


MultiplicationMatrix(Ring,Ideal,Thing):= (C,I,a) -> (
	B:=basis(C/I);
	r:=degree I;
	Q:=(a)*(transpose B)*B;
	toVector := q -> last coefficients(q,Monomials=>B);
	Matrep := q -> (M:=toVector(q*B_(0,0));i=1;while i<r do
	    (M=M|(toVector (q*B_(0,i))) ; i=i+1); M);
	lift(Matrep a, coefficientRing C))

AlgebraicTrace=method()

AlgebraicTrace(Ring,Thing) := (C,a) -> (
	M=MultiplicationMatrix(C,a);
	trace M)
    
AlgebraicTrace(Ring,Ideal,Thing) := (C,I,a) -> (
	M=MultiplicationMatrix(C,I,a);
	trace M)

AlgebraicNorm=method()

AlgebraicNorm(Ring,Thing) := (C,a) -> (
	M=MultiplicationMatrix(C,a);
	det M)

 AlgebraicNorm(Ring,Ideal,Thing) := (C,I,a) -> (
	M=MultiplicationMatrix(C,I,a);
	det M)
