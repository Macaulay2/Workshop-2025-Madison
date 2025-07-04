restart    

--We define a method to represent multiplication by an element in a finite dimensional K-algebra as a matrix indexed by basis elements of its base field

MultiplicationMatrix=method()

--This applies the method by accepting a K-algebra C and an element a as inputs

MultiplicationMatrix(Ring,Thing):= (C,a) -> (
	B:=basis(C);
	r:=degree C;
	Q:=(a)*(transpose B)*B;
	toVector := q -> last coefficients(q,Monomials=>B);
	Matrep := q -> (M:=toVector(q*B_(0,0));i:=1;while i<r do
	    (M=M|(toVector (q*B_(0,i))) ; i=i+1); M);
	lift(Matrep a, coefficientRing C))

--This applies the method by accepting a polynomial ring C, an ideal I and an element a as input to find matrix representation of multiplication by the element over the corresponding quotient ring 

MultiplicationMatrix(Ring,Ideal,Thing):= (S,I,b) -> (
	B:=basis(S/I);
	r:=degree I;
	Q:=(b)*(transpose B)*B;
	toVector := q -> last coefficients(q,Monomials=>B);
	Matrep := q -> (M:=toVector(q*B_(0,0));i:=1;while i<r do
	    (M=M|(toVector (q*B_(0,i))) ; i=i+1); M);
	lift(Matrep b, coefficientRing S))

--We use the implemented matrix representation to calclate the algebraic trace
    
AlgebraicTrace=method()

AlgebraicTrace(Ring,Thing) := (C,a) -> (
	M:=MultiplicationMatrix(C,a);
	trace M)
    
AlgebraicTrace(Ring,Ideal,Thing) := (S,I,b) -> (
	M:=MultiplicationMatrix(S,I,b);
	trace M)

--We use the implemented matrix representation to calclate the algebraic norm 
    
AlgebraicNorm=method()

AlgebraicNorm(Ring,Thing) := (C,a) -> (
	M:=MultiplicationMatrix(C,a);
	det M)

 AlgebraicNorm(Ring,Ideal,Thing) := (S,I,b) -> (
	M:=MultiplicationMatrix(S,I,b);
	det M)
