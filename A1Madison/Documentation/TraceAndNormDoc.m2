document {
	Key => {MultiplicationMatrix, (MultiplicationMatrix, Ring, Thing), (MultiplicationMatrix, Ring, Ideal, Thing) )},
	Headline => "Computes the matrix over a ",TEX///$\mathbb{K}$///," -basis for multiplication by an element in a finite dimensional ",TEX///$\mathbb{K}$///,"-algebra",
	Usage => "MultiplicationMatrix(C,a)
              MultiplicationMatrix(,I,a)",
	Inputs => {
	    Ring=> "C" => {"a finite dimensional K-algebra"},
	    Thing=> "a" => {"an element in C"},
	    Ring=> "S" => {"a polynomial ring"},
	    Ideal=> "I" => {"an ideal in the polynomial ring S"},
	    Thing=> "b" => {"an element in S/I"}
	    },
	Outputs => {
	    Matrix => {"the matrix representation over a basis for multiplication by the given element"}
	    },
	PARA {"Given an algebra over a field or a polynomial ring and an ideal contained within it, this function generates a matrix with entries in the field representing multiplication by the user generated element"},
	EXAMPLE lines ///
		 M = matrix(GF(17), {{7,9},{9,6}});
		 diagonalizeViaCongruence M
	 	 ///,
	SeeAlso => {"AlgebraicTrace", "AlgebraicNorm"}
     	}
