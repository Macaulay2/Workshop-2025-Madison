document {
	Key => {getMultiplicationMatrix, (getMultiplicationMatrix, Ring, Thing), (getMultiplicationMatrix, Ring, Ideal, Thing) },
	Headline => "Computes the matrix over a ",TEX///$\mathbb{K}$///," -basis for multiplication by an element in a finite dimensional ",TEX///$\mathbb{K}$///,"-algebra",
	Usage => "getMultiplicationMatrix(C,a)
              getMultiplicationMatrix(S,I,b)",
	Inputs => {
	    Ring=> "C" => {"a finite dimensional ",TEX///$\mathbb{K}$///,"-algebra"},
	    Thing=> "a" => {"an element in C"},
	    Ring=> "S" => {"a polynomial ring"},
	    Ideal=> "I" => {"an ideal in the polynomial ring S"},
	    Thing=> "b" => {"an element in S/I"}
	    },
	Outputs => {
	    Matrix => {"the matrix representation over a basis for multiplication by the given element"}
	    },
	PARA {"For an algebra C over a field ",TEX///$\mathbb{K}$///," or a polynomial ring S and an ideal I, this function generates a matrix with entries in the ",TEX///$\mathbb{K}$///," representing multiplication by the user prescribed element in C or S/I respectively"},
	EXAMPLE lines ///
		 L = QQ[x]/(x^6+x^5+x^4+x^3+x^2+x+1)
         F = toField L
		 N=getMultiplicationMatrix(F[a,b,c], ideal(a^2,b^2,c^2),1+a*b+b*c+c*a)
	 	 ///,
	EXAMPLE lines ///
		QQ[x,y]
        L = QQ[x,y]/(x^2+y^2+1)
		F = frac L	
		A=getMultiplicationMatrix(F[z], ideal(z^2+1), 1+y*x^2*z)
		///,
	SeeAlso => {"getTrace", "getNorm"}
     	}

document {
	Key => {getTrace, (getTrace, Ring, Thing), (getTrace, Ring, Ideal, Thing) },
	Headline => "Computes the algebraic trace over ",TEX///$\mathbb{K}$///," for an element in a finite dimensional ",TEX///$\mathbb{K}$///,"-algebra",
	Usage => "getTrace(C,a)
              getTrace(S,I,b)",
	Inputs => {
	    Ring=> "C" => {"a finite dimensional ",TEX///$\mathbb{K}$///,"-algebra"},
	    Thing=> "a" => {"an element in C"},
	    Ring=> "S" => {"a polynomial ring"},
	    Ideal=> "I" => {"an ideal in the polynomial ring S"},
	    Thing=> "b" => {"an element in S/I"}
	    },
	Outputs => {
	    RingElement => {"the algebraic trace over ",TEX///$\mathbb{K}$///," for an element in the algebra"}
	    },
	PARA {"For an element in an algebra C over a field ",TEX///$\mathbb{K}$///," or a polynomial ring S and an ideal I, this function computes the algebraic trace over ",TEX///$\mathbb{K}$///,""},
	EXAMPLE lines ///
		 L = QQ[x]/(x^6+x^5+x^4+x^3+x^2+x+1)
         F = toField L
		 N=getTrace(F[a,b,c], ideal(a^2,b^2,c^2),1+a*b+b*c+c*a)
	 	 ///,
	EXAMPLE lines ///
		QQ[x,y]
        L = QQ[x,y]/(x^2+y^2+1)
		F = frac L	
		A=getTrace(F[z], ideal(z^2+1), 1+y*x^2*z)
		///,
	SeeAlso => {"getMultiplicationMatrix", "getNorm"}
     	}

		document {
	Key => {getNorm, (getNorm, Ring, Thing), (getNorm, Ring, Ideal, Thing) },
	Headline => "Computes the algebraic norm over ",TEX///$\mathbb{K}$///," for an element in a finite dimensional ",TEX///$\mathbb{K}$///,"-algebra",
	Usage => "getNorm(C,a)
              getNorm(S,I,b)",
	Inputs => {
	    Ring=> "C" => {"a finite dimensional ",TEX///$\mathbb{K}$///,"-algebra"},
	    Thing=> "a" => {"an element in C"},
	    Ring=> "S" => {"a polynomial ring"},
	    Ideal=> "I" => {"an ideal in the polynomial ring S"},
	    Thing=> "b" => {"an element in S/I"}
	    },
	Outputs => {
	    RingElement => {"the algebraic norm over ",TEX///$\mathbb{K}$///," for an element in the algebra"}
	    },
	PARA {"For an element in an algebra C over a field ",TEX///$\mathbb{K}$///," or a polynomial ring S and an ideal I, this function computes the algebraic norm over ",TEX///$\mathbb{K}$///,""},
	EXAMPLE lines ///
		 L = QQ[x]/(x^6+x^5+x^4+x^3+x^2+x+1)
         F = toField L
		 N=getNorm(F[a,b,c], ideal(a^2,b^2,c^2),1+a*b+b*c+c*a)
	 	 ///,
	EXAMPLE lines ///
		QQ[x,y]
        L = QQ[x,y]/(x^2+y^2+1)
		F = frac L	
		A=getNorm(F[z], ideal(z^2+1), 1+y*x^2*z)
		///,
	SeeAlso => {"getMultiplicationMatrix", "getTrace"}
     	}
