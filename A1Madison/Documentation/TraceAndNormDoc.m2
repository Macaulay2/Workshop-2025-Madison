document {
	Key => {multiplicationMatrix, (multiplicationMatrix, Ring, Thing), (multiplicationMatrix, Ring, Ideal, Thing) },
	Headline => "Computes the matrix over a ",TEX///$\mathbb{K}$///," -basis for multiplication by an element in a finite dimensional ",TEX///$\mathbb{K}$///,"-algebra",
	Usage => "multiplicationMatrix(C,a)
              multiplicationMatrix(S,I,b)",
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
		 N=multiplicationMatrix(F[a,b,c], ideal(a^2,b^2,c^2),1+a*b+b*c+c*a)
	 	 ///,
	EXAMPLE lines ///
		QQ[x,y]
        L = QQ[x,y]/(x^2+y^2+1)
		F = frac L	
		A=multiplicationMatrix(F[z], ideal(z^2+1), 1+y*x^2*z)
		///,
	SeeAlso => {"algebraicTrace", "algebraicNorm"}
     	}

document {
	Key => {algebraicTrace, (algebraicTrace, Ring, Thing), (algebraicTrace, Ring, Ideal, Thing) },
	Headline => "Computes the algebraic trace over ",TEX///$\mathbb{K}$///," for an element in a finite dimensional ",TEX///$\mathbb{K}$///,"-algebra",
	Usage => "algebraicTrace(C,a)
              algebraicTrace(S,I,b)",
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
		 N=algebraicTrace(F[a,b,c], ideal(a^2,b^2,c^2),1+a*b+b*c+c*a)
	 	 ///,
	EXAMPLE lines ///
		QQ[x,y]
        L = QQ[x,y]/(x^2+y^2+1)
		F = frac L	
		A=algebraicTrace(F[z], ideal(z^2+1), 1+y*x^2*z)
		///,
	SeeAlso => {"multiplicationMatrix", "algebraicNorm"}
     	}

		document {
	Key => {algebraicNorm, (algebraicNorm, Ring, Thing), (algebraicNorm, Ring, Ideal, Thing) },
	Headline => "Computes the algebraic norm over ",TEX///$\mathbb{K}$///," for an element in a finite dimensional ",TEX///$\mathbb{K}$///,"-algebra",
	Usage => "algebraicNorm(C,a)
              algebraicNorm(S,I,b)",
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
		 N=algebraicNorm(F[a,b,c], ideal(a^2,b^2,c^2),1+a*b+b*c+c*a)
	 	 ///,
	EXAMPLE lines ///
		QQ[x,y]
        L = QQ[x,y]/(x^2+y^2+1)
		F = frac L	
		A=algebraicNorm(F[z], ideal(z^2+1), 1+y*x^2*z)
		///,
	SeeAlso => {"multiplicationMatrix", "algebraicTrace"}
     	}
