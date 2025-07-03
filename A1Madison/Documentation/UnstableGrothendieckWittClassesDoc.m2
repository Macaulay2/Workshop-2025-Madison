document{
    Key => {UnstableGrothendieckWittClass, (net, UnstableGrothendieckWittClass), (texMath, UnstableGrothendieckWittClass)},
    Headline => "a new type, intended to capture the isomorphism class of an element of the unstable Grothendieck-Witt group of a base field",
    PARA {"An ", TT "UnstableGrothendieckWittClass" ," object is a type of ", TO2(HashTable, "HashTable"), " encoding the data of a ", TO2(GrothendieckWittClass, "GrothendieckWittClass"), " over a field ", TEX///$k$///, " and a scalar ", TEX///$b$///, " over the same field."},
	
    PARA{"An ", TT"UnstableGrothendieckWittClass", " object can be built in one of several way using the ", TO2(makeGWuClass, "makeGWuClass"), " method."},

	PARA{"We can construct an ", TT"UnstableGrothendieckWittClass", " over a field ", TEX///$k$///, " object by applying ", TT"makeGWuClass", " to a pair ", TT"(M,a)", " where ", TT"M", " is either a symmetric matrix of full rank or a ", TO2(GrothendieckWittClass, "GrothendieckWittClass"), " and ", TT"a", " is a scalar which agrees with the determinant of ", TT"M", " in the quotient ", TEX///$k^{\times}/(k^{\times})^{2}$///, "."},
	EXAMPLE lines///
	M = matrix(QQ, {{0,1},{1,0}})
	alpha = makeGWuClass(M, 1)
	class alpha
	///,
    EXAMPLE lines///
    beta0 = makeGWClass matrix(QQ, {{0,1},{1,0}})
    beta = makeGWuClass(beta0, 1)
	class beta
    ///,
	PARA{"Alternatively, we can construct an ", TT"UnstableGrothendieckWittClass", " object by applying ", TT"makeGWuClass", " to", TT"M", " where ", TT"M", " is either a symmetric matrix of full rank or a ", TO2(GrothendieckWittClass, "GrothendieckWittClass"), ", and the scalar ", TT"a", "is assumed to be the determinant of (a matrix representative of)", TT"M", "."},
	EXAMPLE lines///
	N = matrix(QQ, {{0,1},{1,0}})
	gamma = makeGWuClass(N)
	class gamma
	///,
    EXAMPLE lines///
    delta0 = makeGWClass matrix(QQ, {{0,1},{1,0}})
    delta = makeGWuClass(delta0)
	class delta
    ///,
    PARA{"The underlying matrix representative, base field, and scalar can be recovered by the commands ", TO2(getMatrix, "getMatrix"), ", ", TO2(getBaseField, "getBaseField"), ", and ", TO2(getScalar, "getScalar"), ", respectively."},
    EXAMPLE lines///
    getMatrix alpha
    getBaseField alpha
	getScalar alpha
    ///,
    SeeAlso => {"makeGWuClass", "getBaseField", "getMatrix", "getScalar"},
    }

document {
    Key => {makeGWuClass, (makeGWuClass, Matrix)},
	Headline => "constructor for unstable Grothendieck-Witt classes",
	Usage => "makeGWuClass(alpha)
			  makeGWuClass(alpha, a)
			  makeGWuClass(M)
			  makeGWuClass(M, a)",
	Inputs => {
		GrothendieckWittClass => "alpha" => {"a ", TO2(GrothendieckWittClass, "GrothendieckWittClass"), "."},
	    Matrix => "M" => {"a non-singular symmetric matrix defined over an arbitrary field of characteristic not 2"},
		Number => "a" => {"a nonzero element", TEX///$a\in k^{\times}$///, "."}
	    },
	Outputs => {
	    UnstableGrothendieckWittClass => {"the unstable Grothendieck-Witt class represented by the input data."}
	    },
	PARA {"The unstable Grothendieck-Witt class represented by the data ", TT"(M,a)", " or ", TT"(alpha, a)", " where ", TT"a", " is assumed to be the determinant of (a matrix representative of)", TT"M", " or ", TT"alpha", " if not otherwise specified."},
	EXAMPLE lines///
	M = matrix(QQ, {{0,1},{1,0}})
	alpha = makeGWuClass(M, 1)
	class alpha
	///,
    EXAMPLE lines///
    beta0 = makeGWClass matrix(QQ, {{0,1},{1,0}})
    beta = makeGWuClass(beta0, 1)
	class beta
    ///,
	EXAMPLE lines///
	N = matrix(QQ, {{0,1},{1,0}})
	gamma = makeGWuClass(N)
	class gamma
	///,
    EXAMPLE lines///
    delta0 = makeGWClass matrix(QQ, {{0,1},{1,0}})
    delta = makeGWuClass(delta0)
	class delta
    ///,
	PARA {"Over the complex numbers, real numbers, rational numbers, or finite fields, the constructor ", TT"makeGWuClass", " verifies that the scalar ", TT"a", ", if provided, agrees with the determinant of ", TT"M", " as an element of ", TEX///$k^{\times}/(k^{\times})^{2}$///, ". Over arbitrary fields, the constructor will permit the construction of an unstable Grothendieck-Witt class with any nonzero scalar and warn the user to verify that the scalar agrees with the determinant of the matrix representative of the stable part up to squares."},
	EXAMPLE lines///
	N = matrix(QQ, {{0,1},{1,0}})
	makeGWuClass(N, 1)
	makeGWuClass(N, -1)
	///,
    EXAMPLE lines///
    R = QQ[x]/(x^2 + 1);
	K = toField R
	P = matrix(K, {{1,0}, {0,x}})
	makeGWuClass(P, x^3)
    ///,

	SeeAlso => {"GrothendieckWittClass", "getMatrix", "getBaseField", "getScalar"}
        }

document {
    Key => {getScalar, (getScalar, UnstableGrothendieckWittClass)},
	Headline => "the underlying scalar of an unstable Grothendieck-Witt class",
	Usage => "getScalar beta",
	Inputs => {
	    UnstableGrothendieckWittClass => "beta" => {"an unstable Grothendieck-Witt class over a field of characteristic not 2"}
	    },
	Outputs => {
	    Ring => {"the underlying scalar of the unstable Grothendieck-Witt class ", TT "beta"}
	    },
	PARA {"Given an unstable Grothendieck-Witt class, ", TT "beta", ", this command outputs the scalar underlying the class."},
	EXAMPLE lines ///
		 beta = makeGWuClass matrix(QQ, {{0,2},{2,0}});
		 getScalar beta
	 	 ///,
    SeeAlso => {"UnstableGrothendieckWittClass", "makeGWuClass"}
        }

document {
    Key => {getMatrix, (getMatrix, GrothendieckWittClass), (getMatrix, UnstableGrothendieckWittClass)},
	Headline => "the underlying matrix of a Grothendieck-Witt class or unstable Grothendieck-Witt class",
	Usage => "getMatrix beta",
	Inputs => {
	    GrothendieckWittClass => "beta" => {"the isomorphism class of a non-degenerate symmetric bilinear form over a field of characteristic not 2"},
		UnstableGrothendieckWittClass => "beta" => {"an unstable Grothendieck-Witt class over a field of characteristic not 2"}
	    },
	Outputs => {
	    Ring => {"the underlying matrix of the (unstable) Grothendieck-Witt class ", TT "beta"}
	    },
	PARA {"Given the isomorphism class of a symmetric bilinear form, ", TT "beta", ", this command outputs the underlying matrix of the form."},
	EXAMPLE lines///
	M = matrix(QQ, {{0,1},{1,0}})
	getMatrix makeGWClass
	getMatrix makeGWuClass(M, -2)
	///,
    SeeAlso => {"GrothendieckWittClass", "makeGWClass", "UnstableGrothendieckWittClass", "makeGWuClass"}
        }

document {
    Key => {getBaseField, (getBaseField, GrothendieckWittClass), (getBaseField, UnstableGrothendieckWittClass)},
	Headline => "the base field of a Grothendieck-Witt class or unstable Grothendieck-Witt class",
	Usage => "getBaseField beta",
	Inputs => {
	    GrothendieckWittClass => "beta" => {"the isomorphism class of a non-degenerate symmetric bilinear form over a field of characteristic not 2"},
		UnstableGrothendieckWittClass => "beta" => {"an unstable Grothendieck-Witt class over a field of characteristic not 2"}
	    },
	Outputs => {
	    Ring => {"the base field of the (unstable) Grothendieck-Witt class ", TT "beta"}
	    },
	PARA {"Given the isomorphism class of a symmetric bilinear form, ", TT "beta", ", this command outputs the base field of the form."},
	EXAMPLE lines///
	M = matrix(QQ, {{0,1},{1,0}})
	getBaseField makeGWClass M
	getBaseField makeGWuClass(M, -2)
	///,
    SeeAlso => {"GrothendieckWittClass", "UnstableGrothendieckWittClass", "makeGWuClass"}
        }

document {
    Key => {addGWuDivisorial, (addGWuDivisorial, List, List)},
	Headline => "the divisorial sum of local degrees of a rational function",
	Usage => "addGWuDivisorial(L1, L2)",
	Inputs => {
	    List => "L1" => {"a list of unstable Grothendieck-Witt classes representing unstable local degrees"},
		List => "L2" => {"a list of elements of the base field corresponding to the roots at which the unstable local degrees of ", TT"L1", " are computed."}
	    },
	Outputs => {
	    UnstableGrothendieckWittClass => {"the divisorial sum of the unstable Grothendieck-Witt classes in ", TT "L1", " with respect to the divisor determined by the points in ", TT"L2", "."}
	    },
	PARA {"Let ", TEX///$\frac{f(x)}{g(x)}:\mathbb{P}^{1}_{k}\to\mathbb{P}^{1}_{k}$///, " be a pointed rational function with zeroes ", TEX///$\{r_{1},\dots,r_{n}\}$///, " and ", TEX///$\{\beta_{1},\dots,\beta_{n}\}$///, " the unstable local ", TEX///$\mathbb{A}^{1}$///, "-degrees at the ", TEX///$r_{i}$///, ". The unstable global ", TEX///$\mathbb{A}^{1}$///, "-degree of the rational function is not computed as the ", TO2(addGWu, "addGWu"), " of the local unstable degrees, but as the divisorial sum [IMS+24, Thm. 1.2]."},
	PARA {"The following example computes the divisorial sum of the rational function ", TEX///$\frac{x^{2}+x-2}{3x+5}$///, " over ", TEX///$\mathbb{Q}$/// , " where the lists of unstable Grothendieck-Witt classes are given by ", TEX///$\{(\langle \frac{1}{3}\rangle, \frac{1}{3}), (\langle \frac{8}{3}\rangle, \frac{8}{3})\}$///, " and ", TEX///$\{-2, 1\}$///, "." },
	EXAMPLE lines///
	M1 = matrix(QQ, {{1/3}})
	alpha = makeGWuClass(M1)
	M2 = matrix(QQ, {{8/3}})
	beta = makeGWuClass(M2)
	addGWuDivisorial({alpha, beta}, {-2, 1})
	///,
	 PARA{EM "Citations:"},
    UL{
	{"[IMS+24] Igieobo, et. al., ", EM "Motivic configurations on the line,", TT"arXiv: 2411.15347", "."},
    },
    SeeAlso => {"UnstableGrothendieckWittClass"} -- to add local and global degree functions
        }