doc ///
    Key
        getSignature
        (getSignature, GrothendieckWittClass)
    Headline
        computes the signature of a symmetric bilinear form over the real numbers or rational numbers
    Usage
        getSignature(beta)
    Inputs
        beta: GrothendieckWittClass
            a symmetric bilinear form defined over $\mathbb{Q}$ or $\mathbb{R}$
    Outputs
        :ZZ
            the signature of the symmetric bilinear form $\beta$
    Description
        Text
            Given a symmetric bilinear form, after diagonalizing it, we can consider the number of positive entries minus the number of negative entries appearing along the diagonal. This is the @ITALIC("signature")@ of a symmetric bilinear form, and is one of the primary invariants used to classify forms. For more information, see @TO2(isIsomorphicForm, "isIsomorphicForm")@. 
        Example
            M = matrix(RR, {{0,0,1},{0,1,0},{1,0,0}});
            beta = makeGWClass M;
            getSignature beta
    SeeAlso
        isIsomorphicForm
        getHilbertSymbol
        getSumDecomposition
        getSumDecompositionString
///

doc ///  
    Key 
        getIntegralDiscriminant
        (getIntegralDiscriminant, GrothendieckWittClass)
    Headline 
        computes the integral discriminant for a rational symmetric bilinear form
    Usage 
        getIntegralDiscriminant(beta)
    Inputs
        beta: GrothendieckWittClass
            denoted by $\beta \in \text{GW}(\mathbb{Q})$
    Outputs
        :ZZ
            an integral square class representative of $\text{disc}(\beta)$
    Description
        Example 
            beta = makeGWClass matrix(QQ, {{1,4,7},{4,3,-1},{7,-1,5}});
            getIntegralDiscriminant beta
            getDiagonalClass beta
    SeeAlso
        isIsomorphicForm
///


doc /// 
    Key
        getHasseWittInvariant
        (getHasseWittInvariant, GrothendieckWittClass, ZZ)
    Headline
        computes the Hasse-Witt invariant at a prime $p$ for the quadratic form of the Grothendieck-Witt class
    Usage
        getHasseWittInvariant(beta, p)
    Inputs
        beta: GrothendieckWittClass
            denoted by $\beta \in \text{GW}(\mathbb{Q})$
        p: ZZ
            a prime number  
    Outputs
        :ZZ
            the Hasse-Witt invariant of $\beta$ at the prime $p$
    Description
        Text
            The Hasse-Witt invariant of a diagonal form $\langle a_1,\ldots,a_n\rangle$ over a field $K$ is defined to be the product $\prod_{i<j} \left(a_i,a_j\right)_p$ where $(-,-)_p$ is the @TO2(getHilbertSymbol, "Hilbert symbol")@. 
            
            The Hasse-Witt invariant of a form will be equal to 1 for almost all primes.  In particular, after diagonalizing a form $\beta \cong \left\langle a_1,\ldots,a_n\right\rangle$, the Hasse-Witt invariant at a prime $p$ will be 1 automatically if $p\nmid a_i$ for all $i$. Thus we only have to compute the Hasse-Witt invariant at @TO2(getRelevantPrimes, "primes dividing diagonal entries")@.
        Example
            beta = makeGWClass matrix(QQ, {{1,4,7},{4,3,-1},{7,-1,5}});
            getHasseWittInvariant(beta, 7)
    SeeAlso
        isIsomorphicForm
        getRelevantPrimes
///       




document{
    Key => {getRelevantPrimes, (getRelevantPrimes, GrothendieckWittClass)},
    Headline => "outputs a list containing all primes p where the Hasse-Witt invariant of a symmetric bilinear form is nontrivial",
    Usage => "getRelevantPrimes beta",
    Inputs => {
	GrothendieckWittClass => "beta" => {"denoted by ", TEX///$\beta\in\text{GW}(\mathbb{Q})$///},
	},
    Outputs => {
        List => {"a finite list of primes ", TEX///$(p_1,\ldots,p_r)$///, " containing all primes ",
	    TEX///$p$///, " where the ",TO2(getHasseWittInvariant,"Hasse-Witt invariant "), TEX///$\phi_p(\beta)$///," is nontrivial"},
	},
    PARA{"It is a classical result that the ", TO2(getHasseWittInvariant,"Hasse-Witt invariants"), " of a quadratic form are equal to 1 for all but finitely many primes (see e.g. [S73, IV Section 3.3]). As the Hasse-Witt invariants are computed as a product of ", TO2(getHilbertSymbol,"Hilbert symbols") , " of the pairwise entries appearing on a diagonalization of the symbol, it suffices to consider primes dividing diagonal entries."},
   EXAMPLE lines ///
   beta = makeDiagonalForm(QQ, (6,7,22));
   getRelevantPrimes beta
   ///,
   PARA{EM "Citations:"},
    UL{
	{"[S73] J.P. Serre, ", EM "A course in arithmetic,", " Springer-Verlag, 1973."},
    },
    SeeAlso => {"getHasseWittInvariant"}
}

document {
        Key => {getRank, (getRank, GrothendieckWittClass), (getRank, Matrix)},
        Headline => "calculates the rank of a symmetric bilinear form",
        Usage => "getRank beta",
        Inputs => {
            GrothendieckWittClass => "beta" => {"denoted by ",  TEX///$\beta \in  GW(\mathbb{Q})$///, " or a symmetric matrix ", TEX///$\beta$///}
            },
        Outputs => {
            ZZ => {"the rank of the form ", TEX///$\beta$///}
           },
        PARA {"This computes the rank of the form ", TEX///$\beta$/// },
        EXAMPLE lines ///
                 beta = makeDiagonalForm(QQ, (3,5,7,11))
                 getRank beta
		 ///,                
        EXAMPLE lines ///
                 M = matrix(QQ, {{1,4,7},{4,3,-1},{7,-1,5}})
                 getRank M
                 ///,
    SeeAlso => {"isIsomorphicForm", "getSignature"}
}
