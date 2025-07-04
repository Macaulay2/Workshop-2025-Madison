doc ///
    Key
        CpMackeyFunctor
    Headline
        The type of Cp Mackey functors
    Description
        Text
            This documents the @TO2(Type,"type")@ of $C_p$-Mackey functors. A $C_p$-Mackey functor is encoded as a @TO2(HashTable, "hash table")@ with keys through which we can pull various parts of the data.
        Text
            The keys available (internally) in the hash table are as follows, together with their external {\tt get} methods:
            
            @UL {
                (TT ".PrimeOrder", " yields the ", TO2(ZZ,"order"), " of the cyclic group over which we're working."),
                (TT ".Underlying", 
                " encodes the ", TO "module", " ",
                TEX"$M(C_p/e)$"),
                (TT ".Fixed",
                " encodes the ",TO "module", " ",
                TEX"$M(C_p/C_p)$"),
                (TT ".Res",
                " encodes the restriction ", TO2("Macaulay2Doc :: module homomorphisms","homomorphism")," ",
                TEX"$M(C_p/C_p) \\to M(C_p/e)$"),
                (TT ".Tr",
                " encodes the transfer ", TO2("Macaulay2Doc :: module homomorphisms","homomorphism")," ",
                TEX"$M(C_p/e) \\to M(C_p/C_p)$"),
                (TT ".Conj",
                " encodes the conjugation ", TO2("Macaulay2Doc :: module homomorphisms","homomorphism")," ",
                TEX"$M(C_p/e) \\to M(C_p/e)$"),
            }@

            {\bf Constructing new Mackey functors:}

            {\bf Operations with Mackey functor types:} One of the first operations to know about is @TO2((prune,CpMackeyFunctor),"pruning")@ a $C_p$-Mackey functor. There is more information at the prune command, but this basically takes the data of a Mackey functor and attempts to simplify it and make it more reader-friendly. Pruning is the only unary operation on $C_p$-Mackey functors that we consider. In terms of binary (and $n$-ary) operations, there are a few, including:

            @UL {
                ("the ", TO2(boxProduct,"box product"), " of two Mackey functors, via  ", TT "M**N"),
                ("the ", TO2(directSum,"direct sum")," of two Mackey functors, as ", TT"M++N"),
                ("the ", TO2(InternalHom,"internal hom"), " of two Mackey functors")
            }@
    SeeAlso
        "background on Mackey functors"
        CpMackeyFunctor
///

doc ///
    Key
        makeCpMackeyFunctor
        (makeCpMackeyFunctor,ZZ,Matrix,Matrix,Matrix)
    Headline
        make a Cp Mackey functor
    Usage
        makeCpMackeyFunctor(p,R,T,C)
    Inputs
        p : ZZ
            a prime number
        R : Matrix
            a restriction homomorphim
        T : Matrix
            a transfer homomorphim
        C : Matrix
            a conjugation homomorphim
    Outputs
        : CpMackeyFunctor
    Description
        Text
            Builds an instance of the @TO("CpMackeyFunctor")@ type from the data of a prime number and three matrices, encoding restriction, transfer, and conjugation, respectively. The method will verify whether the input data yields a valid Mackey functor.
    SeeAlso
        "background on Mackey functors"
        CpMackeyFunctor
///

doc ///
    Key
        (prune,CpMackeyFunctor)
    Headline
        prune a CpMackeyFunctor
    SeeAlso
        "background on Mackey functors"
        CpMackeyFunctor
///

doc ///
    Key
        isCohomological
        (isCohomological,CpMackeyFunctor)
    Headline
        checks if a CpMackeyFunctor is cohomological
    Usage
        isCohomological M
    Inputs
        M : CpMackeyFunctor
    Outputs
        : Boolean
            whether $M$ is {\em cohomological} or not
    Description
        Text
            A Mackey functor $M$ defined over $C_p$ is said to be {\em cohomological} if the composite of restriction followed by transfer is identical to multiplication by $p$ as a map from the @TO2(Fixed,"fixed module")@ of $M$ to itself.
    SeeAlso
        "background on Mackey functors"
        CpMackeyFunctor
///

doc ///
    Key
        drawVerticalCpMackeyFunctor
        (drawVerticalCpMackeyFunctor,CpMackeyFunctor)
    Headline
        reorients the net of a Mackey functor vertically
    Usage
        drawVerticalCpMackeyFunctor M
    Inputs
        M : CpMackeyFunctor
    Outputs
        : Net
    Description
        Text
            Reorients the @TO2(Net,"net")@ of a Mackey functor vertically
        Example
            M = makeRandomCpMackeyFunctor(2)
            drawVerticalCpMackeyFunctor M
    SeeAlso
        "background on Mackey functors"
        CpMackeyFunctor
///