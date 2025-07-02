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
                (TT ".PrimeOrder", " yields the ", TO2(ZZ,"order"), " of the cyclic group over which we're working. Accesible via the  ", TO "getPrimeOrder", " method"),
                (TT ".Underlying", 
                " encodes the ", TO "module", " ",
                TEX"$M(C_p/e)$",
                ". It can be accessed via the ",
                TO "getUnderlyingModule", " method."),
                (TT ".Fixed",
                " encodes the ",TO "module", " ",
                TEX"$M(C_p/C_p)$",
                ". It can be accessed via the ",
                TO "getFixedModule", " method."),
                (TT ".Res",
                " encodes the restriction ", TO2("Macaulay2Doc :: module homomorphisms","homomorphism")," ",
                TEX"$M(C_p/C_p) \\to M(C_p/e)$",
                ". It can be accessed via the ",
                TO "getRestriction", " method."),
                (TT ".Tr",
                " encodes the transfer ", TO2("Macaulay2Doc :: module homomorphisms","homomorphism")," ",
                TEX"$M(C_p/e) \\to M(C_p/C_p)$",
                ". It can be accessed via the ",
                TO "getTransfer", " method."),
                (TT ".Conj",
                " encodes the conjugation ", TO2("Macaulay2Doc :: module homomorphisms","homomorphism")," ",
                TEX"$M(C_p/e) \\to M(C_p/e)$",
                ". It can be accessed via the ",
                TO "getConjugation", " method."),
            }@

            {\bf Constructing new Mackey functors:}

            {\bf Operations with Mackey functor types:} One of the first operations to know about is @TO2((prune,CpMackeyFunctor),"pruning")@ a $C_p$-Mackey functor. There is more information at the prune command, but this basically takes the data of a Mackey functor and attempts to simplify it and make it more reader-friendly. Pruning is the only unary operation on $C_p$-Mackey functors that we consider. In terms of binary (and $n$-ary) operations, there are a few, including:

            @UL {
                ("the ", TO2(boxProduct,"box product"), " of two Mackey functors, via  ", TT "M**N"),
                ("the ", TO2(directSum,"direct sum")," of two Mackey functors, as ", TT"M++N"),
                ("the ", TO2(InternalHom,"internal hom"), " of two Mackey functors")
            }@
    ///

doc ///
    Key
        getPrimeOrder
    Headline
        getPrimeOrder
///

doc ///
    Key
        getUnderlyingModule
    Headline
        getUnderlyingModule
///

doc ///
    Key
        getFixedModule
    Headline
        getFixedModule
///

doc ///
    Key
        getRestriction
    Headline
        getRestriction
///

doc ///
    Key
        getTransfer
    Headline
        getTransfer
///

doc ///
    Key
        getConjugation
    Headline
        getConjugation
///

doc ///
    Key
        makeCpMackeyFunctor
        (makeCpMackeyFunctor,ZZ,Matrix,Matrix,Matrix)
    Headline
        make a Cp Mackey functor
///

doc ///
    Key
        (prune,CpMackeyFunctor)
    Headline
        prune a CpMackeyFunctor
///