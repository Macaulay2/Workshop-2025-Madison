doc ///
    Key
        CpMackeyFunctor
    Headline
        The type of Cp Mackey functors

    Description
        Text
            This documents the @TO2(Type,"type")@ of $C_p$-Mackey functors. A $C_p$-Mackey functor is encoded as a @TO2(HashTable, "hash table")@ with keys through which we can pull various parts of the data.
        Text
            The data available in the hash table is:
            
            @UL {
                (BOLD "underlying module: ", 
                "the module ", 
                TEX"$M(C_p/e)$",
                ". It can be accessed via the ",
                TO "getUnderlyingModule", " method."),
                (BOLD "fixed module: ",
                "the module ",
                TEX"$M(C_p/C_p)$",
                ". It can be accessed via the ",
                TO "getUnderlyingModule", " method."),
                (BOLD "restriction: ",
                "the restriction homomorphism ",
                TEX"$M(C_p/C_p) \\to M(C_p/e)$",
                ". It can be accessed via the ",
                TO "getRestriction", " method."),
                (BOLD "transfer: ",
                "the transfer homomorphism ",
                TEX"$M(C_p/e) \\to M(C_p/C_p)$",
                ". It can be accessed via the ",
                TO "getTransfer", " method."),
                (BOLD "conjugation: ",
                "the conjugation homomorphism ",
                TEX"$M(C_p/e) \\to M(C_p/e)$",
                ". It can be accessed via the ",
                TO "getConjugation", " method."),
            }@
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