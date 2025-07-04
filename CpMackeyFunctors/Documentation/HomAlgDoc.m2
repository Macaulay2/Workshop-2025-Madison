doc ///
Node
    Key
        Tor
        (Tor,ZZ,CpMackeyFunctor,CpMackeyFunctor)
    Headline
        computes Tor of two Cp-Mackey fuctors
    Usage
        Tor_i(M,N)
    Inputs
        i : ZZ
        M : CpMackeyFunctor
        N : CpMackeyFunctor
    Outputs
        : CpMackeyFunctor
            the $i$th Tor Cp-Mackey Functor of $M$ and $N$
    Description
        Text
            Given two @TO2(CpMackeyFunctor,"Cp-Mackey functors")@ $M$ and $N$, we can compute their $i$th Tor.
        Example
            RU = makeComplexRepresentationMackeyFunctor(3);
            Tor_3(RU,RU)

    SeeAlso
        "the abelian category of Mackey functors"
        (TorCoh,ZZ,CpMackeyFunctor,CpMackeyFunctor)
        (Ext,ZZ,CpMackeyFunctor,CpMackeyFunctor)
        (ExtCoh,ZZ,CpMackeyFunctor,CpMackeyFunctor)
Node
    Key
        TorCoh
        (TorCoh,ZZ,CpMackeyFunctor,CpMackeyFunctor)
    Headline
        computes Tor of two cohomological Cp-Mackey fuctors
    Usage
        TorCoh(i,M,N)
    Inputs
        i : ZZ
        M : CpMackeyFunctor
        N : CpMackeyFunctor
    Outputs
        : CpMackeyFunctor
            the $i$th cohomological Tor Cp-Mackey Functor of $M$ and $N$
    Description
        Text
            Given two cohomological @TO2(CpMackeyFunctor,"Cp-Mackey functors")@ $M$ and $N$, we can compute their $i$th Tor in the category of cohomological Mackey functors.
        Example
            N = cokernel(matrix({{3}}));
            M = makeZeroOnUnderlyingMackeyFunctor (3,N);
            TorCoh(3,M,M)
    SeeAlso
        "the abelian category of Mackey functors"
        (Tor,ZZ,CpMackeyFunctor,CpMackeyFunctor)
        (Ext,ZZ,CpMackeyFunctor,CpMackeyFunctor)
        (ExtCoh,ZZ,CpMackeyFunctor,CpMackeyFunctor)
///


doc ///
Node
    Key
        Ext
        (Ext,ZZ,CpMackeyFunctor,CpMackeyFunctor)
    Headline
        computes Ext of two Cp-Mackey fuctors
    Usage
        Ext^i(M,N)
    Inputs
        i : ZZ
        M : CpMackeyFunctor
        N : CpMackeyFunctor
    Outputs
        : CpMackeyFunctor
            the $i$th Ext Cp-Mackey Functor of $M$ and $N$
    Description
        Text
            Given two @TO2(CpMackeyFunctor,"Cp-Mackey functors")@ $M$ and $N$, we can compute their $i$th Ext.
        Example
            RU = makeComplexRepresentationMackeyFunctor(3);
            prune Ext^4(RU,RU)

    SeeAlso
        "the abelian category of Mackey functors"
        (ExtCoh,ZZ,CpMackeyFunctor,CpMackeyFunctor)
        (Tor,ZZ,CpMackeyFunctor,CpMackeyFunctor)
        (TorCoh,ZZ,CpMackeyFunctor,CpMackeyFunctor)
Node
    Key
        ExtCoh
        (ExtCoh,ZZ,CpMackeyFunctor,CpMackeyFunctor)
    Headline
        computes Ext of two cohomological Cp-Mackey fuctors
    Usage
        Ext(i,M,N)
    Inputs
        i : ZZ
        M : CpMackeyFunctor
        N : CpMackeyFunctor
    Outputs
        : CpMackeyFunctor
            the $i$th cohomological Ext Cp-Mackey Functor of $M$ and $N$
    Description
        Text
            Given two cohomological @TO2(CpMackeyFunctor,"Cp-Mackey functors")@ $M$ and $N$, we can compute their $i$th Ext Mackey functor in cohomological Mackey functors.
        Example
            N = cokernel(matrix({{3}}));
            M = makeZeroOnUnderlyingMackeyFunctor (3,N);
            prune ExtCoh(3,M,M)
    SeeAlso
        "the abelian category of Mackey functors"
        (Ext,ZZ,CpMackeyFunctor,CpMackeyFunctor)
        (Tor,ZZ,CpMackeyFunctor,CpMackeyFunctor)
        (TorCoh,ZZ,CpMackeyFunctor,CpMackeyFunctor)
///