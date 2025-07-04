doc ///
Node
    Key
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
Node
    Key
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

///


doc ///
    Key
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
///