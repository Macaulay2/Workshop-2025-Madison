doc ///
    Key
        getUnstableGlobalA1Degree
        (getUnstableGlobalA1Degree, RingElement)   
    Headline
        computes the unstable global A1 degree of a rational function
    Usage
        getUnstableGlobalA1Degree(q)
    Inputs 
        q : a rational function @TT("f/g")@ where @TT("f")@ and @TT("g")@ are polynomials in one variable over a field @TT("k")@, with @TT("g")@ not identically zero   
    Outputs
    Description
     As termed by Igieobo, et. al, the unstable $\mathbb{A}^1$-degree is valued in the unstable Grothendieck-Witt ring $GW^u(k)$, where there is a group isomorphism $[\mathbb{P}^{1}_{k},\mathbb{P}^{1}_{k}]\longrightarrow\mathrm{GW}(k)\times_{k^{\times}/(k^{\times})^{2}}k^{\times}.$  As in the stable case, the global unstable $\mathbb{A}^{1}$-degree at rational points is given as a type of B\'{e}zoutian.  Work by Casanave shows that there is a group is a isomorphism  $$\deg^{u}:[\mathbb{P}^{1}_{k},\mathbb{P}^{1}_{k}]\longrightarrow\mathrm{GW}^{u}(k)$$
    by $\deg^{u}(f/g)=(\mathrm{B\acute{e}z}^{\mathrm{mon}}(f/g),\det\mathrm{B\acute{e}z}^{\mathrm{mon}}(f/g))$.  Taking in a rational function @TT("f/g")@, the unstable global $\mathbb{A}^{1}$-degree is computed as the Grothendieck-Witt class of the B\'{e}zoutian matrix of @TT("f/g")@ together with the determinant of that matrix.    
        Text
        Example
///