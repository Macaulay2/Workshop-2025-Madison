doc ///
    Key
        getUnstableGlobalA1Degree
        (getUnstableGlobalA1Degree, RingElement)   
    Headline
        computes the unstable global A1 degree of a pointed rational function f/g
    Usage
        getUnstableGlobalA1Degree(q)
    Inputs 
        q : a rational function @TT("f/g")@ where @TT("f")@ and @TT("g")@ are polynomials in one variable over a field @TT("k")@, with @TT("g")@ not identically zero   
    Outputs
        : UnstableGrothendieckWittClass
            the isomorphism class of a non-degenerate symmetric bilinear form represented by a B\'{e}zoutian matrix @TT("N")@ together with a scalar @TT("t") as its determinant
    Description
        Text
             As termed by Igieobo, et. al, the unstable $\mathbb{A}^1$-degree is valued in the unstable Grothendieck-Witt ring $GW^u(k)$, where there is a group isomorphism $[\mathbb{P}^{1}_{k},\mathbb{P}^{1}_{k}]\longrightarrow\mathrm{GW}(k)\times_{k^{\times}/(k^{\times})^{2}}k^{\times}.$  As in the stable case, the global unstable $\mathbb{A}^{1}$-degree at rational points is given as a type of B\'{e}zoutian.  Work by Casanave shows that there is a group is a isomorphism  $$\deg^{u}:[\mathbb{P}^{1}_{k},\mathbb{P}^{1}_{k}]\longrightarrow\mathrm{GW}^{u}(k)$$
    by $\deg^{u}(f/g)=(\mathrm{B\acute{e}z}^{\mathrm{mon}}(f/g),\det\mathrm{B\acute{e}z}^{\mathrm{mon}}(f/g))$.  Taking in a rational function @TT("f/g")@, this method returns the Grothendieck-Witt class of the B\'{e}zoutian matrix and its determinant as an unstable Grothendieck-Witt class.  
        Example
            q = (x^2 + x - 2)/(3*x + 5);
            getGlobalUnstableA1Degree(q)
           
///
