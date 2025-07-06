doc ///
    Key
        getUnstableGlobalA1Degree
        (getUnstableGlobalA1Degree, frac R)   
    Headline
        computes the unstable global A1 degree of a pointed rational function $f/g:\mathbb{P}^{1}_{k}\to\mathbb{P}^{1}_{k}$
    Usage
        getUnstableGlobalA1Degree(q)
    Inputs 
        q :frac R
             a polynomial or a rational function @TT("f/g")@ where @TT("f")@ and @TT("g")@ are polynomials in one variable over a field @TT("k")@, with @TT("g")@ not identically zero   
    Outputs
        : Sequence
             a B\'{e}zoutian matrix @TT("N")@ and its scalar @TT("t") determinant as an UnstableGrothendieckWittClass
    Description
        Text
             As termed by Igieobo, et. al [I24], the unstable $\mathbb{A}^1$-degree is valued in the unstable Grothendieck-Witt ring $GW^u(k)$, where there is a group isomorphism $[\mathbb{P}^{1}_{k},\mathbb{P}^{1}_{k}]\longrightarrow\mathrm{GW}(k)\times_{k^{\times}/(k^{\times})^{2}}k^{\times}.$  As in the stable case  [BMP23], the global unstable $\mathbb{A}^{1}$-degree at rational points is given as a type of B\'{e}zoutian.  There is a group isomorphism  [C09]  $$\deg^{u}:[\mathbb{P}^{1}_{k},\mathbb{P}^{1}_{k}]\longrightarrow\mathrm{GW}^{u}(k)$$
    by $\deg^{u}(f/g)=(\mathrm{B\acute{e}z}^{\mathrm{mon}}(f/g),\det\mathrm{B\acute{e}z}^{\mathrm{mon}}(f/g))$.  Taking in a rational function @TT("f/g")@, this method returns the Grothendieck-Witt class of the B\'{e}zoutian matrix and its determinant as an unstable Grothendieck-Witt class.  
        Example
            q = (x^2 + x - 2)/(3*x + 5);
            getGlobalUnstableA1Degree(q)
    References
        [I24] Igieobo et. al, "Motivic Configurations on the Line", arXiv:2411.15347 [math.AT], 2024.
        [BMP23] T. Brazelton, S. McKean, S. Pauli, Bezoutians and the A1-Degree, Algebra & Number Theory, 2023.
        [C09] C. Casanave, Algebraic Homotopy Classes of Rational Functions", arXiv: 0912.2227 [math.AT], 2009.
    SeeAlso 
        getLocalUnstableA1Degree
///

doc///
    Key
        getLocalUnstableA1Degree
        (getLocalUnstableA1Degree, RingElement, RingElement)
    Headline
        computes a local unstable A1-Brouwer degree of a pointed rational function at a root
    Usage
        getLocalUnstableA1Degree(q, r)
    Inputs
        q : a rational function @TT("f/g")@ where @TT("f")@ and @TT("g")@ are polynomials in one variable over a field @TT("k")@, with @TT("g")@ not identically zero
        r : a root of the polynomial @TT("f")@ in the ring of fractions of the polynomial ring over @TT("k")
    Outputs
        : UnstableGrothendieckWittClass`
\\\