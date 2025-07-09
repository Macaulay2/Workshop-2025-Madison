doc ///
    Key
        addGWu
        (addGWu, UnstableGrothendieckWittClass,UnstableGrothendieckWittClass)     
    Headline
        the direct sum for  two unstable Grothendieck-Witt Classes
    Usage
        addGWu(beta, gamma)
    Inputs
        beta : UnstableGrothendieckWittClass
            the isomorphism class of a non-degenerate symmetric bilinear form represented by a matrix @TT("M")@ together with a scalar @TT("s")@
        gamma : UnstableGrothendieckWittClass
            the isomorphism class of a non-degenerate symmetric bilinear form represented by a matrix @TT("N")@ together with a scalar @TT("t")@

    Outputs
        : UnstableGrothendieckWittClass
            the isomorphism class of the direct sum of the bilinear forms represented by the matrices @TT("M")@  and @TT("N")@, and the scalar is the product @TT("st")@
    Description
        Text
            This computes the direct sum of the Grothendieck-Witt classes @TT("beta")@ and @TT("gamma")@.
        Example
            M = matrix(QQ, {{2,1},{1,2}})
            N = matrix(QQ, {{1,2},{2,6}})
            beta = makeGWuClass M
            gamma = makeGWuClass N
            addGWu(beta, gamma)
///