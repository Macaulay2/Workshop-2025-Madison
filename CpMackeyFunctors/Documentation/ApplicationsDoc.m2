doc ///
    Key
        "explicit applications of the CpMackeyFunctors package"
    Headline
        info about how this package can be applied to real-life examples
    Description
        Text
            This package can be used to perform Tor and Ext computations for $C_p$-Mackey functors.  It also implements Tor and Ext in the abelian subcategory of cohomological $C_p$-Mackey functors.
            Computations of Tor and Ext groups for (cohomological) Mackey functors are important for work in equivariant homotopy theory.  As a demonstration, we present some computations which appear in the work of Mingcong Zeng [Zeng18].
            As setup, let $B_1$ be the $C_p$-Mackey functor with fixed module $\mathbb{Z}/p$ and underlying module $0$.  Let $Z$ denote the @TO2(makeFixedPointMackeyFunctor,"fixed point Mackey functor")@ for $\ZZ$ as a trivial $C_p$-module (called $\square$ on the @TO2("list of common Mackey functors","list here")@).
            
            {\bf Proposition:} For any prime number $p$ we have  
            \[
                \mathrm{Ext}^i_{Z}(B_1,Z) = \begin{cases} B_1 & i=3 \\ 0. & \textrm{else.} \end{cases}
            \]
            Note that by a theorem of Arnold (see [Arn] and also [BSW]) the category of cohomological Mackey functors has global dimension $3$, so it suffices to compute the first four Ext Mackey functors. For concreteness, we will do the computation at the prime 11. 
        Example
            B = cokernel(matrix({{11}}));
            B1 = makeZeroOnUnderlyingMackeyFunctor(11,B)
            Z = makeFixedPointMackeyFunctor(11,id_(ZZ^1))
            for i to 3 do (print (prune ExtCoh(i,B1,Z)))
        Text
            We can also recover some Tor computations which appear in Zeng's work.  Let $DZ$ denote the Mackey functor which is $\mathbb{Z}$ at fixed and underlying, transfer and conjugation are the identity, and restriction is multiplication by $p$. 

            {\bf Proposition:} For any prime number $p$ we have 
            \[
                \mathrm{Tor}^{Z}_i(B_1,B_1) = \begin{cases} B_1 & i=0,3 \\ 0 & \textrm{else.} \end{cases}
            \]
            and 
            \[
                \mathrm{Tor}^{Z}_i(DZ,DZ) = \begin{cases} DZ & i=0 \\ B_1 & i=1 \\ 0 & \textrm{else.} \end{cases}
            \]
            Again, we will perform these computations at the prime $p=11$. 
        Example
            B = cokernel(matrix({{11}}));
            B1 = makeZeroOnUnderlyingMackeyFunctor(11,B)
            DZ = makeOrbitMackeyFunctor(11,id_(ZZ^1))
            for i to 3 do (print (prune TorCoh(i,B1,B1)))
            for i to 3 do (print (prune TorCoh(i,DZ,DZ)))
    SeeAlso
        "background on Mackey functors"
        "the abelian category of Mackey functors"
        "constructing examples of Mackey functors"
        "list of common Mackey functors"
///
