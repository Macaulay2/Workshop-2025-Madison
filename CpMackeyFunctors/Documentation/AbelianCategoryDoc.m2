doc ///
    Key
        "the abelian category of Mackey functors"
    Headline
        background on the category of Mackey functors and how to code within it
    Description
        Text
            {\bf The category of Mackey functors:} Fixing a finite group $G$ (always cyclic of prime order in this package), we obtain a {\em category} $\text{Mack}_G$, whose objects are @TO2(CpMackeyFunctor,"Mackey functors")@ and whose morphisms are @TO2(MackeyFunctorHomomorphism,"Mackey functor homomorphisms")@. Given any Mackey functor $M$ we can take its {\em identity homomorphism} as follows:
        
        Example
            M = makeRandomCpMackeyFunctor(3,GenBound=>2)
            id_M
        
        Text
            {\bf The abelian category structure:} Given any two Mackey functors $M$ and $N$, we can form the set $\text{Hom}_{\text{Mack}_G}(M,N)$ of homomorphisms from $M$ to $N$. Since we can @TO2((symbol ++, MackeyFunctorHomomorphism,MackeyFunctorHomomorphism),"add")@ any two homomorphisms, we have that $\text{Hom}_{\text{Mack}_G}(M,N)$ is really an @TO2((Hom,CpMackeyFunctor,CpMackeyFunctor) ,"abelian group of homomorphisms")@. Composition is moreover $\ZZ$-linear, hence we can say that the category $\text{Mack}_G$ is {\em pre-additive}.

        Example
            Hom(makeRandomCpMackeyFunctor(2,GenBound=>2),makeRandomCpMackeyFunctor(2,GenBound=>2))
            prune oo
        
        Text
            Given any two Mackey functors, we can take their @TO2(directSum,"direct sum")@ and obtain another Mackey functor. This is a {\em finite biproduct} on the category $\text{Mack}_G$, in other words it is both a product and a coproduct. There is a zero object for this biproduct, namely the @TO2(makeZeroMackeyFunctor,"zero Mackey functor")@. This gives $\text{Mack}_G$ what is often called a {\em semiadditive} structure. 
            
            Moreover, any Mackey functor homomorphism admits a @TO2((kernel,MackeyFunctorHomomorphism),"kernel")@ and a @TO2((cokernel,MackeyFunctorHomomorphism),"cokernel")@. All this data we've collected so far is often called a {\em pre-abelian} structure. Finally, one can prove that images and coimages agree, hence $\text{Mack}_G$ is an {\em abelian category}.

            {\bf The monoidal structure:} The category $\text{Mack}_G$ furthermore has a {\em symmetric monoidal structure}, given by the @TO2(boxProduct,"box product")@ of two Mackey functors. Taking the box product admits a right adjoint, in other words there is an @TO2(InternalHom,"internal hom")@. More explicitly, there is a natural isomorphism of abelian groups, natural in any three Mackey functors $M$, $N$, and $P$:
            \[\text{Hom}_{\text{Mack}_G}(M \square N, P) \cong \text{Hom}_{\text{Mack}_G}(M, [N,P]).\]

            {\bf Free resolutions:} todo

            {\bf Ext and Tor:} todo

    SeeAlso
        "background on Mackey functors"
        "constructing examples of Mackey functors"
        "explicit applications of the CpMackeyFunctors package"

///

doc ///
    Key
        (symbol ++,CpMackeyFunctor,CpMackeyFunctor)
    Headline
        direct sum of Cp-Mackey functors
    Usage
        M ++ N
    Inputs
        M : CpMackeyFunctor
        N : CpMackeyFunctor
    Outputs
        : CpMackeyFunctor
            the direct sum of $M$ and $N$
    Description
        Text
            Given two @TO2((CpMackeyFunctor),"Cp-Mackey Functors")@ $M$ and $N$, we can form their {\em direct sum}, defined by taking the levelwise direct sum of the underlying and fixed modules.
///

doc ///
    Key
        (kernel, MackeyFunctorHomomorphism)
        [kernel,DegreeLimit]
        [kernel,Strategy]
        [kernel,SubringLimit]
    Headline
        kernel of a Mackey functor homomorphism
    Usage 
        kernel f
    Inputs
        f : MackeyFunctorHomomorphism
    Outputs
        : CpMackeyFunctor
            the kernel of the homomorphism f
///

doc ///
    Key
        (cokernel, MackeyFunctorHomomorphism)
    Headline
        cokernel of a Mackey functor homomorphism
    Usage
        cokernel f
    Inputs
        f : MackeyFunctorHomomorphism
    Outputs
        : CpMackeyFunctor
            the cokernel of the homomorphism f
///