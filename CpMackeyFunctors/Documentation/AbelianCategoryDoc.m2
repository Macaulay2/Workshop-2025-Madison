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