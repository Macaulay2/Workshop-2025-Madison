doc ///
    Key
        (InternalHom, CpMackeyFunctor, CpMackeyFunctor)
    Headline
        returns the internal hom Mackey functor between two Mackey functors.
    Usage
        InternalHom(N,M)
    Inputs
        M : CpMackeyFunctor
        N : CpMackeyFunctor
    Outputs
        : CpMackeyFunctor
            the internal hom Mackey functor from M to N.
    Description
        Text
            Mackey functors form a closed symmetric monoidal category.  This method returns the internal hom Mackey functor.
        Example
            InternalHom(makeRealRepresentationMackeyFunctor 5, makeComplexRepresentationMackeyFunctor 5)
///

doc ///          
    Key
        (InternalHom, CpMackeyFunctor, MackeyFunctorHomomorphism)
    Headline
        returns the induced map on an internal hom.
    Usage
        InternalHom(F,M)
        InternalHom(M,F)
    Inputs
        M : CpMackeyFunctor
        F : MackeyFunctorHomomorphism
    Outputs
        : MackeyFunctorHomomorphism
            the induced map on internal hom Mackey functors
    Description
        Text
            The internal hom of Mackey functors is functorial in each variable.  This method returns the induced maps.
        Example
            InternalHom(makeRealRepresentationMackeyFunctor(5), complexLinearizationMap(5))
///

--this last doesn't work

-- doc ///          
--     Key
--         (InternalHom, MackeyFunctorHomomorphism, CpMackeyFunctor)
--     Headline
--         returns the induced map on an internal hom.
--     Usage
--         InternalHom(F,M)
--         InternalHom(M,F)
--     Inputs
--         M : CpMackeyFunctor
--         F : MackeyFunctorHomomorphism
-- ///