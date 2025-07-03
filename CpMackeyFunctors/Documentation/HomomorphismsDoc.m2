doc ///
    Key
        MackeyFunctorHomomorphism
    Headline
        the class of Mackey functor homomorphisms
    Subnodes
        FixedMap
        UnderlyingMap
    Description
        Text
            Common ways to construct a Mackey functor homomorphism:
        Text
            • @ TO2((map, CpMackeyFunctor, CpMackeyFunctor, Matrix, Matrix), "map") @  
        Text
            • @ TT "id_M" @
///
-- ^^ There must be a better way!

doc ///
    Key
        FixedMap
    Headline
        the fixed-point level of a Mackey functor homomorphism
    Usage
        F.FixedMap
///

doc ///
    Key
        UnderlyingMap
    Headline
        the underlying level of a Mackey functor homomorphism
    Usage
        F.UnderlyingMap
///

doc ///
    Key
        (map,CpMackeyFunctor, CpMackeyFunctor, Matrix, Matrix)
        [map,Degree]
        [map,DegreeLift]
        [map,DegreeMap]
    Headline
        constructs a map between two Mackey functors
    Usage
        map(N,M,U,F)
    Inputs
        N : CpMackeyFunctor
        M : CpMackeyFunctor
        U : Matrix
        F : Matrix
    Outputs
        : MackeyFunctorHomomorphism
            the morphism of Mackey functors which is F at the fixed point level and U at the underlying level.
    Description
        Text
            A morphisms of $C_p$-Mackey functors consists of a group homomorphism on the fixed and underlying levels which commutes with the transfer, restriction, and conjugation morphisms.  This method will throw an error if the morphism is not well defined.
        Example
            map(makeComplexRepresentationMackeyFunctor 3, makeBurnsideMackeyFunctor 3, matrix {{1}}, matrix {{1,1},{0,1},{0,1}})
///

doc ///
    Key
        (source, MackeyFunctorHomomorphism)
    Headline
        returns the source of a Mackey functor homomorphism
    Usage
        source(F)
    Inputs
        F : MackeyFunctorHomomorphism
    Outputs
        : CpMackeyFunctor
            the source of the homomorphism F.
    Description
        Text
            A homomorphism between Mackey functors has a source and a target. This method returns the source.
        Example
            source(complexLinearizationMap(5))
///

doc ///
    Key
        (target, MackeyFunctorHomomorphism)
    Headline
        returns the target of a Mackey functor homomorphism
    Usage
        target(F)
    Inputs
        F : MackeyFunctorHomomorphism
    Outputs
        : CpMackeyFunctor
            the target of the homomorphism F.
    Description
        Text
            A homomorphism between Mackey functors has a source and a target. This method returns the target.
        Example
            target(complexLinearizationMap(5))
///

doc ///
    Key
        complexLinearizationMap
        (complexLinearizationMap, ZZ)
    Headline
        returns the complex linearization map for the prime p
    Usage
        complexLinearizationMap(p)
    Inputs
        p : ZZ
    Outputs
        : MackeyFunctorHomomorphism
            the linearization map from the Burnside Mackey functor to the complex representation functor for the prime p.
    Description
        Text
            Every $G$-set $X$ determines a complex representation with basis $X$, called the permutation representation of $X$.  This induces a map from the @TO2(makeBurnsideMackeyFunctor,"Burnside Mackey functor")@ to the @TO2(makeComplexRepresentationMackeyFunctor,"complex representation Mackey functor")@ called the linearization map.
        Example
            complexLinearizationMap(5)
    SeeAlso
        realLinearizationMap
///

doc ///
    Key
        realLinearizationMap
        (realLinearizationMap, ZZ)
    Headline
        returns the real linearization map for the prime p
    Usage
        realLinearizationMap(p)
    Inputs
        p : ZZ
    Outputs
        : MackeyFunctorHomomorphism
            the linearization map from the Burnside Mackey functor to the real representation functor for the prime p.
    Description
        Text
            Every $G$-set $X$ determines a real representation with basis $X$, called the permutation representation of $X$.  This induces a map from the Burnside Mackey functor to the real representation Mackey functor called the linearization map.
        Example
            realLinearizationMap(5)
    SeeAlso
        complexLinearizationMap
///

doc ///
    Key
        (symbol |, MackeyFunctorHomomorphism,MackeyFunctorHomomorphism)
    Headline
        Horizontal concatenation of Mackey functor homomorphisms
    Usage
        f | g
    Outputs
        : MackeyFunctorHomomorphism
            the horizontal concatenation of the homomorphisms $f$ and $g$.
    Inputs
        f : MackeyFunctorHomomorphism
        g : MackeyFunctorHomomorphism
    Description
        Text
            Given two Mackey functor homomorphisms $f : A \to C$ and $g : B \to C$ with the same codomain, $f | g$ is the homomorphism $(f,g) : A \oplus B \to C$.
        Example
            A = makeBurnsideMackeyFunctor 2;
            id_A | id_A
///

doc ///
    Key
        (symbol ||, MackeyFunctorHomomorphism,MackeyFunctorHomomorphism)
    Headline
        Vertical concatenation of Mackey functor homomorphisms
    Usage
        f || g
    Outputs
        : MackeyFunctorHomomorphism
            the vertical concatenation of the homomorphisms $f$ and $g$.
    Inputs
        f : MackeyFunctorHomomorphism
        g : MackeyFunctorHomomorphism
    Description
        Text
            Given two Mackey functor homomorphisms $f : A \to B$ and $g : A \to C$ with the same domain, $f || g$ is the homomorphism $(f,g) : A \to B \oplus C$.
        Example
            A = makeBurnsideMackeyFunctor 2;
            id_A || id_A
///

doc ///
    Key
        isTrivialMackeyFunctor
        (isTrivialMackeyFunctor,CpMackeyFunctor)
    Headline
        checks if a Mackey functor is trivial or not
    Usage
        isTrivialMackeyFunctor M
    Inputs
        M : CpMackeyFunctor
    Outputs
        : Boolean
            whether $M$ is trivial or not
    Description
        Text
            Checks whether $M$ is {\em trivial} (meaning isomorphic to the @TO2(makeZeroMackeyFunctor,"zero Mackey functor")@).

    SeeAlso
        "background on Mackey functors"
        makeZeroMackeyFunctor
///