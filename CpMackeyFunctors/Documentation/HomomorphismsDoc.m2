document{
    Key => {MackeyFunctorHomomorphism},
    Headline => "the class of Mackey functor homomorphisms",
    PARA{
        "Common ways to construct a Mackey functor homomorphism:",
    },
    UL {
        TO2((map, CpMackeyFunctor, CpMackeyFunctor, Matrix, Matrix), "map"),
        TT "id_M",
    },
    Subnodes => {
        TO FixedMap,
        TO UnderlyingMap,
    }
}


document{
    Key => {FixedMap},
    Headline => "the fixed-point level of a Mackey functor homomorphism",
    Usage => "f.FixedMap",
}

document{
    Key => {UnderlyingMap},
    Headline => "the underlying level of a Mackey functor homomorphism",
    Usage => "f.UnderlyingMap",
}

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
            Every $G$-set $X$ determines a complex representation with basis $X$, called the permutation representation of $X$.  This induces a map from the Burnside Mackey functor to the complex representation Mackey functor called the linearization map.
        Example
            complexLinearizationMap(5)
///

doc ///
    Key
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
///

doc ///
    Key
        (makeUniversalMapFixed, CpMackeyFunctor, Vector)
    Headline
        constructs a homomorphism out of the Burnside Mackey functor
    Usage
        makeUniversalMapFixed(M,x)
    Inputs
        M : CpMackeyFunctor
        x : Vector
    Outputs
        : MackeyFunctorHomomorphism
            the homomorphism $\underline{A}\to M$ induced by the vector $x$.
    Description
        Text
            The Burnside Mackey functor $\underline{A}$ is the free Mackey functor on a fixed generator, meaning there is a natural isomorphism $\text{Hom}_{\text{Mack}_{C_p}}(\underline{A},M)\cong M(C_p/C_p)$ for any $C_p$-Mackey functor $M$. This method implements this correspondence: for an element $x\in M(C_p/C_p)$ it returns the induced homomorphism of Mackey functors $\underline{A}\to M$. The homomorphism is determined by sending $1\mapsto x$ where $1\in\underline{A}(C_p/C_p)$ represents the singleton $C_p$-set.
        Example
            makeUniversalMapFixed(makeRealRepresentationMackeyFunctor(5), vector (matrix {{1},{2},{3}}))
///

-- document{
--     Key => {makeUniversalMapFixed, (makeUniversalMapFixed, CpMackeyFunctor, Vector)},
--     Headline => "constructs a homomorphism out of the Burnside Mackey functor",
--     Usage => "makeUniversalMapFixed(M,x)",
--     Inputs => {
--         CpMackeyFunctor => "M" => {""},
--         Vector => "x" => {"a vector in the fixed module ", TEX///$x\in M(C_p/C_p)$///}
--     },
--     Outputs => {
--         MackeyFunctorHomomorphism => {"the homomorphism ", TEX///$\underline{A}\to M$///, " induced by the vector ", TEX///$x$///}
--     },
--     PARA{
--         "The Burnside Mackey functor ", TEX///$\underline{A}$///, " is the free Mackey functor on a fixed generator, meaning there is a natural isomorphism ", TEX///$\text{Hom}_{\text{Mack}_{C_p}}(\underline A,M)\cong M(C_p/C_p)$///, " for any ", TEX///$C_p$///, "-Mackey functor ", TEX///$M$///, ". This method implements this correspondence: for an element ", TEX///$x\in M(C_p/C_p)$///, ", it returns the induced homomorphism of Mackey functors ", TEX///$\underline A\to M$///, ". The homomorphism is determined by sending ", TEX///$1\mapsto x$///, ", where ", TEX///$1\in\underline A(C_p/C_p)$///, " represents the singleton ", TEX///$C_p$///, "-set."
--     },
--     EXAMPLE lines///
--     RO = makeRealRepresentationMackeyFunctor 5;
--     makeUniversalMapFixed(RO, vector (matrix {{1},{2},{3}}))
--     ///,
--     SeeAlso => {"makeBurnsideMackeyFunctor", "makeFixedFreeMackeyFunctor"}
-- }

-- document{
--     Key => {makeUniversalMapUnderlying, (makeUniversalMapUnderlying, CpMackeyFunctor, Vector)},
--     Headline => "constructs a homomorphism out of the free Mackey functor on an underling generator",
--     Usage => "makeUniversalMapUnderlying(M,x)",
--     Inputs => {
--         CpMackeyFunctor => "M" => {""},
--         Vector => "x" => {"a vector in the underlying module ", TEX///$x\in M(C_p/e)$///}
--     },
--     Outputs => {
--         MackeyFunctorHomomorphism => {"the homomorphism ", TEX///$\underline{B}\to M$///, " induced by the vector ", TEX///$x$///}
--     },
--     PARA{
--         "The free Mackey functor ", TEX///$\underline{B}$///, " on an underlying generator has the property there is a natural isomorphism ", TEX///$\text{Hom}_{\text{Mack}_{C_p}}(\underline B,M)\cong M(C_p/e)$///, " for any ", TEX///$C_p$///, "-Mackey functor ", TEX///$M$///, ". This method implements this correspondence: for an element ", TEX///$x\in M(C_p/e)$///, ", it returns the induced homomorphism of Mackey functors ", TEX///$\underline B\to M$///, ". The homomorphism is determined by sending ", TEX///$1\mapsto x$///, ", where ", TEX///$1\in\underline B(C_p/e)$///, " represents the identity element ", TEX///$1\in C_p/e = \{1,\gamma,\ldots,\gamma^{p-1}\}$///, "."
--     },
--     EXAMPLE lines///
--     B = makeUnderlyingFreeMackeyFunctor 3;
--     makeUniversalMapUnderlying(B, vector (matrix {{1},{2},{3}}))
--     ///,
--     SeeAlso => {"makeUnderlyingFreeMackeyFunctor"}
-- }
