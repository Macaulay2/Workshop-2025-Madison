document{
    Key => {map, (map, CpMackeyFunctor, CpMackeyFunctor, Matrix, Matrix)},
    Headline => "constructs a map between two Mackey functors",
    Usage => "map(N,M,U,F)",
    Inputs => {
	    CpMackeyFunctor => "N" => {""},
        CpMackeyFunctor => "M" => {""},
        Matrix => "U" => {"represents the map on the underlying level"},
        Matrix => "F" => {"represents the map on the fixed-point level"},
	},
    Outputs => {
	    MackeyFunctorHomomorphism => {"the morphism of Mackey functors which is F at the fixed point level and U at the underlying level."},
	},
    PARA{
        "A morphisms of ", TEX///$C_p$///, "-Mackey functors consists of a group homomorphism on the fixed and underlying levels which commutes with the transfer, restriction, and conjugation morphisms.  This method will throw an error if the morphism is not well defined."
    },
    EXAMPLE lines///
    map(makeComplexRepresentationMackeyFunctor 3, makeBurnsideMackeyFunctor 3, matrix {{1}}, matrix {{1,1},{0,1},{0,1}})
    ///
}

document{
    Key => {source, (source, MackeyFunctorHomomorphism)},
    Headline => "returns the source of a Mackey functor homomorphism",
    Usage => "map(f)",
    Inputs => {
        MackeyFunctorHomomorphism => "f" => {""},
	},
    Outputs => {
	    CpMackeyFunctor => {"the source of the homomorphism f."},
	},
    PARA{
        "A homomorphism between Mackey functors has a source and a target. This method returns the source."
    },
    EXAMPLE lines///
    source(complexLinearizationMap(5))
    ///
}

document{
    Key => {target, (target, MackeyFunctorHomomorphism)},
    Headline => "returns the source of a Mackey functor homomorphism",
    Usage => "map(f)",
    Inputs => {
        MackeyFunctorHomomorphism => "f" => {""},
	},
    Outputs => {
	    CpMackeyFunctor => {"the target of the homomorphism f."},
	},
    PARA{
        "A homomorphism between Mackey functors has a source and a target. This method returns the target."
    },
    EXAMPLE lines///
    target(complexLinearizationMap(5))
    ///
}

document{
    Key => {complexLinearizationMap, (complexLinearizationMap, ZZ)},
    Headline => "returns the complex linearization map for the prime p",
    Usage => "complexLinearizationMap(p)",
    Inputs => {
        ZZ=> "p" => {"a prime"},
	},
    Outputs => {
	    MackeyFunctorHomomorphism => {"the linearization map from the Burnside Mackey functor to the complex representation functor for the prime p."},
	},
    PARA{
        "Every ", TEX///$G$///,"-set X determines a complex representation with basis X, called the permutation representation of X.  This induces a map from the Burnside Mackey functor to the complex representation Mackey functor called the linearization map."
    },
    EXAMPLE lines///
    complexLinearizationMap(5)
    ///
}

document{
    Key => {realLinearizationMap, (realLinearizationMap, ZZ)},
    Headline => "returns the real linearization map for the prime p",
    Usage => "realLinearizationMap(p)",
    Inputs => {
        ZZ=> "p" => {"a prime"},
	},
    Outputs => {
	    MackeyFunctorHomomorphism => {"the linearization map from the Burnside Mackey functor to the real representation functor for the prime p."},
	},
    PARA{
        "Every ", TEX///$G$///,"-set X determines a real representation with basis X, called the permutation representation of X.  This induces a map from the Burnside Mackey functor to the real representation Mackey functor called the linearization map."
    },
    EXAMPLE lines///
    realLinearizationMap(5)
    ///
}