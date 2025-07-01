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