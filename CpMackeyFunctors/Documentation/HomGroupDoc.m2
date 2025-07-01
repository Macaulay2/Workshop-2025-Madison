document{
    Key => {(Hom, CpMackeyFunctor,CpMackeyFunctor),[Hom,DegreeLimit],[Hom,MinimalGenerators],[Hom,Strategy]},
    Headline => "computes the Hom group between two Cp Mackey functors",
    Usage => "Hom(M,N)",
    Inputs => {
	CpMackeyFunctor => "M" => {"a ", TEX///$C_p$///," Mackey functor M"},
    CpMackeyFunctor => "N" => {"a ", TEX///$C_p$///," Mackey functor N"},
	},
    Outputs => {
	Module => {"The hom group as a ", TEX///$\mathbb{Z}$///,"-module"},
	},
    PARA{
        "Given two Mackey functors ", TEX"$M$", " and ", TEX"$N$"," the set of Mackey functors homomorphisms from ",TEX"$N$"," to ",TEX"$M$"," forms an abelian group under pointwise addition.",
    },
    EXAMPLE lines///
    M = makeBurnsideMackeyFunctor(2);
    N = makeUnderlyingFreeMackeyFunctor(2);
    Hom(M,N)
    ///
}