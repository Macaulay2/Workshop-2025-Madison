document {
    Key=> {
        boxProduct,
        (boxProduct,CpMackeyFunctor,CpMackeyFunctor),
        (symbol ⊠,CpMackeyFunctor,CpMackeyFunctor),
    },
    Headline => "box product of Cp-Mackey functors",
    Usage => "M ⊠ N\nboxProduct(M,N)",
    Inputs =>{
        CpMackeyFunctor => "M" => {""},
        CpMackeyFunctor => "N" => {""},
    },
    Outputs =>{
        CpMackeyFunctor => {"the box product of ",TEX"$M$"," and ", TEX"$N$", "."},
    },
    PARA{
        "Given two ", TO2(CpMackeyFunctor,"Cp-Mackey functors"), " ", TEX"$M$"," and ", TEX"$N$", ", we can form their box product, defined by taking the levelwise box product of the underlying and fixed modules."
    },
}