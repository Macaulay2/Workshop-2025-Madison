document {
    Key=> {
        (symbol ++,CpMackeyFunctor,CpMackeyFunctor),
    },
    Headline => "direct sum of Cp-Mackey functors",
    Usage => "M ++ N",
    Inputs =>{
        CpMackeyFunctor => "M" => {""},
        CpMackeyFunctor => "N" => {""},
    },
    Outputs =>{
        CpMackeyFunctor => {"the direct sum of ",TEX"$M$"," and ", TEX"$N$", "."},
    },
    PARA{
        "Given two ", TO2(CpMackeyFunctor,"Cp-Mackey functors"), " ", TEX"$M$"," and ", TEX"$N$", ", we can form their direct sum, defined by taking the levelwise direct sum of the underlying and fixed modules."
    },
}

document {
    Key => {
        (kernel, MackeyFunctorHomomorphism),
    },
    Headline => "kernel of a Mackey functor homomorphism",
    Usage => "kernel f",
    Inputs => {
        MackeyFunctorHomomorphism => "f" => {""},
    },
    Outputs => {
        CpMackeyFunctor => {"the kernel of the homomorphism f."},
    },
}

document {
    Key => {
        (cokernel, MackeyFunctorHomomorphism),
    },
    Headline => "cokernel of a Mackey functor homomorphism",
    Usage => "cokernel f",
    Inputs => {
        MackeyFunctorHomomorphism => "f" => {""},
    },
    Outputs => {
        CpMackeyFunctor => {"the cokernel of the homomorphism f."},
    },
}