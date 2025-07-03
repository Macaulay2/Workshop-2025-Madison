computeHomology = method()
computeHomology(MackeyFunctorHomomorphism,MackeyFunctorHomomorphism) := CpMackeyFunctor => (f,g) -> (
    i := inducedMap(source(f),ker(f));
    phi := inducedMap(ker(f),source(g),g);

    coker phi
) 

