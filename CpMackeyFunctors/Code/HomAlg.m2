-- given two composable morphisms f,g such that f*g == 0, return ker(f)/im(g)
computeHomology = method()
computeHomology(MackeyFunctorHomomorphism,MackeyFunctorHomomorphism) := CpMackeyFunctor => (f,g) -> (
    comp := f*g;
    if (comp != map(target comp, source comp, 0)) then (
        error "inputted maps do not form a complex"
    );
    i := inducedMap(source(f),ker(f));
    phi := inducedMap(ker(f),source(g),g);
    coker phi
) 

-- will resolve second argument
Tor(ZZ,CpMackeyFunctor,CpMackeyFunctor) := CpMackeyFunctor => opts -> (i,M,N) -> (
    d := getResolution(N,i+1);
    if i == 0 then (
        coker (M**(d_1))
    ) else (
        computeHomology (M**d_i), (M**d_(i+1))
    )
)
