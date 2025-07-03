-- most naive approach, produces really big free modules
makeFreeModuleSurjection = method()
makeFreeModuleSurjection(CpMackeyFunctor) := MackeyFunctorHomomorphism => (M) -> (
    return makeUniversalMap(M, gens(M.Underlying), gens(M.Fixed))
)

-- get resolution up to F_(n-1) <-- F_n
getResolution = method()
getResolution(CpMackeyFunctor,ZZ) := List => (M,n) -> (
    if n < 0 then return {};
    k := 0;
    if M.cache#?"ProjRes" then (
        k = length M.cache#"ProjRes";
    ) else (
        M.cache#"ProjRes" = {makeFreeModuleSurjection M};
        k = 1
    );

    while k <= n do (
        -- add differential d_k to cache
        d := (M.cache#"ProjRes")#(k-1);
        M.cache#"ProjRes" |= {inducedMap(source d, kernel d) * (makeFreeModuleSurjection (kernel d))};
        k += 1
    );

    M.cache#"ProjRes"_(for i to n list i)
)
