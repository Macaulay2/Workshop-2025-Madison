-- most naive approach, produces really big free modules
makeFreeModuleSurjection = method()
makeFreeModuleSurjection(CpMackeyFunctor) := MackeyFunctorHomomorphism => (M) -> (
    return makeUniversalMap(M, gens(M.Fixed), gens(M.Underlying))
)
