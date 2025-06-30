needsPackage "CpMackeyFunctors"

x := Hom(makeBurnsideMackeyFunctor 3, makeBurnsideMackeyFunctor 3)
assert(isFreeModule prune x and rank x == 2)

y := Hom(makeUnderlyingFreeMackeyFunctor 5, makeUnderlyingFreeMackeyFunctor 5)
assert(isFreeModule prune y and rank y == 5)