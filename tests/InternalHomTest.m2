needsPackage "CpMackeyFunctors"

p := 5;
A := makeBurnsideMackeyFunctor p;
B := makeUnderlyingFreeMackeyFunctor p;
RO := makeRealRepresentationMackeyFunctor p;

assert(rank getUnderlyingModule InternalHom(A,B) == rank Hom(getUnderlyingModule A, getUnderlyingModule B))
assert(rank getFixedModule InternalHom(A,B) == rank Hom(A, B))

assert(rank getUnderlyingModule InternalHom(B,RO) == rank getUnderlyingModule (B ** RO))
assert(rank getFixedModule InternalHom(B,RO) == rank getFixedModule (B ** RO))