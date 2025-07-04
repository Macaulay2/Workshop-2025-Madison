------------------------------------------------------
-- Demo for homological algebra over exterior algebras
-- 2025 Macaulay2 workshop in Madison, WI
------------------------------------------------------

-- Team members:
--  Andrew Karsten
--  Boyana Martinova
--  Caitlin M. Davis
--  Jiaheng Li
--  Jianuo Zhou
--  Mahrud Sayrafi
--  Penelope Beall
--  Sreehari Suresh Babu
--
-- Team leads:
--  Michael K. Brown
--  Gregory G. Smith

restart
needsPackage "ExteriorResolutions"

-----------------------------------------------
-- Injective resolutions over exterior algebras
-----------------------------------------------

-- Let's start with defining an exterior algebra
E = ZZ/11[e_0..e_3, SkewCommutative => true]
e_1*e_0
e_0^2
k = coker vars E -- the residue field is an E-module

-- We compute the injective resolution of k over E
elapsedTime P = injectiveResolution(k, LengthLimit => 5)
augP = coaugmentationMap P
isQuasiIsomorphism augP
-- In fact, this is the Priddy complex on {e_3,e_2,e_1,e_0}
-- whose terms are P'_i = E ** (S_i)^*
elapsedTime P' = priddyComplex(matrix{{e_3,e_2,e_1,e_0}}, LengthLimit => 5)
P == P'


-- A non-cyclic Example
f = random(E^1, E^{-1,-1})
M = k ++ image f
P = injectiveResolution(M, LengthLimit => 5)

augP = coaugmentationMap P
isQuasiIsomorphism augP


-- Resolving a Complex
I = ideal random(E^1, E^{-1,-1})
C = Hom(freeResolution(I, LengthLimit => 5), comodule I)
P = injectiveResolution(C, LengthLimit => 5)

augP = coaugmentationMap P;
isQuasiIsomorphism augP

------------------------------------------
-- Koszul duality functors between E and S
------------------------------------------

-- Let's define a Koszul pair S, E
(S, E) = koszulPair(n = 2, ZZ/101)
isCommutative S -- the polynomial ring
isCommutative E -- the exterior algebra
isSkewCommutative E

-- RR(k) = E
koszulRR coker vars S

-- RR(S) = k
C = koszulRR(S^1, Concentration => (-5,0))
prune HH_0 C == coker vars E
-- Moreover, it coincides with the Priddy complex on {e_0,e_1,e_2}
C == priddyComplex(vars E, LengthLimit => 5)

-- LL(k) = S
koszulLL coker vars E

-- LL(E(n+1)) = k
C = koszulLL E^{n+1}
-- In fact, this is the Priddy complex on {x_0,x_1,x_2} with a twist and shift
C == (priddyComplex(-vars S) ** S^{-3})[-3]
-- Which is also isomorphic to the Koszul complex on {x_0,x_1,x_2}
f = randomComplexMap(C, koszulComplex vars S, Cycle => true)
isQuasiIsomorphism f

------------------------------------------
-- Composition of LL and RR
------------------------------------------

-- LL(RR(S)) = S, but because of bounds
-- we get a finite length quotient of it
koszulRR(S^1, Concentration => (-5,0))
koszulLL oo
prune HH oo
oo == complex comodule truncate(5+1, S)

-- LL(RR(finite length S-module)) is identity
M = comodule truncate(n+2, S)
koszulRR oo
koszulLL oo
complex M == prune HH oo

-- RR(LL(E)) = E, but because of bounds
-- the ends of the complex may not be exact
koszulLL E^1
C = koszulRR(oo, Concentration => (-5,0))
prune HH C
complex E == prune HH canonicalTruncation(C, (-4,0))

-- RR(LL(perfect complex of E-modules)) is identity
N = coker matrix {{e_0}}
F = freeResolution(N, LengthLimit => n)
F.dd
koszulLL F
koszulRR(oo, Concentration => (-n-1,n))
prune canonicalTruncation(oo, (-n, n-1))
complex N == prune HH oo

------------------------------------------
-- Functoriality of LL and RR
------------------------------------------

-- LL is functorial (and RR too)
f = random(E^{2:1}, E^2)
g = random(E^2, E^{2:-1})
koszulLL f * koszulLL g == koszulLL(f * g)
koszulLL f

C = complex {map(S^{1}, S^1, S_0)}
D = koszulRR(C, Concentration=>(-5,1))
D.dd_1
D.dd_0
