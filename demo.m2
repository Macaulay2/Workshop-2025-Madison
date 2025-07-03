needsPackage "ExteriorResolutions"

(S, E) = koszulPair(4, ZZ/101)
isCommutative S
isSkewCommutative E

f = map(S^{-1}^4, S^4, matrix { { S_0, S_1, S_2, S_3 }, {S_3, S_0, S_1, S_2}, {S_2, S_3, S_0, S_1}, {S_1, S_2, S_3, S_0} })

g = map(S^4, S^{1}^4, matrix { { S_0, S_1, S_2, S_3 }, {S_1, S_2, S_3, S_0}, {S_2, S_3, S_0, S_1}, {S_3, S_0, S_1, S_2} })

f * g

koszulRR(f, Concentration => (-5,5)) * koszulRR(g, Concentration => (-5,5)) == koszulRR(f * g, Concentration => (-5,5))

koszulRR(complex(coker vars S), Concentration => (-5,5))

C = koszulRR(complex S, Concentration => (-5,5))

P = priddyComplex(vars E, S, LengthLimit => 5)

betti P

betti C

betti C == betti P

koszulLL(complex(coker vars E), Concentration => (-5,5))

C = prune koszulLL(complex E, Concentration => (-5,5))

koszulComplex(vars S)

betti C == betti (koszulComplex(vars S)[5]**S^{5})

----

(S,E) = koszulPair(2, ZZ/101)

C = complex {map(S^{1}, S^1, S_0)}

D = koszulRR(C, Concentration=>(-5,5))

D.dd_1
D.dd_0

----

C = koszulLL(E^1, Concentration => (-3, 0))
D = koszulRR(S^1, Concentration => (-5, 0))

complex E == naiveTruncation(prune HH koszulRR(C, Concentration => (-5,0)), -4, 0)
complex comodule truncate(6, S) == prune HH koszulLL(D, Concentration => (-2, 4))