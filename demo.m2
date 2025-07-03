needsPackage "ExteriorResolutions"

(S, E) = koszulPair(4, ZZ/101)

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

(S, E) = koszulPair(2, ZZ/101)

C = koszulLL(complex E, Concentration => (-5,5))

koszulComplex(vars S)

betti C == betti (koszulComplex(vars S)[5]**S^{5})

D = koszulRR(C, Concentration => (-5,5))

prune HH D == complex E