(S, E) = koszulPair(4, ZZ/19937)
M = matrix { { S_0, S_1, S_2, S_3 }, {S_3, S_0, S_1, S_2}, {S_2, S_3, S_0, S_1}, {S_1, S_2, S_3, S_0} }
f = map(S^{-1}^4, S^4, M)
N = matrix { { S_0, S_1, S_2, S_3 }, {S_1, S_2, S_3, S_0}, {S_2, S_3, S_0, S_1}, {S_3, S_0, S_1, S_2} }
g = map(S^4, S^{1}^4, N)
koszulRR(f, Concentration=>(-5,5)) * koszulRR(g, Concentration=>(-5,5)) == koszulRR(f * g, Concentration=>(-5,5))