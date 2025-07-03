-- tests
load "seedGenerator.m2";
R = QQ[x_1..x_n]
p = 3
W = matrix{{1,0,0,1,1,0},{1,0,1,0,1,0}}

elementaryAbelianPGroupSeedGenerator(W,3,R)