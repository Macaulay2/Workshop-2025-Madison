-- Macaulay2 demo for the elementary abelian p groups.

load "invariantSeedGeneration.m2";


-- seed generation
-- define the ring
R = QQ[x_0..x_2];
-- define our dimensions for the problem:
W = matrix{{0,1,1},{1,0,1}}; -- here it's easiest to work when all the dimensions are the same
ZList = {2,2,2};
theseeds = genseeds(R,W,ZList);

-- expanding the seeds
expandseeds(theseeds, ZList);


-- then combine the two functions into a composite function:

-- Initial 3 variable example
R = QQ[x_0..x_2];
W = matrix{{0,1,1},{1,0,1}};
ZList = {2,2,2};
invargenset(R, W, ZList)

-- 5 variable example
R = QQ[v,w,x,y,z];
ZList = {4,4,4,4,4};
W = matrix{{1,1,1,1,1},{1,1,1,1,0}};
invargenset(R,W,ZList);