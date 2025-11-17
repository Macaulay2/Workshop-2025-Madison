uninstallPackage RationalPolytopes

restart
installPackage "RationalPolytopes"


P = convexHull matrix "-1, 2/3"
EP = ehrhartQP P
ES = ehrhartSeries P
value ES
numerator ES
denominator ES

hStarVector P

-- period collapse phenomenon
P = convexHull transpose matrix "1,0; -1,0; 0,1/2; 0,-1/2"
EP = ehrhartQP P
ES = ehrhartSeries P
value ES
numerator ES
denominator ES

hStarVector P


-------------------------

P = convexHull transpose matrix {
    {1,1,1, 1,1,1, 1,1,1},
    {0,1,1, 1,1,1, 1,1,2},
    {0,1,1, 0,1,1, 1,2,2},
    {0,0,1, 1,1,2, 1,1,2},
    {0,0,1, 0,1,2, 1,2,2},
    {0,1/2,1, 1/2,3/2,3/2, 1,3/2,3/2},
    {1/2,1/2,1, 1/2,1/2,3/2, 1,3/2,2}
    }

ehrhartQP P

Q = convexHull transpose matrix {
    {0,0,0,0},
    {1,0,0,0},
    {0,1,0,0},
    {0,0,1,0},
    {1,1,1,2}
    }

ehrhartQP Q

ESP = ehrhartSeries P
ESQ = ehrhartSeries Q

numerator ESP

factor value ESP
factor value ESQ

faces P
vertices P


-- middle vertices
M = (vertices P)_{0,1,2,3,4} - (vertices P)_{0,0,0,0,0}

rank M

--
--           5
--          / \
--         /   \
--        /     \
--       0-1-2-3-4  <-- 3 dimensional subspace 'M'
--        \     /
--         \   /
--          \ /
--           6
--

reducedRowEchelonForm transpose M

T = (vertices P)_{5,6} - (vertices P)_{0,0}

T % M



P' = convexHull transpose matrix {
    {-1,0,0,0},
    {1,1,0,0},
    {0,-1/2,0,0}
    }

vertices P'
ehrhartQP P'
factor value ehrhartSeries P'

E = convexHull transpose matrix {
    {0,0,0,0},
    {0,0,0,1},
    {0,0,1,0}
    }
ehrhartSeries E

dim P
dim P'
dim E

P'' = convexHull(P', E)
vertices P''
dim P''

factor value ehrhartSeries P''
factor value ehrhartSeries P


P3 = convexHull transpose matrix {
    {-1,0},
    {1,1},
    {0,-1/2}
    }
factor value ehrhartSeries P3

-- idea: want to prove P and P'' are mutation equivalent
-- then show that Q and P'' are not mutation equivalent
vertices P''
vertices P
vertices Q
