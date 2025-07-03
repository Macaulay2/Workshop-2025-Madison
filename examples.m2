restart
path = append(path, "./Durham_2024/")
needsPackage "EquivariantEhrhart"

-- Ex 1

M = matrix{{1,0,0},{0,1,0},{0,0,1}}
P = convexHull M
vertices P

gList = {matrix{{0,1,0},{1,0,0},{0,0,1}}}

equivariantEhrhartSeries(P,gList)

-- Ex 2
-- Have to double check why this one has a negative coefficient in highest degree
M = matrix{{1,1,0,0},{0,1,1,0},{0,0,1,1},{1,0,0,1}}
P = convexHull M
vertices P

gList = {matrix{{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,0,0,0}}}

equivariantEhrhartSeries(P,gList)


--
P = convexHull matrix "-1, 2/3"
vertices P
ehrhartQP P
displayQP ehrhartQP P

ehrhartSeries P

-----------------

end

-- other examples 


Example [Reeve tetrahedra; 2.4 https://arxiv.org/pdf/2307.10852]

P_q := Convex hull of {(0,0,0), (1,0,0), (0,1,0), (1,q,q+1)}

if q = 12 then the Ehrhart polynomial has a negative coefficient

check that:
i_(P_q) = (q+1)/6 d^3 + d^2 + (11-q)/6 d + 1
h* = q t^2 + 1


Example Standard reflexive simplex
P_d := convex hull e_1 .. e_d, -(e_1 + .. + e_d)
h* = t^d + .. + t + 1


Example HKN18: In RR^5:
P = convexHull 0, e_1 .. e_4, (5,5,5,5,8)
h* = 1 + t + 2t^2 + t^3 + 2t^4 + t^5
hStarVector = {1,1,2,1,2,1}


Example:[Section 2, https://arxiv.org/pdf/2110.10204]
P = [-1, 2/3] in RR^1



