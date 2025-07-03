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
