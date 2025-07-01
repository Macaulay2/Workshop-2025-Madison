restart
path = append(path, "./Durham_2024/")
needsPackage "EquivariantEhrhart"

-- Ex 1

M = matrix{{1,0,0},{0,1,0},{0,0,1}}
P = convexHull M
vertices P

gList = {matrix{{1,0,0},{0,1,0},{0,0,1}}}
