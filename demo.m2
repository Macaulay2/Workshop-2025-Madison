restart

-- The ehrhart method doesn't give a warning when it's applied to a non-integral rational polytope, and sometimes its give the wrong polynomial.

P = convexHull transpose matrix {{0,0}, {1,0}, {0,1/2}}
#latticePoints(3*P) -- => 6
substitute(ehrhart(P), x => 3) -- => 7

needsPackage("RationalPolytopes", FileName => "Durham_2024/RationalPolytopes.m2")
apply(1..3, t -> #latticePoints(t*P)) -- => (2, 4, 6)
f = ehrhartQP P
apply(1..3, t -> (ehrhartQP P) t) -- => (2, 4, 6)

-- Example of a rational polytope with a denominator larger than the period.

restart

needsPackage("RationalPolytopes", FileName => "Durham_2024/RationalPolytopes.m2")

P = convexHull transpose matrix {{-1,0}, {0,1/2}, {0,-1/2}, {1,0}}
denominator P -- => 2
period ehrhartQP P -- => 1
