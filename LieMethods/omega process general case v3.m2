
--******** omega constant c_p,n*************


cayleyConstant = method()
cayleyConstant (ZZ, ZZ) := (p,n) -> (
    RRRR = QQ[x_(1,1)..x_(n,n)];
    M = genericMatrix(RRRR, x_(1,1), n, n);
    detM = det M;
    cpn = detM^p;
    
    for i from 1 to p do (
        cpn = diff(detM, cpn)
    );
    sub(cpn, ZZ)
)

cayleyConstant(2,2)

--******** reynolds GLn (4.5.27)*************

reynoldsGLnMap = method()
reynoldsGLnMap RingElement := f -> (
    R = ring f;
    n = floor(sqrt(numgens R));
    
    d = (degree f)#0;
    if d % n != 0 then (
        error "degree of f must be divisible by the number of variables"
    );
    M = genericMatrix(R, n, n);
    p = floor(d / n);
    detM = det M;
    omegaf = f * (detM)^p;
    for i from 1 to p do (
        omegaf = diff(detM, omegaf)
    );
    cpn = cayleyConstant(p, n);
    omegaf / sub(cpn, R)
)

n = 2
r = 2
R = QQ[x_(1,1)..x_(n,n)];
f = random (r*n, R)
reynoldsGLnMap(f)


--******** reynolds SLn (4.5.28)*************

reynoldsSLnMap = method()
reynoldsSLnMap (RingElement, ZZ) := (f, p) -> (
    R = ring f;
    n = floor(sqrt(numgens R));
    d = (degree f)#0;
    if d % n != 0 then (
        error "degree of f must be divisible by the number of variables"
    );
    r = floor(d / n);
    M = genericMatrix(R, n, n);
    detM = det M;
    omegaf = f * (detM)^p;
    for i from 1 to r do (
        omegaf = diff(detM, omegaf)
    );
    crn = cayleyConstant(r,n);
    phi = map(ring omegaf,R, gens(ring omegaf));
    phi(sub(detM^(r-p), R)) * omegaf / phi(sub(crn, R))
)

n = 2
r = 2
p = 1
R = QQ[x_(1,1)..x_(n,n)];
f = random (r*n, R)
reynoldsSLnMap(f, p)




--******** SL2 invariants of Sym2 (4.5.31)*************

actionBySLn = method()

actionBySLn (ZZ, ZZ) := (n, d) -> (
    N = binomial(n+d-1, d); -- number of coefficients in the polynomial

    -- Building the set of general coefficients
    U = QQ[X_1..X_n, MonomialOrder => Lex];
    allMonomials = (entries basis(d, U))#0;
    aVars = apply(exponents(sum(allMonomials)), p -> a_p);

    -- Building the polynomial ring with the coefficients
    S = QQ[aVars, X_1..X_n];

    -- Mapping the set of monomials into the polynomial ring
    fromUtoS = map(S, U, take(gens S, -n));
    allMonomials = allMonomials / fromUtoS;

    -- Building the universal form of the polynomial of degree d in n variables
    allMonomialsMatrix = matrix {allMonomials};
    aVarsMatrix = matrix {(gens S)_{0..N-1}};
    universaldForm = (allMonomialsMatrix * transpose(aVarsMatrix))_0_0;

    -- ****** linear change of variables & extract coefficients ******

    zVars = z_(1,1)..z_(n,n);

    R = QQ[aVars, zVars, X_1..X_n];
    R' = QQ[aVars, zVars][X_1..X_n];

    fromStoR = map(R, S, join(take(gens R, N),take(gens R, -n)));
    universaldForm = fromStoR(universaldForm);
    aVarsMatrix = matrix {(gens S)_{0..N-1} / fromStoR};
    fromRtoR' = map(R', R, gens coefficientRing R' | gens R');

    zVarsMatrix = genericMatrix(R, (z_(1,1))_R, n, n);
    XVarsMatrix = genericMatrix(R, (X_1)_R, n, 1);

    SLnLinearChangeofVars = map(R,R, aVarsMatrix | matrix mutableMatrix(R, 1, n^2) | transpose(zVarsMatrix * XVarsMatrix));
    transformeddForm = SLnLinearChangeofVars(universaldForm)
)

g = actionBySLn(2, 2)




a0image = coefficient(X_1^2,fromRtoR'(transformeddForm)) -- see page 197
a1image = coefficient(X_1*X_2,fromRtoR'(transformeddForm))
a2image = coefficient(X_2^2,fromRtoR'(transformeddForm))

-- ** computing R(a1^2)

a1sqimage = a1image^2

omegaf = a1sqimage
for i from 1 to 2 do (
    omegaf = diff(det zVarsMatrix, sub(omegaf, R))
);
cayleyConstant(2,2)
omegaf / sub(cayleyConstant(2,2), R)
-- need the map from R to QQ[...]