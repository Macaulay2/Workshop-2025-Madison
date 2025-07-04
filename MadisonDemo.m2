-- Madison Workshop Demo

-- Cayley Omega Process

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

n = 2;
r = 2;
R = QQ[x_(1,1)..x_(n,n)];
f = random (r*n, R);
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

n = 2;
r = 2;
p = 1;
R = QQ[x_(1,1)..x_(n,n)];
f = random (r*n, R);
reynoldsSLnMap(f, p)