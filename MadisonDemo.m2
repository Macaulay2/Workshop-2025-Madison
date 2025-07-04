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

-- Casimir and Reynolds Operators for the Lie algebra sln
-- Written by Jeong-Hoon Ju
--******** Lie Algebra sln (Example 4.5.14) ********

Eijm = (i0,j0,m) -> ( matrix apply(m, i -> apply(m, j -> if i==i0 and j==j0 then 1/1 else 0/1)) );
Hin = (i,n) -> ( Eijm(i,i,n) - Eijm(i+1,i+1,n));
HDualin= (i,n) -> (((n-i)/n)*(sum for j from 0 to i-1 list Eijm(j,j,n))-(i/n)*(sum for j from 1 to n-i list Eijm(n-j,n-j,n)))

slnBasis = (n) -> (
    B:={};
    Hbasis := apply(n-1, i -> Hin(i,n));
    Xbasis := flatten apply(n, i -> delete(null,apply(n, j -> if i<j then Eijm(i,j,n))));   
    Ybasis := flatten apply(n, i -> delete(null,apply(n, j -> if i>j then Eijm(i,j,n)))); 
    flatten {Hbasis, Xbasis, Ybasis}
)

slnHbasis = (n) -> (
    B:={};
    Hbasis := apply(n-1, i -> Hin(i,n))
);

slnXbasis = (n) -> (
    B:={};
    Xbasis := flatten apply(n, i -> delete(null,apply(n, j -> if i<j then Eijm(i,j,n))))
);

slnYbasis = (n) -> (
    B:={};
    Ybasis := flatten apply(n, i -> delete(null,apply(n, j -> if i>j then Eijm(i,j,n))))
);

slnDualBasis = (n) -> (
    B:={};
    HDualbasis := apply(n-1, i -> HDualin(i+1,n));
    XDualbasis := flatten apply(n, i -> delete(null,apply(n, j -> if i<j then Eijm(j,i,n))));   
    YDualbasis := flatten apply(n, i -> delete(null,apply(n, j -> if i>j then Eijm(j,i,n))));
    L:= flatten {HDualbasis, XDualbasis, YDualbasis};
    apply(L, i -> (1/(2*n))*i)
)

slnHDualbasis = (n) -> (
    B:={};
    HDualbasis := apply(n-1, i -> HDualin(i+1,n))
);

slnXDualbasis = (n) -> (
    B:={};
    XDualbasis := flatten apply(n, i -> delete(null,apply(n, j -> if i<j then Eijm(j,i,n))))
);

slnYDualbasis = (n) -> (
    B:={};
    YDualbasis := flatten apply(n, i -> delete(null,apply(n, j -> if i>j then Eijm(j,i,n))))
);


--******** Casimir Operator for the represention Sym^d(QQ^n) of sln ********

slnElementActionMatrixEntry = (n,d,k,i,j) -> (
    B = slnBasis(n);
    U=QQ[X_1..X_n];
    allMonomials = (entries basis(d, U))#0;
    f=0;
    for p from 1 to n do(
	for q from 1 to n do(
	   f=f+((B)_k)_(p-1,q-1)*X_p*(diff(X_q,allMonomials_j));
	);
    );
    (M,C) = coefficients(f,Monomials => allMonomials);
    sub(C_(i,0),QQ)
);

slnDualActionMatrix = (n,d,k) -> (
    N = binomial(n+d-1, d);
    M = matrix table(N,N,(i,j) -> slnElementActionMatrixEntry(n,d,k,i,j));
    M = transpose M;
    M = -M
);

slnCoordiAction = (n,d,k,S,F) -> (
    M = slnDualActionMatrix(n,d,k);
    N = binomial(n+d-1, d);
    aVars = gens S;
    actionF=0;
    for i from 0 to N-1 do (
	for j from 0 to N-1 do (
	    actionF = actionF + M_(i,j)*aVars_i*diff(aVars_j,sub(F,S));
	);
    );
    actionF
)

slnDualElementActionMatrixEntry = (n,d,k,i,j) -> (
    B = slnDualBasis(n);
    U=QQ[X_1..X_n];
    allMonomials = (entries basis(d, U))#0;
    f=0;
    for p from 1 to n do(
	for q from 1 to n do(
	   f=f+((B)_k)_(p-1,q-1)*X_p*(diff(X_q,allMonomials_j));
	);
    );
    (M,C) = coefficients(f,Monomials => allMonomials);
    sub(C_(i,0),QQ)
);

slnDualDualActionMatrix = (n,d,k) -> (
    N = binomial(n+d-1, d);
    M = matrix table(N,N,(i,j) -> slnDualElementActionMatrixEntry(n,d,k,i,j));
    M = transpose M;
    M = -M
);

slnDualCoordiAction = (n,d,k,S,F) -> (
    M = slnDualDualActionMatrix(n,d,k);
    N = binomial(n+d-1, d);
    actionF=0;
    aVars = gens S;
    for i from 0 to N-1 do (
	for j from 0 to N-1 do (
	    actionF = actionF + M_(i,j)*aVars_i*diff(aVars_j,sub(F,S));
	);
    );
    actionF
)

slnCasimir = (n,d,S,F) -> (
    N = binomial(n+d-1, d);
    actionF=0;
    for k from 0 to N-1 do (
	actionF = actionF + slnCoordiAction(n,d,k,S,slnDualCoordiAction(n,d,k,S,F));
    );
    actionF
)

--******** Example 4.5.21 ********

n = 2;
d = 2;

N = binomial(n+d-1, d);
U=QQ[X_1..X_n];
allMonomials = (entries basis(d, U))#0;
aVars = apply(exponents(sum(allMonomials)), p -> a_p);
S = QQ[aVars]; -- S is the coordinate ring for the space of homogeneous polynomials of degree d with n variables

f1= a_{1,1}^2;
slnCasimir(2,2,S,f1) -- n=2, d=2, F is a polynomial in S
slnCasimir(2,2,S,slnCasimir(2,2,S,f1))

f2= a_{1,1}^4;
slnCasimir(2,2,S,f2) -- n=2, d=2, F is a polynomial in S


--******** Reynolds Operator for the represention Sym^d(QQ^n) of sln (Algorithm 4.5.20) ********

--ver1: with Gordie Novak
slnReynolds = (n,d,S,f) -> (
    fList = {f};
    gList = fList;
    
    Rf = 0;

    b = new MutableList;

    fin = 0;
    for l when (l < #gList and gList_l != 0) do (
        
        b = append(b, {});

        for i from 0 to l-1 do(
            b#l = b#l | {0};
        );
        b#l = b#l | {1};

        a = leadCoefficient(gList_(l));
        for i when i<l do(
            if (leadMonomial(gList_(l)) == leadMonomial(gList_(i))) then (
                gList = replace(l, gList_(l)-a*(gList_(i)), gList);
                for j from 0 to i do ( 
                    b#l = replace(j, ((b#l_j) - (a*(b#i_j))), b#l); 
                );
                a=leadCoefficient(gList_(l));
            ); 
        );
        if a != 0 then(
            gList = replace(l, (gList_(l))/a, gList);

            for j from 0 to l do( 
                b#l = replace(j, (b#l_j) / a, b#l); 
            );
            fList = append(fList, slnCasimir(n,d,S,fList_(l)));
            gList = append(gList, fList_(l+1));
        )
        else (
            if (gList_l == 0) then (
                fin = l;
                break;
            );
            l = l-1;
        );
    );

    if gList_(fin) == 0 then ( 
        if (b#fin_0) != 0 then (
            Rf = 0;
        ) else (
            Rf = ( sum for j from 1 to fin list 
                (b#fin_j)*fList_(j-1) ) 
                / (b#fin_1); 
        );
    );

    print "\n--------";
    print "Printing Rf;";
    return Rf;
)

--ver2: with Kisun Lee
slnReynolds = (n,d,S,f) -> (
    fList = {f};
    gList = fList;
    l = 0;
    Rf = 0;
    while last gList != 0 do (
	if gList_(l) !=0 then(
	    b_(l,l) = 1;
	    apply(l, i -> b_(l,i) = 0);
	    a = leadCoefficient(gList_(l));
	    apply(l, i -> (
		    if  leadMonomial(gList_(l)) == leadMonomial(gList_(i))
		    then (
			gList = replace(l, gList_(l)-a*(gList_(i)), gList);
			apply(i+1, j -> b_(l,j) = b_(l,j)-a*b_(i,j));
			a=leadCoefficient(gList_(l));
			); 
		    ));
	    if a != 0 then(
		gList = replace(l, (gList_(l))/a, gList);
		apply(l+1, j -> b_(l,j) = b_(l,j)/a);
		l = l+1;
		fList = append(fList, slnCasimir(n,d,S,fList_(l-1)));
		gList = append(gList, fList_(l));
		);
	    );
    );
    if gList_(l) == 0 then ( if b_(l,0) != 0 then Rf = 0 else Rf = ( sum for j from 1 to l list b_(l,j)*fList_(j-1) ) / b_(l,1); );
Rf
)

--******** Example 4.5.21 ********

slnReynolds(2,2,S,f1)


-- Documentation of the permutation groups methods
installPackage "OrbitSum" -- runs all checks
viewHelp "OrbitSum" -- opens documentation in browser
