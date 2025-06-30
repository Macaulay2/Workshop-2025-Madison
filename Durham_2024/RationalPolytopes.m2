newPackage(
    "RationalPolytopes",
    Version => "0.1",
    Date => "",
    Headline => "A package for Ehrhart theory of rational polytopes",
    Authors => {
        {Name => "Oliver Clarke", Email => "oliver.clarke@ed.ac.uk", HomePage => "https://www.oliverclarkemath.com/"},
        {Name => "Alex Milner", Email => "A.J.C.Milner@sms.ed.ac.uk", HomePage => ""},
        {Name => "Victoria Schleis", Email => "victoria.m.schleis@durham.ac.uk", HomePage => "https://victoriaschleis.github.io/"},
        {Name => "Vincenzo Reda", Email => "redav@tcd.ie", HomePage => ""},
        {Name => "Benoît Guerville-Ballé", Email => "benoit.guerville-balle@math.cnrs.fr", HomePage => "https://www.benoit-guervilleballe.com"}
        },
    AuxiliaryFiles => false,
    DebuggingMode => false,
    PackageExports => {"Polyhedra", "Normaliz"}
    )

export {
    "hStar",
    "Ehrhart",
    "EhrhartQP",
    "isPeriod",
    "quasiPolynomial",
    "period",
    "displayQP",
    "coefficientMonomial",
    "ehrhartSeries",
    "ReturnDenominator",
    "Backend",
		"latticePointsFromHData"
    }


-- define the denominator of polyhedron to be the lcm of all denominators of
-- coordinates of its vertices
denominator Polyhedron := P -> (
    if not P.cache#?"denominator" then (
        P.cache#"denominator" = lcm for j in flatten entries vertices P list denominator promote(j,QQ);
        );
    P.cache#"denominator"
    )

-* QuasiPolynomial Type *-


-- isPeriod(M, q)
-- M: Matrix
-- q: ZZ
-- returns true if for every j:
-- M^{j*q .. (j+1)*q-1} == M^{0 .. q-1}
isPeriod = method()
isPeriod(Matrix, ZZ) := (M, q) -> (
    if numRows M % q != 0 then false
		else (
        result := for j from 0 to numRows M // q - 1 do (
            if M^{0 .. q-1} != M^{j*q .. (j+1)*q-1} then break false;
            );
				if result === null then true else false
        )
    )


-- collapseMatrix M
-- M: Matrix
-- returns the smallest matrix that represents the
-- same quasi-polynomial as M
collapseMatrix = method()
collapseMatrix Matrix := M -> (
    q := for p from 1 to (numRows M)//2 do (
        if isPeriod(M, p) then break p;
        );
		if q =!= null then (
        M^{0 .. q-1}
        ) else M
    )

-- Definition of the Type QuasiPolynomial

QuasiPolynomial = new Type of HashTable

quasiPolynomial = method()
quasiPolynomial Matrix := M -> (
    MCollapsed:= collapseMatrix M;
    new QuasiPolynomial from {
        period => numRows MCollapsed,
        coefficients => MCollapsed,
        cache => new CacheTable,
        }
    )

net QuasiPolynomial := QP -> (
    description := "QuasiPolynomial of degree " | net(numColumns(QP#coefficients)-1) | " and of period " | net(QP#period);
    display := net displayQP(QP, Truncate => 5);
    stack {description, display}
    )


quasiPolynomial List := L -> (
    if not isMember(false, for l in L list instance(l,List)) then (
        D:=max for p in L list length p;
        L1:=for p in L list ((for i in 0..D-length p -1 list 0)|p);
        quasiPolynomial(matrix(L1))
        )
    else if not isMember(false, for l in L list instance(class l,PolynomialRing)) then(
        if not isMember(false, for l in L list numgens class l==1) then (
            D1:=max for p in L list (degree p)#0;
            lM:=for p in L list sub( (coefficients(p, Monomials=>for d in 0..D1 list ((generators class p)#0)^(D1-d)))#1, QQ);
            M:=transpose fold((a,b) -> a|b, lM);
            quasiPolynomial(M)
            )
        )
    )

-- note that we can borrow === function to determine equality of
-- QuasiPolynomials because we collapse the coefficients matrix
QuasiPolynomial == QuasiPolynomial := (QP1, QP2) -> QP1 === QP2


-- QuasiPolynomial as a function.
QuasiPolynomial ZZ := (QP, v) -> (
    internalQuasiPolynomial(QP,v)
    )

internalQuasiPolynomial = method()
internalQuasiPolynomial(QuasiPolynomial, ZZ) := (QP,t) -> (
    r := (QP#coefficients)^{t%QP#period};
    T := matrix for i in 0..(numColumns QP#coefficients - 1) list {t^(numColumns QP#coefficients - i - 1)};
    (r*T)_(0,0)
    )

-- Various methods associated to a QuasiPolynomial

-- display QuasiPolynomial in a friendly way
displayQP = method(
    Options => {
        Truncate => null -- null or integer
        }
    )
displayQP(QuasiPolynomial) := opts -> QP -> (
    R := monoid[getSymbol "t"];
    t := R_0;
    if opts.Truncate === null or QP#period <= opts.Truncate or opts.Truncate < 3 then (
        sum for d in 0 .. (numColumns QP#coefficients)-1 list (
            -- the leading coefficient is a constant so we could use this instead:
            --(if d > 0 then (QP#coefficients)_{d} else (QP#coefficients)_{d}^{0}
            --    ) expression if d < (numColumns QP#coefficients)-1 then (
            --    t^(numColumns (QP#coefficients)-d-1)

            (QP#coefficients)_{d} expression if d < (numColumns QP#coefficients)-1 then (
                t^(numColumns (QP#coefficients)-d-1)
                )
            else ""
            )
        )
    else (
        sum for d in 0 .. (numColumns QP#coefficients)-1 list (
            M := (QP#coefficients)_{d};
            Mcut := M^{0 .. opts.Truncate - 3, numRows QP#coefficients - 1};
            Mnet := net Mcut;
            netRows := (for i from 0 to opts.Truncate - 3 list Mnet#i) | {
                "| " | concatenate(width Mnet - 4: ".") | " |"} | {Mnet#(-1)};
            (stack netRows) expression if d < (numColumns QP#coefficients)-1 then (
                t^(numColumns (QP#coefficients)-d-1)
                )
            else ""
            )
        )
    )


degree QuasiPolynomial := QP -> (
    numColumns(QP#coefficients)-1
    )

period = method()
period QuasiPolynomial := QP -> (
    QP#period
    )

coefficients QuasiPolynomial := QP -> (
		QP#coefficients
		)

-- coefficientMonomial(QP, i)
-- QP: QuasiPolynomial
-- i : ZZ
-- returns the coefficient of t^i in the QuasiPolynomial
coefficientMonomial = method()
coefficientMonomial(QuasiPolynomial,ZZ) := (QP,i) -> (
    if i < degree(QP)+1 then M:=QP#coefficients_{degree(QP)-i};
    if i > degree(QP) then M=0;
    M
    )

-* Ehrhart Polynomial part *-

Ehrhart=method(TypicalValue=>RingElement)
Ehrhart (Polyhedron,ZZ):=(P, i) -> (
    n:=dim P;
    k:=denominator P;
    R:=QQ[getSymbol "x"];
    x:=R_"x";
    S:=for j from 0 to n list i+j*k;
    if n==0 and (not isEmpty P) then return 1+0*x;
    if isEmpty P then return 0+0*x;
    v:=matrix apply(S,h->(
            if h == 0 then {0}
            else {-1+#latticePoints(h*P)}
            )
        );
    v=promote(v,QQ);
    M:=promote(matrix apply(S,a->reverse apply(n+1,j->( a^j ))),QQ);
    M=flatten entries((inverse M)*v);


    1+sum apply(n+1,a->M_(a)*x^(n-a))
    )


EhrhartQP = method(
    Options => {
        Backend => "Normaliz" -- Normaliz or M2
        }
    )

EhrhartQP Polyhedron := opts -> P -> (
    if not P#cache#?"EhrhartQP" then (
        QP := if opts.Backend == "Normaliz" then (
            EhrhartQPNormaliz P
            )
        else if opts.Backend == "M2" then (
            EhrhartQPM2 P
            )
        else error("Unknown Backend");
        QP#cache#"OriginalPolyhedron" = P;
        P#cache#"EhrhartQP" = QP;
        );
    P#cache#"EhrhartQP"
    )

EhrhartQPM2 = method()
EhrhartQPM2 Polyhedron := P -> (
    k := denominator P;
    quasiPolynomial(for i from 0 to k-1 list Ehrhart(P,i))
    )

-- use the EhrharhtSeries of P (computed with Normaliz)
-- to construct the Quasi-Polynomial

EhrhartQPNormaliz = method()
EhrhartQPNormaliz Polyhedron := P -> (
    ES := ehrhartSeries(P, Backend => "Normaliz");
    R := ring ES;
    t := R_0;
    n := dim P;
    k := denominator P;
    latticePointCounts := for i from 0 to (n+1)*k -1 list (
        numberOfLatticePoints := ES(0);
        ES = (ES - numberOfLatticePoints)/t;
        numberOfLatticePoints
        );
    R' := QQ[getSymbol "x"];
    x := R'_0;
    QuasiPolyList := for i from 0 to k-1 list (
        S:=for j from 0 to n list i+j*k;
        if n==0 and (not isEmpty P) then return 1+0*x;
        if isEmpty P then return 0+0*x;
        v := matrix(QQ, apply(S,h -> {-1+latticePointCounts_h}));
        M := matrix(QQ, apply(S,a -> reverse apply(n+1,j ->  a^j )));
        M = flatten entries((inverse M)*v);
        1 + sum apply(n+1 , a -> M_a*x^(n-a))
        );
    quasiPolynomial(QuasiPolyList)
    )



-- hStar currently uses M2 code to compute
-- maybe we should add a Normaliz version too

-- For the M2 and Normaliz versions have different denominators - can they
-- be uniformised?
hStar = method(
    Options => {
        ReturnDenominator => false, --returns a pair of polys (h, d) s.t. Ehrhart series is h/d
        Backend => "Normaliz" -- either Normaliz or M2
        })

hStar(Polyhedron, Ring) := opts -> (P, R) -> (
    if numgens R < 1 then error("ring must have at least one generator");
    if not P#cache#?"EhrhartSeriesNumerator" then (
        if opts.Backend == "M2" then (
            hStarM2(P, R);
            )
        else if opts.Backend == "Normaliz" then (
            hStarNormaliz(P, R);
            )
        else error("unknown Backend option: " | toString opts.Backend);
        );
    if opts.ReturnDenominator then (
        P#cache#"EhrhartSeriesNumerator",
        P#cache#"EhrhartSeriesDenominator"
        )
    else P#cache#"EhrhartSeriesNumerator"
    )

hStar(Polyhedron) := opts -> P -> (
    R:=QQ[getSymbol "t"]; -- potentially redundant if hStar has already been computed
    hStar(P, R, opts)
    )

-- M2 version of hStar polynomial
-- once computed, it updates the cache
hStarM2 = method()
hStarM2(Polyhedron, Ring) := (P, R) -> (
    n:=dim P;
    dnom := lcm for i in flatten entries vertices P list denominator promote(i,QQ);
    p:=1;
    t:=R_0;
    for i from 1 to (n+1)*dnom do (p=p + #latticePoints(i*P) * t^i);
    r:=(1-t^dnom)^(n+1);
    f := (p*r) % t^((n+1)*dnom);
    P#cache#"EhrhartSeriesNumerator" = f;
    P#cache#"EhrhartSeriesDenominator" = r;
    (f, r)
    )


-- Normaliz version of hStar polynomial
-- once computed, it updates the cache

-- WARNING: the denominator of the hstar might not be the usual one!
hStarNormaliz = method()
hStarNormaliz(Polyhedron, Ring) := (P, R) -> (
    t := R_0;
    C := normaliz(transpose vertices P, "polytope"); -- Maybe all of this data can be stored for later use
    numeratorCoefficients := C#"inv"#"hilbert series num";
    denominatorFactors := C#"inv"#"hilbert series denom";
    f := sum for i from 0 to #numeratorCoefficients -1 list (numeratorCoefficients#i) * t^i;
    r := product for i from 0 to #denominatorFactors -1 list 1 - t^(denominatorFactors#i);
    P#cache#"EhrhartSeriesNumerator" = f;
    P#cache#"EhrhartSeriesDenominator" = r;
    (f, r)
    )


ehrhartSeries = method(
    Options => {
        Backend => "Normaliz" -- Normaliz or M2
        }
    )

ehrhartSeries(Polyhedron, Ring) := opts -> (P, R) -> (
    if not P#cache#?"EhrhartSeries" then (
        (h, d) := hStar(P, R, ReturnDenominator => true, Backend => opts.Backend);
        R' := ring h; -- if R' =!= R then we previously constructed R' so we should ignore R
        F := frac R';
        h = h_F;
        d = d_F;
        P#cache#"EhrhartSeries" = h/d;
        );
    P#cache#"EhrhartSeries"
    )

ehrhartSeries Polyhedron := opts -> P -> (
    ehrhartSeries(P, QQ[getSymbol "t"], opts)
    )

---------------------------------------
---------------------------------------
-- Temporary version of Normaliz
-- to use with computing rational polytopes
--------------------------------------

debug Normaliz

-- writes the given data in a normaliz input file
-- doWriteNmzData = method()
-- writes several matrices in a normaliz input file
doWriteNmzData List := matrices -> (
    checkNmzFile("doWriteNmzData");
    outf := nmzFile | ".in" << "";
    for p in matrices do (
        sgr := p#0;
        nmzMode := p#1;
        outf << numRows sgr << endl;
        outf << numColumns sgr << endl;
        if ring sgr =!= ZZ and ring sgr =!= QQ then error("matrix with non-rational entries");
        for i from 0 to numRows sgr - 1 do (
            s:= "";
            for j from 0 to numColumns sgr - 1
            do s = s | toString(sgr_(i,j)) | " "; -- MODIFIED: this handles ZZ and QQ entries
            outf << s << endl;
            );
        --Until version 3.9.4, input type normal_toric_ideal was called lattice_ideal
        if normalizProgram#"version" < "3.10" and nmzMode == "normal_toric_ideal" then nmzMode = "lattice_ideal";
        outf << nmzMode << endl);
    outf << close
    )


-- lattice points of a rational polytope using Normaliz
-- returns a matrix whose columns are the lattice points
latticePointsFromHData = method()
latticePointsFromHData(Matrix, Matrix) := (I, v) -> (
		-- polytope is given by Ix <= v
		M := -I | v;
		normalizOutput := normaliz(M, "inhom_inequalities");
		n := numColumns normalizOutput#"gen";
		transpose normalizOutput#"gen"_{0 .. n-2}
		)

latticePointsFromHData(Matrix, Matrix, Matrix, Matrix) := (I, v, E, w) -> (
		-- polytope is given by Ix <= v, Ex = w
		M   := -I |  v;
		M'  := -E |  w;
		M'' :=  E | -w;
		normalizOutput := normaliz(M || M' || M'', "inhom_inequalities");
		n := numColumns normalizOutput#"gen";
		transpose normalizOutput#"gen"_{0 .. n-2}
		)


---------------------------------------
-* Documentation section *-
beginDocumentation()

doc ///
  Key
    RationalPolytopes
  Headline
    A package for Ehrhart theory of rational polytopes
///

doc ///
  Key
    Ehrhart
  Headline
    a function
  Usage
    p = Ehrhart(P,i)
  Inputs
    P : Polyhedron
    i : ZZ
      (P,i) is a polyhedron and integer to calculate the ith polynomial piece contributing to the Ehrhart quasipolynomial for P
  Outputs
    p : QQ[x]
      the ith polynomial contributing to the Ehrhart quasipolynomial of P
  Description
    Text
      it calculates the ith polynomial piece contributing to the Ehrhart quasipolynomial of a polyhedron P
    Example
      Ehrhart(convexHull transpose matrix "0,0;1/2,0;0,1/2",0)
      Ehrhart(convexHull transpose matrix "0,0;1/2,0;0,1/2",1)

  SeeAlso
    RationalPolytopes
///


doc ///
  Key
    EhrhartQP
  Headline
    a function
  Usage
    L = EhrhartQP(P)
  Inputs
    P : Polyhedron
      the polyhedron for which we calculate the Ehrhart quasipolynomial
  Outputs
    L : List
      a list of polynomial pieces contributing to the Ehrhart quasipolynomial of P
  Description
    Text
      it calculates the Ehrhart quasipolynomial of polyhedron P
    Example
      EhrhartQP(convexHull transpose matrix "0,0;1/2,0;0,1/2")
///



doc ///
  Key
    hStar
    (hStar, Polyhedron)
    (hStar, Polyhedron, Ring)
    [hStar, Backend]
    [hStar, ReturnDenominator]
  Headline
    the $h^*$-polynomial of a polytope
  Usage
    hStar P
    hStar(P,R)
  Inputs
    P : Polyhedron
      A convex polyhedron which must be compact
    R : Ring
      A ring in at least one variable
    ReturnDenominator => Boolean
      whether to return the denominator of the Ehrhart series
    Backend => String
      either "Normaliz" or "M2", selects the method for computing
      the Ehrhart series
  Outputs
    f : RingElement
      the $h^*$-polynomial (in the ring R) of $P$
  Description
    Text
      If $P$ is a lattice polytope then this function returns the
      $h^*$-polynomial of $P$.
    Example
      hStar convexHull transpose matrix "-1,0; 0,-1; 1,0; 0,1"
    Text
      Otherwise if $P$ is a rational polytope then the behavior
      may give unexpected results because the calculation of
      the Ehrhart series (see @TO ehrhartSeries@) may use
      a denominator that is different to the expected one.
    Example
      hStar(convexHull transpose matrix "0; 1/2", Backend => "Normaliz")
      hStar(convexHull transpose matrix "0; 1/2", Backend => "M2")
    Text
      To return the denominator of the Ehrhart series, set the
      optional argument @TO ReturnDenominator@ to @TT "true"@.
      In this case, the result is a pair that forms the numerator
      and denominator of the Ehrhart series.
    Example
      hStar(convexHull transpose matrix "0; 1/2",
          Backend => "Normaliz", ReturnDenominator => true)
      hStar(convexHull transpose matrix "0; 1/2",
          Backend => "M2", ReturnDenominator => true)
  SeeAlso
    RationalPolytopes
    ehrhartSeries
///

doc ///
  Key
    isPeriod
  Headline
    a function
  Usage
    result = isPeriod(M,q)
  Inputs
    M : Matrix
    q : ZZ
  Outputs
    result : Boolean
      result is true if q is the period and false otherwise
  Description
    Text
      It returns true if q is the period and false otherwise
    Example
      isPeriod(matrix "1,1;2,2;1,1;2,2",2)
      isPeriod(matrix "1,1;2,2;1,1;3,2",2)
      isPeriod(matrix "1,2,5;2,4,3;1,2,5;3,4,5;1,2,5;2,4,3;1,2,5;3,4,5",3)
  SeeAlso
    RationalPolytopes
///

doc ///
  Key
    period
  Headline
    a function
  Usage
    p = period(QP)
  Inputs
    QP : QuasiPolynomial
      A quasipolynomial of which we want to compute the period
  Outputs
    p : ZZ
      The period of QP
  Description
    Text
      Computes the period of a quasipolynomial
    Example
      period(quasiPolynomial(matrix{{1,2,3},{1,4,5}}))
      period(quasiPolynomial(matrix{{1,1},{1,2},{1,1},{1,2}}))
      period(quasiPolynomial(matrix{{3,6,7,2},{3,4,4,2},{3,2,5,6}}))
  SeeAlso
    RationalPolytopes
///

doc ///
  Key
    coefficientMonomial
  Headline
    a function
  Usage
    L = coefficientMonomial(QP,i)
  Inputs
    QP : QuasiPolynomial
      A quasipolynomial of which we want to know the coefficients
    i : ZZ
      The degree of the monomials of QP of which we want to know the coefficients
  Outputs
    L : List
      The coefficients of the monomials of degree i appearing in QP
  Description
    Text
      Computes the coefficients of the monomials of degree i appearing in QP
    Example
      coefficientMonomial(quasiPolynomial(matrix{{1,2,3},{1,4,5}}),2)
      coefficientMonomial(quasiPolynomial(matrix{{1,1},{1,2},{1,1},{1,2}}),0)
      coefficientMonomial(quasiPolynomial(matrix{{3,6,7,2},{3,4,4,2},{3,2,5,6}}),2)
  SeeAlso
    RationalPolytopes
///



-* Test section *-
TEST /// -- (1)
R=QQ[t]
assert(1_R == hStar(convexHull transpose matrix "0,0,0;1,0,0;0,1,0;0,0,1",R))
assert(t^5+3*t^4+4*t^3+4*t^2+3*t+1 == hStar(convexHull transpose matrix "1,0;-1,0;0,1/2;0,-1/2",R))
assert(t+1 == hStar(convexHull transpose matrix "0; 1/2",R))
assert(t^5+t^3+t^2+1 == hStar(convexHull transpose matrix "1/4; 1/2",R))
///



end--

----------------------------------
-- Plans for future development --
----------------------------------

-- To-do list --

-- check exported functions work with easy examples
-- that can be computed by hand


-- implement a method for internalQuasiPolynomial that implements the following procedure:
-- 1) check the cache for a stored list of polynomials
-- 2) if there is no list, use the coefficients matrix to produce a list of polynomials and cache them
-- 3) take the input i modulo the period to obtain j, and return polynomial number j evaluated at i


-- decide what should be done if we try to create a quasi polynomial of period 1.
-- it's just a polynomial! So should we return a genuine polynomial or not?


-- if a quasi polynomial is made from a polytope, then store a reference to that polytope in the cache of
-- the quasi polynomial


-- cache the quasi-polynomial in the polyhedron and avoid recomputing the quasi-polynomial if it is already cached
-- note that the Polyhedron type is just a hashtable with a single entry: cache


-- check the definition of hStar polynomial in literature and check whether the denominator of the Ehrhart series is:
-- (1 - t^(denominator P))^(dim P)  or  (1 - t^(period P))^(dim P)
-- Answer:
--   in the literature, authors typically do not define a denominator / hstar polynomial for rational polytopes
--


-- simplify the names of functions: E.g. EhrhartQP -> ehrhart (overriding the one in Polyhedra)
-- or, if we don't want to override, then choose a name without abbreviations: e.g. ehrhartQuasiPolynomial
-- function names and variables should start with lower case
-- periodQP -> period (may need to change the key in the QuasiPolynomial type)


-- think about how a user might interact with the package and what would make life easier for them.
-- E.g. A user comes along with a polytope in mind: either they know the vertices or a half-space description
--      the user want to compute the Ehrhart quasi-polynomial, Ehrhart series, hStar poly, delta-vector (coefficients of hStar poly)


-- Whenever we perform a computation, e.g. computing the EhrhartQP, store the result in the cache
-- and before performing computations, check if we have already computed it by checking the cache
-- a useful piece of code is:
C = new CacheTable from {1 => "hi"}
C#?1 -- 1 is a key of the hash table
C#?2 -- but 2 is not


--------------------------------------------------

-* Development section *-
restart

uninstallPackage "RationalPolytopes"

restart
installPackage "RationalPolytopes"

viewHelp "RationalPolytopes"
debug needsPackage "RationalPolytopes"

check "RationalPolytopes"


P=convexHull transpose matrix "0;1/2"
EQP=EhrhartQP(P)

P=convexHull transpose matrix "1,0;-1,0;0,1/2;0,-1/2"
displayQP EhrhartQP(P)

P=convexHull transpose matrix "-1/2; 1/2"
EhrhartQP(P)



-- Test of the constructor of the Type QuasiPolynomial


S=matrix({{1,2,3},{1,4,5}})
M=matrix({{1,2,3},{0,1,0},{1,2,3},{0,1,0}})

QP=quasiPolynomial(M)
QP#period
print QP
displayQP QP
degree QP
periodQP QP
coefficientMonomial(QP,0)

R=QQ[x]
R1=QQ[t]
L={x^2+3,2*t}
quasiPolynomial(L)

L={{1,0,3},{2,0}}
quasiPolynomial(L)

-----------------------------
restart
needsPackage "RationalPolytopes"

debug RationalPolytopes

P = convexHull transpose matrix "1,0;-1,0;0,1/20;-1,11/20"
latticePoints P
EhrhartQP(P)
displayQP(EhrhartQP P, Truncate => 1000)

QP1 = EhrhartQPM2(P)
QP2 = EhrhartQPNormaliz(P)

QP1 === QP2

hStar P
ehrhartSeries P

hStar(P, Backend => "M2")
ehrhartSeries(P, Backend => "M2")

f = EhrhartQP(P)
displayQP(f)


-------------------

P = convexHull transpose matrix "0,0;0,1;1,0;1,1"
f = EhrhartQP(P)
displayQP(f)
f 100

C = normaliz(transpose vertices P, "polytope")


-------------------------

restart
needsPackage "RationalPolytopes"

-- Ax + v >= 0
-- M = [A | v]
M = matrix "-8,2,500; 1, -1, 0; 2, 7, 3"
s = "inhom_inequalities"

debug Normaliz
elapsedTime runNormaliz(M, s);
oo#"gen" -- list of lattice points as row vectors of a matrix + 1


P = polyhedronFromHData(- matrix "-8, 2; 1, -1; 2, 7", matrix "500; 0; 3")
elapsedTime (latticePoints P);


setNmzFile()
writeNmzData(M, s)
-- runNormaliz(M, s)

-- created a temp file with data
-- manually edited it to say inhom_inequalities + LatticePoints


dir := select(".*/", nmzFile);
runDir := if dir != {} then dir#0 else null;
elapsedTime runProgram(normalizProgram, getNmzExec(), collectNmzOptions() | baseFilename nmzFile,
		RunDirectory => runDir, Verbose => debugLevel > 0);

    -- return nothing if .gen is not generated
    if not nmzGen then ( if nmzFilename == "" then rmNmzFiles(); return );

    if not opts.allComputations then (
	nmzData := readNmzData "gen";
	rc := new RationalCone from { "gen" => nmzData, "inv" => getNumInvs() };
	if nmzFilename == "" then rmNmzFiles();
	return rc);

    -- read all files written
    files := { "inv" => getNumInvs() };
    suffixes := { "gen","egn","esp","tri","typ","ht1","ext","tgn" };
    for s in suffixes do if fileExists(nmzFile | "." | s) then files = append(files, s => readNmzData s);

    L := readMultipleNmzData "cst";
    files = append(files, "sup" => L#0);
    files = append(files, "equ" => L#1);
    files = append(files, "cgr" => L#2);

    C := new RationalCone from files;

    if nmzFilename == "" then rmNmzFiles();
    C

-- y <= x, x >= 0, y = 0
latticePointsFromHData(matrix "-1, 1;1, 0", matrix "0; 10", matrix "0, 1", matrix "0")

-- strange output if the polyhedron is unbounded
latticePointsFromHData(matrix "-1; 1", matrix "-1; 10")
