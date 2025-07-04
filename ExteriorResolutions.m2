newPackage("ExteriorResolutions",
    Version => "1.1",
    Date => "24 June 2025",
    Headline => "Injective resolutions over exterior algebras",
    Authors => {
	{Name => "Penelope Beall",       Email => "pbeall@ucdavis.edu",          HomePage => "https://pbeall.github.io" },
	{Name => "Michael K. Brown",     Email => "mkb0096@auburn.edu",          HomePage => "http://webhome.auburn.edu/~mkb0096/" },
	{Name => "Caitlin M. Davis",     Email => "cmdavis22@wisc.edu",          HomePage => "https://sites.google.com/wisc.edu/caitlindavis/home" },
	{Name => "Andrew Karsten",       Email => "akk0071@auburn.edu",          HomePage => "https://www.auburn.edu/cosam/departments/math/students/grad/graduate-students.htm" },
	{Name => "Jiaheng Li",           Email => "jli3440@gatech.edu",          HomePage => "" },
	{Name => "Jianuo Zhou",          Email => "jzhou632@gatech.edu",         HomePage => "https://math.gatech.edu/people/jianuo-zhou"},
	{Name => "Boyana Martinova",     Email => "boyana.martinova@gmail.com",  HomePage => "https://sites.google.com/view/bmartinova/home"},
	{Name => "Mahrud Sayrafi",       Email => "mahrud@umn.edu",              HomePage => "https://mahrud.github.io/" },
	{Name => "Gregory G. Smith",     Email => "ggsmith@mast.queensu.ca",     HomePage => "https://mast.queensu.ca/~ggsmith/" },
	{Name => "Sreehari Suresh Babu", Email => "sreeharisbabu183@gmail.com",  HomePage => "https://sreehari183.github.io/" }
	},
    Keywords => {"Commutative Algebra"},
    PackageExports => { "Complexes", "SimplicialComplexes" },
    AuxiliaryFiles => true,
    DebuggingMode  => true
    )

export {
    --methods
    "injectiveResolution",
    "injectiveResolutionMap",
    "coaugmentationMap",
    "priddyComplex",
    "priddyDifferential",
    "koszulRR",
    "koszulLL",
    "koszulPair",
    "koszulDual",
    "exteriorStanleyReisner",
}

--------------------------------------------------

-- TODO: move to Core, c.f. https://github.com/Macaulay2/M2/issues/3844
-- assume 'test' is a monotonic function, i.e. false for all i < n then true for i >= n
binarySearch = method()
-- return the first index of an element in L such that test(L#i) is true
binarySearch(List,        Function) := (L,          test) -> binarySearch(0, #L-1, i -> test(L#i))
-- shorthand for search in [0, n)
binarySearch(         ZZ, Function) := (      high, test) -> binarySearch(0, high-1, test)
-- shorthand for when a lower bound isn't known
binarySearch(Nothing, ZZ, Function) := (null, high, test) -> (
    dist := 1;
    while true do if test(high - dist)
    then (high, dist) = (high - dist, dist * 2)
    else break binarySearch(high - dist + 1, high, test));
-- shorthand for when an upper bound isn't known
binarySearch(ZZ, Nothing, Function) := (low, null, test) -> (
    dist := 1;
    while true do if test(low + dist)
    then break binarySearch(low, low + dist, test)
    else (low, dist) = (low + dist + 1, dist * 2))
-- standard binary search
binarySearch(ZZ, ZZ, Function) := (low, high, test) -> (
    -- TODO: are the first two lines standard?
    if     test(low)  then return low;
    --if not test(high) then return high + 1;
    while high - low > 1 do (
	mid := (high + low) // 2;
	if test(mid) then high = mid else low = mid);
    high)

-- TODO: move to Complexes
hilbertPolynomial Complex := o -> C -> sum(pairs C.module, (i, M) -> (-1)^i * hilbertPolynomial(M, o))

--------------------------------------------------
--- Injective resolutions
--------------------------------------------------

injectiveResolution = method(Options => options freeResolution)
injectiveResolution Module := Complex => opts -> M -> M.cache.injectiveResolution ??= (
    E := ring M;
    if not isSkewCommutative E then error "expected underlying ring to be skew-commutative";
    P := Hom(freeResolution(Hom(M, E), opts), E);
    P.cache.injectiveResolution = M;
    P)

injectiveResolution Complex := Complex => opts -> C -> C.cache.injectiveResolution ??= (
    E := ring C;
    if not isSkewCommutative E then error "expected underlying ring to skew-commutative";
    D := Hom(freeResolution(Hom(C, E), opts), E);
    D.cache.injectiveResolution = C;
    D
    )


injectiveResolutionMap = method(Options => options freeResolution)
injectiveResolutionMap Module := ComplexMap => opts -> M -> (
    C := injectiveResolution(M, opts);
    map(C, complex M, i -> if i === 0 then map(C_0, M, transpose syz transpose presentation M))
    )
injectiveResolutionMap Complex := ComplexMap => opts -> C -> (
    D := injectiveResolution(C, opts);
    W := ring C;
    CDoubleDual:= Hom(Hom(C,W^1),W^1);
    DDualMap := map(CDoubleDual, C, i -> (
	    gensMatrix := gens C_i;
	    h :=map(C_i, source gensMatrix, id_(source gensMatrix));
	    ddh := Hom(Hom(h,W^1),W^1);
	    map(CDoubleDual_i, C_i, matrix ddh)
	    )
	);
    Hom(resolutionMap(Hom(C,W^1),opts),W^1)* DDualMap
    --Found this code in Divisor package, seems to work, need to cache hfC in injRes --Sreehari
    --(lo,hi) := concentration C;
    --tempHash := new MutableHashTable;
    --tempHash#hi = map(D_hi, C_hi, transpose syz transpose presentation C_hi);
    --for i in reverse toList(lo..(hi-1)) do (
        --tempHash#i = (dd^C_(i + 1))\\(dd^D_(i+1) * tempHash#(i + 1));
    --);
   -- map(D, C, tempHash)
    -*
    map(D, C, i -> (
	    if isFreeModule C_i then map(D_i,C_i, id_(C_i))
	    else map(D_i, C_i, transpose syz transpose presentation C_i)
	    )
	)
    *-
    )

coaugmentationMap = method()
coaugmentationMap Complex := ComplexMap => P ->  (
    if not P.cache.?injectiveResolution then
	error "expected input to be constructed as an injective resolution";
    C := P.cache.injectiveResolution;
    injectiveResolutionMap C
    )

------------------------------------------------
--- Priddy complex
--------------------------------------------------

priddyDifferential = method(TypicalValue => Matrix)
priddyDifferential(ZZ, Matrix) := (i, m) -> m.cache#(symbol priddyDifferential, i) ??= (
    R := ring m;
    n := numgens R - 1;
    -- TODO: what's the right way to make this step uniform?
    (S, E) := koszulPair(numcols m - 1, coefficientRing R,
	Variables => {x := symbol x, e := symbol e});
    A := if isSkewCommutative R then S else E;
    --L := flatten (degrees m)_1; --degrees of the forms on which we are taking the Priddy complex
    --assert(all(L, i -> odd L_i));
    srcmons := basis(-i, A);
    tarmons := basis(-i+1, A);
    --
    srcexps := apply(first entries srcmons, mon -> flatten exponents mon);
    tarexps := apply(first entries tarmons, mon -> flatten exponents mon);
    --
    -- TODO: what's the right way to make this step uniform?
    W := if isSkewCommutative A then R^1 else R^{n+1};
    --
    src := if srcmons == 0 then R^0 else directSum apply(srcexps,
	e -> W ** R^{sum(numcols m, i -> e_i * flatten degree(m_i))});
    tar := if tarmons == 0 then R^0 else directSum apply(tarexps,
	e -> W ** R^{sum(numcols m, i -> e_i * flatten degree(m_i))});
    --
    -- TODO: can this work on the whole complex at once rather than term by term?
    f := matrix table(numcols tarmons, numcols srcmons,
	(r,c) -> tarmons_{r} // srcmons_{c});
    -- TODO: can we avoid sub?
    map(tar, src, if f == 0 then 0 else sub(f, m)))

priddyComplex = method(TypicalValue => Complex, Options => { LengthLimit => null })
priddyComplex Matrix := opts -> m -> m.cache#(symbol priddyComplex, opts.LengthLimit) ??= (
    (lo, hi) := if instance(opts.LengthLimit, ZZ) and opts.LengthLimit > 0 then (-opts.LengthLimit, 0)
    else try degreeSupport module koszulDual ring m else error "expected positive length limit";
    complex hashTable apply(lo+1..hi, i -> i => priddyDifferential(i, m)))

--------------------------------------------------
--- Koszul duality helpers
--------------------------------------------------

-- returns a pair S = Sym^* K^(n+1) and E = Wedge^* K^(n+1)
koszulPair = method(Options => { Variables => {"x", "e"}})
koszulPair(ZZ, Ring) := opts -> (n, K) -> (
    x := if instance(opts.Variables#0, String) then getSymbol opts.Variables#0 else opts.Variables#0;
    e := if instance(opts.Variables#1, String) then getSymbol opts.Variables#1 else opts.Variables#1;
    S := K[x_0..x_n];
    E := K[e_0..e_n, SkewCommutative => true];
    S.cache.koszulDual = E;
    E.cache.koszulDual = S;
    (S, E))

-- given a polynomial ring returns an exterior algebra and vice versa
koszulDual = method(Options => { Variables => {"x", "e"}})
koszulDual Ring := opts -> A -> A.cache.koszulDual ??= (
    K := coefficientRing A;
    x := if instance(opts.Variables#0, String) then getSymbol opts.Variables#0 else opts.Variables#0;
    e := if instance(opts.Variables#1, String) then getSymbol opts.Variables#1 else opts.Variables#1;
    if isSkewCommutative A
    then K[x_0..x_(numgens A - 1)]
    else K[e_0..e_(numgens A - 1), SkewCommutative => true])

-- given a finite length module, finds (lo,hi)
-- such that M_d == 0 for d < lo and hi < d
degreeSupport = method()
degreeSupport Module := M -> if M == 0 then (0,0) else (
    if hilbertPolynomial M != 0 then error "expected a range provided as Concentration => (lo, hi)";
    - binarySearch(first max degrees M, , i -> hilbertFunction(i, M) == 0) + 1,
    - binarySearch(, first min degrees M, i -> hilbertFunction(i, M) != 0))
degreeSupport Complex := C -> if C == 0 then (0,0) else (
    if hilbertPolynomial C != 0 then error "expected a range provided as Concentration => (lo, hi)";
    supports := apply(pairs C.module, (i, M) -> {i,i} + toList degreeSupport M);
    supports / first // min,
    supports / last  // max)
degreeSupport ComplexMap := f -> if f == 0 then (0,0) else (
    supports := apply({source f, target f}, degreeSupport);
    supports / first // min,
    supports / last  // max)

--------------------------------------------------
--- Koszul duality functors
--------------------------------------------------

-- RR: Com(S) -> Com(E) is the right-adjoint functor
koszulRR = method(Options => { Concentration => null })

-- LL: Com(E) -> Com(S) is the left-adjoint functor
koszulLL = method(Options => options koszulRR)

-- the implementations of koszulRR and koszulLL are identical
-- except in one step, so we combine their implementations

koszulDualityFunctorModule = functor -> opts -> M -> M.cache#(functor, opts) ??= (
    try (lo, hi) := opts.Concentration
    else (lo, hi) = degreeSupport M;
    if hi < lo then error "expected nonempty concentration";

    A := ring M;
    n := numgens A - 1;
    A' := koszulDual A; -- this is A^!
    ev := map(A', A, vars A'); -- not a map of algebras, just substitutes variables

    -- TODO: what's the right way to make this step uniform?
    W := if isSkewCommutative A then A'^1 else A'^{n+1};
    modules := hashTable apply(lo..hi,
	i -> i => W ** A'^{-i} ** (A' ** part(-i, M)));

    if lo == hi
    then complex(modules#lo, Base => lo)
    else complex hashTable apply(lo+1..hi, i -> i => (
	    -- This is borrowed bgg in the BGG package
	    src := basis(-i, M);
	    tar := basis(-i+1, M);

	    -- we construct the differential by factoring it through (vars A) ** src
	    g := ((vars A) ** src) // tar; -- g is a map from source (vars A) ** src to source tar
	    -- which makes the traingle involving (vars A) ** src and tar commute
	    b := (ev g) * ((transpose vars A') ** (ev source src));

	    map(modules#(i-1), modules#i, b))
	)
    )

-- RR(M)^i = E^*(i) \otimes_k M_i
-- E^* = Hom_k(E, k) = E(n+1)
koszulRR Module := Complex => opts -> (koszulDualityFunctorModule koszulRR) opts
-- LL(N)^i = S(i) \otimes_k N_i
koszulLL Module := Complex => opts -> (koszulDualityFunctorModule koszulLL) opts

--------------------------------------------------

koszulDualityFunctorMatrix = functor -> opts -> f -> f.cache#(functor, opts) ??= (
    R := koszulDual ring f;
    src := functor(source f, opts);
    tar := functor(target f, opts);
    map(tar, src, i -> map(tar_i, src_i, R ** part(-i, f))))

-- RR(y**s) = \sum_{l=0}^n y*e_l ** s*x_l
koszulRR Matrix := ComplexMap => opts -> (koszulDualityFunctorMatrix koszulRR) opts
-- LL(s**y) = (-1)^i \sum_{l=0}^n x_l*s \otimes y*e_l
koszulLL Matrix := ComplexMap => opts -> (koszulDualityFunctorMatrix koszulLL) opts

--------------------------------------------------

koszulDualityFunctorComplex = functor -> opts -> C -> (
    try (lo, hi) := opts.Concentration -- bounds for homological degrees i
    else (lo, hi) = degreeSupport C;
    (inf, sup) := concentration C; -- bounds for homological degrees k in C

    terms := hashTable apply(inf..sup,
	k -> k => functor(C_k,    Concentration => (lo-k, hi-k)));
    diffs := hashTable apply((inf+1)..sup,
	k -> k => functor(C.dd_k, Concentration => (lo-k, hi-k)));

    modules := hashTable apply(lo..hi,
	i -> i => hashTable apply(inf..sup,
	    k -> k => (terms#k)_(-k+i)));

    if lo == hi then return complex(directSum values modules#lo, Base => lo);

    complex hashTable apply((lo+1)..hi,
	i -> i => matrix table(sup-inf+1, sup-inf+1,
	    (r,c) -> map(modules#(i-1)#(inf+r), modules#i#(inf+c),
		if r == c   then (-1)^r * dd^(terms#(inf+r))_(-inf-r+i) else
		if r == c-1 then             (diffs#(inf+c))_(-inf-c+i) else 0)
	    )
	)
    )

-- RR(C)^i = \bigoplus_{j\in\ZZ} Hom_k(E(-j), C^{i-j}_j)
--         = \bigoplus_{j\in\ZZ} (E(-j))^* \otimes_k C^{i-j}_j
--         = \bigoplus_{j\in\ZZ} E^*(j)    \otimes_k C^{i-j}_j
koszulRR Complex := Complex => opts -> (koszulDualityFunctorComplex koszulRR) opts
-- LL(D)^i = \bigoplus_{j\in\ZZ} S(j) \otimes_k D^{i-j}_j
koszulLL Complex := Complex => opts -> (koszulDualityFunctorComplex koszulLL) opts

--------------------------------------------------

koszulDualityFunctorComplexMap = functor -> opts -> f -> (
    src := functor(C := source f, opts);
    tar := functor(D := target f, opts);

    try (lo, hi) := opts.Concentration
    else (lo, hi) = degreeSupport f;

    concentrations := apply({C, D}, concentration);
    inf := concentrations / first // min;
    sup := concentrations / last  // max;

    maps := hashTable apply(inf..sup,
	k -> k => functor(f_k, Concentration => (lo-k, hi-k)));

    map(tar, src, hashTable apply(lo..hi,
	    i -> i => map(tar_i, src_i,
		matrix table(sup-inf+1, sup-inf+1,
		    (r,c) -> if r == c then (maps#(inf+r))_(-inf-r+i) else 0)
		)
	    )
	)
    )

-- RR(y**s) = (-1)^i RR(y**s) + y ** dd^C(s) for s \in C^{i-j}_j
koszulRR ComplexMap := ComplexMap => opts -> (koszulDualityFunctorComplexMap koszulRR) opts
-- LL(s**y) = (-1)^i LL(s**y) + s ** dd^D(y) for y \in D^{i-j}_j
koszulLL ComplexMap := ComplexMap => opts -> (koszulDualityFunctorComplexMap koszulLL) opts

--------------------------------------------------
--- Stanley Reisner
--------------------------------------------------

exteriorStanleyReisner = method()
exteriorStanleyReisner(SimplicialComplex, Ring) := (S, K) -> (

    n := dim S;

    R := ring S;
    E := K[gens R, SkewCommutative => true];

    L := {};
    for m in facets S do (
        L = append(L, sub(m, E));
    );

    ideal L

)

--------------------------------------------------
--- Documentation
--------------------------------------------------

beginDocumentation()


doc ///
Node
   Key
      ExteriorResolutions
   Headline
      Package that implements homological constructions over the exterior algebra
--   Description
--      Text
--   SeeAlso
--   References
--      Text

Node
   Key
       priddyComplex
      (priddyComplex, Matrix)
   Headline
       computes the Priddy complex of several elements of an exterior algebra
   Usage
       priddyComplex m
   Inputs
       m:Matrix
	   consisting of a single row and with entries in an exterior algebra
   Outputs
       :Complex
   Description
       Text
           Let $S = k[x_0, \dots, x_n]$ be a polynomial ring over a field.
	   Let $f_1, \dots, f_c$ be odd degree elements of an exterior algebra
	   $E = \Lambda_k(e_0, \dots, e_n)$. Then, we get a complex
	   $E\otimes_k S_0\rightarrow E\otimes_k S_1\rightarrow \cdots$
	   whose differentials are left multiplication by $\sum_{i=0}^n f_i\otimes x_i$.

           Below we obtain the Priddy complex where $n=2$ and $f_1=e_0$ and $f_2=e_0*e_1$ and $k=\mathbb{Q}$.
       Example
	   E = QQ[e_0,e_1, SkewCommutative=>true]
	   m = matrix{{e_0, e_0*e_1}}
	   C = priddyComplex(m, LengthLimit => 3)
	   C.dd
	   prune HH C
Node
   Key
     priddyDifferential
    (priddyDifferential, ZZ, Matrix)

Node
   Key
     injectiveResolution
Node
   Key
     coaugmentationMap

Node
   Key
     koszulRR
Node
   Key
     koszulLL
Node
   Key
       exteriorStanleyReisner
      (exteriorStanleyReisner, SimplicialComplex, Ring)
   Headline
       computes the exterior algebra version of the Stanley Reisner face ideal of a simplicial complex
   Usage
       exteriorStanleyReisner(S, K)
   Inputs
       S:SimplicialComplex
	   given as an ideal in a polynomial ring
   Outputs
       :Ideal
   Description
       Text
           The exterior algebra version of the Stanley Reisner face ideal is defined using the same generators as the ordinary Stanley Reisner face ideal, but over an exterior algebra instead of a polynomial ring.
///

--------------------------------------------------
--- Tests
--------------------------------------------------

-*
XXX
restart
needsPackage "ExteriorResolutions"
check ExteriorResolutions
*-
TEST ///
    E = ZZ/101[e_0..e_3,SkewCommutative => true]
    M = E^2
    Res = injectiveResolution(M, LengthLimit => 4)
    assert(Res == complex M)
    assert isWellDefined Res
    f = coaugmentationMap Res
    assert isWellDefined f
    assert isQuasiIsomorphism f
    assert(f == injectiveResolutionMap M)

    I = ideal vars E
    k = E^1/I
    P = injectiveResolution(k, LengthLimit => 10)
    assert isWellDefined P
    assert isFree P
    f = coaugmentationMap P
    assert isWellDefined f
    assert isQuasiIsomorphism f
    assert isComplexMorphism f
    S = ZZ/101[x_0..x_3]
    for i to 10 list hilbertFunction(i,S)
    assert all (11,i -> rank P_(-i) === hilbertFunction(i,S))
///

TEST ///
    needsPackage "HyperplaneArrangements"
    A=typeA(2)
    I=orlikSolomon(A)
    E=ring I
    B=typeA(2)
    J=orlikSolomon(B,E)
    M=(ring I)^1/I ++ (ring J)^1/J
    P = injectiveResolution(M, LengthLimit => 7)
    assert isWellDefined P
    assert isFree P
    f = coaugmentationMap P
    assert isWellDefined f
    assert isQuasiIsomorphism f
    assert isComplexMorphism f

    A=typeA(2)
    I=orlikSolomon(A)
    E=ring I
    M=random(E^2,E^{-1,-2,-2})
    J=image M
    N=(ring I)^1/I ++ J
    P = injectiveResolution(N, LengthLimit => 7)
    assert isWellDefined P
    assert isFree P
    f = coaugmentationMap P
    assert isWellDefined f
    assert isQuasiIsomorphism f
    assert isComplexMorphism f
///

TEST ///
    --testing on complexes--
    E = ZZ/101[e_0..e_3,SkewCommutative => true]
    C= complex E
    Res = injectiveResolution(C, LengthLimit => 4)
    assert(Res == C)
    assert isWellDefined Res
    f = injectiveResolutionMap(C,LengthLimit=>5)
    assert isWellDefined f
    assert isQuasiIsomorphism f

    I = ideal (e_0+e_1,e_0*e_2)
    C = Hom(freeResolution(I, LengthLimit => 5), comodule I);
    P = injectiveResolution(C, LengthLimit => 5)
    assert isWellDefined P
    assert isFree P
    f = injectiveResolutionMap C -- bug  --new code fixed this bug
    assert isWellDefined f
    assert isQuasiIsomorphism f
    assert isComplexMorphism f
    source f==C
    target f==P --target f can be longer than P, so this can be false

    Mat=random(E^3,E^{-1,-2})
    prune ker Mat
    CMat=complex Mat
    P=injectiveResolution(CMat,LengthLimit=>5)
    f = injectiveResolutionMap CMat
    assert isQuasiIsomorphism f
    assert isComplexMorphism f
    
///

TEST ///
    E = QQ[e_0, e_1,e_2, SkewCommutative=>true]
    m = matrix{{e_0, e_0*e_1*e_2}}
    D = priddyDifferential(-2, m)
    C = priddyComplex(m, LengthLimit => 3)
    assert(D == map(E^{{6}, {8}, {10}, {12}}, E^{{5}, {7}, {9}},
	    {{e_0, 0, 0}, {e_0*e_1*e_2, e_0, 0}, {0, e_0*e_1*e_2, e_0}, {0, 0, e_0*e_1*e_2}}))
    assert isWellDefined C
    assert isHomogeneous C
    assert(C.dd_(-2) == D)
///

TEST ///
    -- Examples tried
    -- This next example doesn't make sense, because one of the forms is even degree.
    E = QQ[e_0, e_1, e_2, e_3, SkewCommutative=>true]
    --m = matrix{{e_0*e_1, e_2*e_3}}
    m = matrix{{e_0*e_1*e_2, e_1*e_2*e_3}}

    priddyDifferential(-2, m)
    priddyComplex(m, LengthLimit => 3)
///

TEST /// -- testing with S^1 and E^1
  (S, E) = koszulPair(n = 2, ZZ/101)

  -- testing koszulRR and koszulLL for Module and Complex
  M = S^1
  C = koszulRR(M, Concentration => (-n,0))
  assert(id_C === koszulRR(id_M, Concentration => (-n,0)))
  P = priddyComplex(vars E, LengthLimit => n)
  assert(C == P)
  F = koszulLL C
  assert(keys F.module == {0,1,2,3})
  assert(prune HH F == complex comodule truncate(n+1, S))

  N = E^{n+1}
  D = koszulLL N
  assert(keys D.module == {0,1,2,3})
  assert(id_D === koszulLL id_N)
  K = koszulComplex vars S
  -- TODO: figure out of it's possible to find the isomorphism
  assert(betti D == betti K)
  F0 = koszulRR(D, Concentration => (-n,0))
  -- Note: we take the canonical truncation to discard the edge homology
  F = canonicalTruncation(F0, (-n+1,0));
  assert(prune HH F == complex N)
///

TEST ///
  (S, E) = koszulPair(n = 2, ZZ/101)

  C = koszulLL E^1
  D = koszulRR(S^1, Concentration => (-n-1,0))
  assert(complex E == prune HH canonicalTruncation(koszulRR(C, Concentration => (-n,0)), -n+1, 0))
  assert(complex comodule truncate(n+2, S) == prune HH koszulLL D)

  -- LL(RR(finite length module)) is identity
  (S, E) = koszulPair(n = 2, ZZ/101)

  M = comodule truncate(n+2, S)
  assert(complex M == prune HH koszulLL koszulRR M)
  -- FIXME: complex M == koszulLL koszulRR M

  -- RR(LL(perfect complex)) is identity
  N = coker matrix {{e_0}}
  D = freeResolution(N, LengthLimit => n)
  C0 = koszulRR(koszulLL D, Concentration => (-n,n))
  C = canonicalTruncation(C0, (-n+1, n-1));
  assert(complex N == prune HH C)
  -- FIXME: complex N == C
///

TEST ///
  (S, E) = koszulPair(n = 2, ZZ/101)

  f = matrix {{x_0}}
  assert isWellDefined koszulRR(f, Concentration => (-n,0))
  assert isWellDefined koszulRR ker f
  assert isWellDefined koszulRR(coker f, Concentration => (-n,0))
  assert isWellDefined koszulRR(image f, Concentration => (-n,0))
  assert isWellDefined koszulRR(complex {f}, Concentration => (-n,0))

  assert isWellDefined koszulLL freeResolution(coker vars E, LengthLimit => 3)
  assert isWellDefined koszulRR(freeResolution(coker vars S, LengthLimit => 3), Concentration => (-5,0))
  assert isWellDefined koszulRR(koszulComplex vars S, Concentration => (-5,0))
  -- Note: koszulComplex vars E does not make sense

  assert isWellDefined(C = priddyComplex vars S)
  assert isWellDefined(koszulLL(C, Concentration => (-5,0)))
  assert isWellDefined(D = priddyComplex(vars E, LengthLimit => 2))
  assert isWellDefined(koszulRR D)
///

TEST /// -- testing with the Koszul complex
  (S, E) = koszulPair(n = 2, ZZ/101)

  C = koszulComplex vars S
  D = koszulRR(C, Concentration => (-n,0))
  C' = koszulLL(D, Concentration => (0,n))

  -- id_(RR(C)) == RR(id_C)
  assert(id_D == koszulRR(id_C, Concentration => (-n,0)))
  -- id_(LL(D)) == LL(id_D)
  assert(id_C' == koszulLL(id_D, Concentration => (0,n)))
///

TEST ///
    (S, E) = koszulPair(4, ZZ/19937)
    M = matrix { { S_0, S_1, S_2, S_3 }, {S_3, S_0, S_1, S_2}, {S_2, S_3, S_0, S_1}, {S_1, S_2, S_3, S_0} }
    f = map(S^{-1}^4, S^4, M)
    N = matrix { { S_0, S_1, S_2, S_3 }, {S_1, S_2, S_3, S_0}, {S_2, S_3, S_0, S_1}, {S_3, S_0, S_1, S_2} }
    g = map(S^4, S^{1}^4, N)
    assert(koszulRR(f, Concentration=>(-5,5)) * koszulRR(g, Concentration=>(-5,5)) == koszulRR(f * g, Concentration=>(-5,5)))
///

TEST ///
    (S, E) = koszulPair(4, ZZ/19937)
    M = matrix { { E_0, E_1, E_2, E_3 }, {E_3, E_0, E_1, E_2}, {E_2, E_3, E_0, E_1}, {E_1, E_2, E_3, E_0} }
    f = map(E^{-1}^4, E^4, M)
    N = matrix { { E_0, E_1, E_2, E_3 }, {E_1, E_2, E_3, E_0}, {E_2, E_3, E_0, E_1}, {E_3, E_0, E_1, E_2} }
    g = map(E^4, E^{1}^4, N)
    assert(koszulLL(f, Concentration=>(-5,5)) * koszulLL(g, Concentration=>(-5,5)) == koszulLL(f * g, Concentration=>(-5,5)))
///

TEST ///
    R = ZZ/19937[a,b,c,d]
    S = simplicialComplex {a*b*c, b*c*d, a*d}
    exteriorStanleyReisner(S, ZZ/19937)
    exteriorStanleyReisner(S, ZZ/101)
///

end--

--------------------------------------------------
--- Development
--------------------------------------------------

restart
needsPackage "ExteriorResolutions"
check ExteriorResolutions

(S,E)= koszulPair(1, ZZ/101)

koszulRR(koszulComplex matrix {{x_0}}, Concentration=>(-5,5))
koszulRR(koszulComplex vars S, Concentration=>(-5,5))

C = koszulComplex matrix {{x_0}}
koszulRR(id_C, Concentration=>(-5,5))

M = E^1

C = koszulComplex vars S

f = koszulRR(id_C, Concentration=>(-5,5))
f == id_(koszulRR(C, Concentration=>(-5,5)))


koszulRR(koszulLL(E^1, Concentration=>(-5,5)), Concentration=>(-5,5))


C = koszulComplex matrix{{x_0}}
koszulRR(C, Concentration=>(-5,5))

assert( f == id_(koszulLL(M, Concentration=>(-5,5))))

C = koszulLL(M, Concentration=>(0,5))

P = priddyComplex(vars E, LengthLimit=>3)


M = coker matrix{{x_0}}
E = ZZ/101[e_0,e_1,e_2, SkewCommutative=>true]

-- the following was pasted and modified from the BGG package
S := ring(M)
numvarsE := numgens E
ev := map(E,S,vars E)
f0 := basis(i,M)
f1 := basis(i+1,M)
-- we construct the differential by factoring it through (vars S)**f0
g := ((vars S)**f0)//f1 -- g is a map from source (vars S)**f0 to source f1 making the traingle involving (vars S)**f0 and f1 commute
b := (ev g)*((transpose vars E)**(ev source f0))
--correct the degrees (which are otherwise wrong in the transpose)
map(E^{(rank target b):i+1},E^{(rank source b):i}, b)

methods koszulComplex
code 0
