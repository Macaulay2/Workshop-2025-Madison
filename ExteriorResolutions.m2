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
    PackageExports => {"Complexes"},
    AuxiliaryFiles => true,
    DebuggingMode  => true
    )

export {
    --methods
    "injectiveResolution",
    "coaugmentationMap",
    "priddyComplex",
    "priddyDifferential",
    "koszulRR",
    "koszulLL",
    "koszulPair",
    "koszulDual",
    "StanleyReisnerExterior",
}

--------------------------------------------------
--- Injective resolutions
--------------------------------------------------

injectiveResolution = method(Options => options freeResolution)
injectiveResolution Module := Complex => opts -> M -> (
    E := ring M;
    n := numgens E;
    if not isSkewCommutative E then error "expected underlying ring to be skew commutative";
    P := Hom(freeResolution(Hom(M, E), opts), E);
    P.cache.Module = M;
    P)
    
injectiveResolution Complex := Complex => opts -> C -> (
    E:=ring C;
    if not isSkewCommutative E then error "expected underlying ring to skew commutative";
    tempHom := Hom(C,E);
    fC := resolutionMap(tempHom,opts);
    P:= Hom(freeResolution(tempHom, opts), E);
    hfC := Hom(fC,E);
    P.cache.resolutionMap = fC;
    P
    )

coaugmentationMap = method()
coaugmentationMap Complex := ComplexMap => C -> C.cache.coaugmentationMap ??= (
    if C.cache.?Module then (
        M := C.cache.Module;
        map(C, complex M, i -> if i === 0 then map(C_0, M, transpose syz transpose presentation M))
    )
    else if C.cache.?resolutionMap then (
        hfC := C.cache.resolutionMap;
        D := source hfC;
        --tempMap := map(D, C, i -> map(D_i, C_i, transpose syz transpose presentation C_i));
        hfC
    )
    else error "Expected an injective resolution"
)

--------------------------------------------------
--- Priddy complex
--------------------------------------------------

priddyDifferential = method(TypicalValue => Matrix)
priddyDifferential(ZZ, Matrix, Ring) := (i, m, S) -> (
    E := ring m;
    n := numgens E - 1;
    --L := flatten (degrees m)_1; --degrees of the forms on which we are taking the Priddy complex
    --assert(all(L, i -> odd L_i));
    monsSrc := basis(-i, S);
    monsTgt := basis(-i+1, S);
    --
    expTgt := apply(first entries monsTgt, n -> flatten exponents n);
    --
    tgt := directSum apply(expTgt, d ->
	E^{{n+1}+sum(numcols m, i -> d_i * flatten degree(m_i))});
    --
    f := matrix table(numcols monsTgt, numcols monsSrc,
	(r,c) -> monsTgt_{r} // monsSrc_{c});
    -- TODO: can we avoid sub?
    map(tgt, , sub(f, m)))

priddyComplex = method(TypicalValue => Complex, Options => { LengthLimit => null })
priddyComplex(Matrix, Ring) := opts -> (m, S) -> (
    complex hashTable apply(opts.LengthLimit, i -> -i => priddyDifferential(-i, m, S)))

--------------------------------------------------
--- Koszul duality Functors
--------------------------------------------------

-- returns a pair S = Sym^* K^(n+1) and E = Wedge^* K^(n+1)
koszulPair = method()
koszulPair(ZZ, Ring) := (n, K) -> (
    x := getSymbol "x";
    e := getSymbol "e";
    S := K[x_0..x_n];
    E := K[e_0..e_n, SkewCommutative => true];
    S.cache.koszulDual = E;
    E.cache.koszulDual = S;
    (S, E))

koszulDual = method()
koszulDual Ring := A -> A.cache.koszulDual ??= (
    K := coefficientRing A;
    x := getSymbol "x";
    e := getSymbol "e";
    if isSkewCommutative A
    then K[x_0..x_(numgens A - 1)]
    else K[e_0..e_(numgens A - 1), SkewCommutative => true])

-- TODO: doesn't work over E yet
degreeSupport = M -> (
    if hilbertPolynomial M != 0 then error "expected a range provided as Concentration => {lo, hi}";
    (min degrees source generators M, max degrees source relations M))

-- RR: Com(S) -> Com(E) is the right-adjoint functor
koszulRR = method(Options => { Concentration => null })
-- RR(M)^i = E^*(i) \otimes_k M_i
-- E^* = Hom_k(E, k) = E(n+1)
koszulRR Module := Complex => opts -> M -> M.cache#(koszulRR, opts) ??= (
    S := ring M;
    n := numgens S - 1;
    E := koszulDual S;
    ev := map(E,S,vars E); -- not a map of algebras, just substiuting variables
    (lo, hi) := try opts.Concentration else degreeSupport M;
    modules := hashTable apply(lo..hi, i -> i => E^{n+1+i} ** (E ** part(i, M)));
    if lo == hi
    then complex(modules#lo, Base => -lo)
    else complex hashTable apply(lo..hi-1,
	i -> -i => (
	    src := basis(i, M);
	    tar := basis(i+1, M);
	    
	    -- we construct the differential by factoring it through (vars S)**src
	    g := ((vars S)**src)//tar; -- g is a map from source (vars S)**src to source tar making the traingle involving (vars S)**src and tar commute
	    b := (ev g)*((transpose vars E)**(ev source src));
	    -- the above 2 lines  were pasted and modified from the BGG package
	    
	    map(modules#(i+1), modules#i, b)
	    )))
-- RR(y**s) = \sum_{l=0}^n y*e_l ** s*x_l
koszulRR Matrix := ComplexMap => opts -> f -> f.cache#(koszulRR, opts) ??= (
    src := koszulRR(source f, opts);
    tar := koszulRR(target f, opts);
    E := koszulDual ring f;
    map(tar, src, i -> map(tar_i, src_i, E**part(-i, f)))
    )

-- RR(C)^i = \bigoplus_{j\in\ZZ} Hom_k(E(-j), C^{i-j}_j)
--         = \bigoplus_{j\in\ZZ} (E(-j))^* \otimes_k C^{i-j}_j
--         = \bigoplus_{j\in\ZZ} E^*(j)    \otimes_k C^{i-j}_j
koszulRR Complex    := Complex    => opts -> C -> (
    (lo,hi) := opts.Concentration; -- bounds for cohomological index i of RR(C)
    (inf,sup) := concentration C; -- bounds for k (homological index)

    RRterms := hashTable apply(inf..sup, k -> -k=>koszulRR(C_k, Concentration=>(lo+k,hi+k)));
    RRdiffs := hashTable apply((inf+1)..sup, k -> -k=>koszulRR(C.dd_k, Concentration=>(lo+k,hi+k)));
    
    modules := hashTable apply(lo..hi, i ->
    	i => hashTable apply(inf..sup, k -> -k=>(RRterms#(-k))^(i+k))
    	);

    complex hashTable apply(lo..(hi-1), i->
	i => matrix table(toList(inf..sup), toList(inf..sup),
	    (r,c)->map(modules#(-i-1)#(-r), modules#(-i)#(-c),
	    if r==c then (-1)^r * dd^(RRterms#(-r))_(-r+i)
	    else if r==c+1 then (RRdiffs#(-r))_(-r+i)
	    else 0
	    )))
    )
-- ???
koszulRR ComplexMap := ComplexMap => opts -> psi -> ()


-- LL: Com(E) -> Com(S) is the left-adjoint functor
koszulLL = method(Options => options koszulRR)
-- LL(N)^i = S(i) \otimes_k N_i
koszulLL Module := Complex  => opts -> N -> N.cache#(koszulLL, opts) ??= (
    E := ring N;
    n := numgens E - 1;
    S := koszulDual E;
    ev := map(S,E,vars S); -- not a map of algebras, just substiuting variables
    (lo, hi) := try opts.Concentration else degreeSupport N;
    modules := hashTable apply(lo..hi, i -> i => S^{i} ** (S ** part(i, N)));
    if lo == hi
    then complex(modules#lo, Base => -lo)
    else complex hashTable apply(lo..hi-1,
	i -> -i => (
	    src := basis(i, N);
	    tar := basis(i+1, N);
	    
	    -- we construct the differential by factoring it through (vars S)**src
	    g := ((vars E)**src)//tar; -- g is a map from source (vars S)**src to source tar making the traingle involving (vars S)**src and tar commute
	    b := (ev g)*((transpose vars S)**(ev source src));
	    -- the above 2 lines  were pasted and modified from the BGG package
	    
	    map(modules#(i+1), modules#i, (-1)^i*b)
	    )))
-- LL(s**y) = (-1)^i \sum_{l=0}^n x_l*s \otimes y*e_l
koszulLL Matrix := ComplexMap => opts -> f -> f.cache#(koszulLL, opts) ??= (
    src := koszulLL(source f, opts);
    tar := koszulLL(target f, opts);
    S := koszulDual ring f;
    map(tar, src, i -> map(tar_i, src_i, S**part(-i, f)))
)

-- LL(D)^i = \bigoplus_{j\in\ZZ} S(j) \otimes_k D^{i-j}_j
koszulLL Complex    := Complex    => opts -> D   -> ()
-- LL(s**y) = (-1)^i LL(s**y) + s ** dd_D(y) for y \in D^{i-j}_j
koszulLL ComplexMap := ComplexMap => opts -> phi -> ()


--------------------------------------------------
--- Stanley Reisner
--------------------------------------------------

needsPackage "SimplicialComplexes"

StanleyReisnerExterior = method()
StanleyReisnerExterior(SimplicialComplex, Ring) := (S, K) -> (

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
      (priddyComplex, Matrix, Ring)
   Headline
       computes the Priddy complex of several elements of an exterior algebra
   Usage
       priddyComplex(m, S)
   Inputs
       m:Matrix
	   consisting of a single row and with entries in an exterior algebra
       S:Ring
   Outputs
       :Complex
   Description
       Text
           Let $S = k[x_0, \dots, x_n]$ be a polynomial ring over a field. Let $f_1, \dots, f_c$ be odd degree elements of an exterior algebra $E = \Lambda_k(e_0, \dots, e_n)$. Then, we get a complex $E\otimes_k S_0\rightarrow E\otimes_k S_1\rightarrow \cdots$ whose differentials are left multiplication by $\sum_{i=0}^n f_i\otimes x_i$.

           Below we obtain the Priddy complex where $n=2$ and $f_1=e_0$ and $f_2=e_0*e_1$ and $k=\mathbb{Q}$.
       Example
	   S = QQ[x_0,x_1]
	   E = QQ[e_0,e_1, SkewCommutative=>true]
	   m = matrix{{e_0, e_0*e_1}}
	   C = priddyComplex(m, S, LengthLimit=>3)
	   C.dd
	   prune HH C
Node
   Key
     priddyDifferential
    (priddyDifferential, ZZ, Matrix, Ring)

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
///

--------------------------------------------------
--- Tests
--------------------------------------------------

-*
restart
needsPackage "ExteriorResolutions"
*-
TEST ///
    E = ZZ/101[e_0..e_3,SkewCommutative => true]
    Res = injectiveResolution(E^2, LengthLimit => 4)
    assert(Res == complex E^2)
    assert isWellDefined Res
    f = coaugmentationMap Res
    assert isWellDefined f
    assert isQuasiIsomorphism f

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

    --testing on complexes--
    E = ZZ/101[e_0..e_3,SkewCommutative => true]
    C= complex E
    Res = injectiveResolution(C, LengthLimit => 4)
    assert(Res == C)
    assert isWellDefined Res
    f = coaugmentationMap Res
    assert isWellDefined f
    assert isQuasiIsomorphism f

    I=ideal (e_0+e_1,e_0*e_2)
    C=Hom(freeResolution(I,LengthLimit=>5),comodule I)
    prune HH C
    P=injectiveResolution(C,LengthLimit=>5)
    assert isWellDefined P
    assert isFree P
    f = coaugmentationMap P
    assert isWellDefined f
    assert isQuasiIsomorphism f
    assert isComplexMorphism f
    
    Mat=random(E^3,E^{-1,-2})
    prune ker Mat
    CMat=complex Mat
    P=injectiveResolution(CMat,LengthLimit=>5)
    f=coaugmentationMap P

///

TEST ///
    S = QQ[x_0,x_1]
    E = QQ[e_0, e_1,e_2, SkewCommutative=>true]
    m = matrix{{e_0, e_0*e_1*e_2}}
    D = priddyDifferential(-2, m, S)
    C = priddyComplex(m, S, LengthLimit=>3)
    assert(D == map(E^{{6}, {8}, {10}, {12}},E^{{5}, {7}, {9}},{{e_0, 0, 0}, {e_0*e_1*e_2, e_0, 0}, {0, e_0*e_1*e_2, e_0}, {0, 0, e_0*e_1*e_2}}))
    assert isWellDefined C
    assert isHomogeneous C
    assert(C.dd_(-2) == D)
///


TEST ///
    -- Examples tried
    -- This next example doesn't make sense, because one of the forms is even degree.
    needsPackage "ExteriorResolutions"
    S = QQ[x_0,x_1]
    E = QQ[e_0, e_1, e_2, e_3, SkewCommutative=>true]
    --m = matrix{{e_0*e_1, e_2*e_3}}
    m = matrix{{e_0*e_1*e_2, e_1*e_2*e_3}}

    priddyDifferential(-2, m, S)
    priddyComplex(m, S, LengthLimit=>3)
///

TEST ///
    (S,E)= koszulPair(2, ZZ/101)
    M = S^1
    
    C = koszulRR(M, Concentration=>(0,3))
    P = priddyComplex(vars E, S, LengthLimit=>3)
    
    assert(C == P)

    F = koszulLL(HH_0 C, Concentration=>(-5,5))
    assert( F_0 == M)
///

TEST ///
    (S,E) = koszulPair(2, ZZ/101)
    M = E^{3}
    
    C = koszulLL(M, Concentration=>(-3,0))
    D = koszulComplex vars S

    assert(betti C == betti D)

    P = priddyComplex(vars S, E, LengthLimit=>3)

    F = koszulRR(HH_0 C, Concentration=>(-5,5))
    assert( F_0 == M)
///
TEST ///
    (S,E)= koszulPair(2, ZZ/101)

    M = S^1

    f = koszulRR(id_M, Concentration=>(-5,5))
    g = koszulRR(map(S^{1}, S^1, S_0), Concentration=>(-5,5))
    
    assert( f == id_(koszulRR(M, Concentration=>(-5,5))))
///
TEST ///
    (S,E)= koszulPair(2, ZZ/101)

    M = E^1

    f = koszulLL(id_M, Concentration=>(-5,5))
    
    assert( f == id_(koszulLL(M, Concentration=>(-5,5))))
///

end--

--------------------------------------------------
--- Development
--------------------------------------------------

restart
needsPackage "ExteriorResolutions"

(S,E)= koszulPair(1, ZZ/101)

M = E^1

f = koszulLL(id_M, Concentration=>(-5,5))

koszulRR(koszulComplex vars S, Concentration=>(-5,5))

C = koszulComplex matrix{{x_0}}
koszulRR(C, Concentration=>(-5,5))
>>>>>>> bf9227e (implemented koszulRR Complex, but did not finish debugging)

assert( f == id_(koszulLL(M, Concentration=>(-5,5))))

C = koszulLL(M, Concentration=>(0,5))

P = priddyComplex(vars E, S, LengthLimit=>3)

)
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
