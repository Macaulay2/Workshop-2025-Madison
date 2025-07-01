newPackage("ExteriorResolutions",
    Version => "1.1",
    Date => "24 June 2025",
    Headline => "Injective resolutions over exterior algebras",
    Authors => {
	{Name => "Penelope Beall",          Email => "pbeall@ucdavis.edu",	        HomePage => "https://pbeall.github.io" },
	{Name => "Michael K. Brown",        Email => "mkb0096@auburn.edu",	        HomePage => "http://webhome.auburn.edu/~mkb0096/" },
	{Name => "Caitlin M. Davis",	    Email => "cmdavis22@wisc.edu",	        HomePage => "https://sites.google.com/wisc.edu/caitlindavis/home" },
	{Name => "Andrew Karsten",	        Email => "akk0071@auburn.edu",	        HomePage => "https://www.auburn.edu/cosam/departments/math/students/grad/graduate-students.htm" },
	{Name => "Jiaheng Li",		        Email => "henryli@gatech.edu",	        HomePage => "" },
	{Name => "Jianuo Zhou",		        Email => "jzhou632@gatech.edu",	        HomePage => "https://math.gatech.edu/people?field_job_type_tid=14"},
	{Name => "Boyana Martinova",	    Email => "boyana.martinova@gmail.com",  HomePage => "https://sites.google.com/view/bmartinova/home"},
	{Name => "Mahrud Sayrafi",	        Email => "mahrud@umn.edu",		        HomePage => "https://mahrud.github.io/" },
	{Name => "Gregory G. Smith",	    Email => "ggsmith@mast.queensu.ca",     HomePage => "https://mast.queensu.ca/~ggsmith/" },
	{Name => "Sreehari Suresh Babu",    Email => "sreeharisbabu183@gmail.com",  HomePage => "https://sreehari183.github.io/" }
    },
    PackageExports => {"Complexes"},
    Keywords => {"Commutative Algebra"}
  )

export {
    --methods
    "injectiveResolution",
    "coaugmentationMap",
    "priddyComplex",
    "priddyDifferential",
    }

--------------------------------------------------
--- Injective resolutions
--------------------------------------------------

--Input: 
--Output:
injectiveResolution = method(Options => options freeResolution);
injectiveResolution(Module) := Complex => opts -> (M) -> (
    E := ring M;
    n := numgens E;
    if not isSkewCommutative E then error "Expected underlying ring to be skew commutative";
    P := Hom(freeResolution(Hom(M,E), opts),E);
    P.cache.Module = M;
    P
)


coaugmentationMap = method()
coaugmentationMap Complex := ComplexMap => 
    (cacheValue symbol coaugmentationMap)(C -> (
            if not C.cache.?Module then error "expected an injective resolution";
            M := C.cache.Module;
            map(C, complex M, i -> if i === 0 then map(C_0, M, generators M // generators C_0))
            )
        )

priddyComplex = method(TypicalValue=>Complex, Options=>{LengthLimit=>null})
priddyComplex(Matrix, Ring) := opts -> (m, S) -> (
    complex hashTable apply(opts.LengthLimit, i -> -i=>priddyDifferential(-i, m, S))
)

priddyDifferential = method(TypicalValue=>Matrix)
priddyDifferential(ZZ, Matrix, Ring) := (i, m, S) -> (
    E := ring m;
    --L := flatten (degrees m)_1; --degrees of the forms on which we are taking the Priddy complex
    --assert(all(L, i -> odd L_i));
    monsSrc := basis(-i, S);
    monsTgt := basis(-i+1, S);

    expTgt := apply(first entries monsTgt, n->flatten exponents n);
    
    tgt := directSum apply(expTgt, n->
	E^{sum(numcols m, i->n_i*flatten degree(m_i))});
    

    f :=  matrix table(numcols monsTgt, numcols monsSrc, (r,c)->(
	    monsTgt_{r}//monsSrc_{c}
    ));

    map(tgt, , sub(f, m))
)


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
       Example
	   S = QQ[x_0,x_1]
	   E = QQ[e_0,e_1, SkewCommutative=>true]
	   m = matrix{{e_0, e_0*e_1}}
	   C = priddyComplex(m, S, LengthLimit=>3)
	   C.dd
	   prune HH C
///


--- TESTS
-* 
restart
needsPackage "ExteriorResolutions"
*-
TEST /// 
    E = ZZ/101[e_0..e_3,SkewCommutative => true]
    Res = injectiveResolution(E^1, LengthLimit => 4)
    assert(Res == complex E)
    assert isWellDefined Res
    f = coaugmentationMap Res
    assert isWellDefined f
    assert isQuasiIsomorphism f

    I = ideal vars E
    k = E^1/I
    P = injectiveResolution(k, LengthLimit => 10)
    P_0
    gens P_0
    prune ker dd^P_0 == k
    gens k
    assert isWellDefined P
    f = coaugmentationMap P
    assert isWellDefined f
    assert isQuasiIsomorphism f
    assert isComplexMorphism f
    assert isCommutative f
    degree f
    prune ker dd^P_0

    g = map(Hom(k,E),Hom(P_0,E),1)
    h = Hom(g,E^1)
    prune source h == k
    target h == P_0
    ambient source h
    ambient k
    Hom(Hom(k,E),E)
    
    dd^P_0*f_0
    f_0
    dd^P_0

    S = ZZ/101[x_0..x_3]
    for i to 10 list hilbertFunction(i,S)
    assert all (11,i -> rank P_(-i) === hilbertFunction(i,S))
    dd^P
    
    J = ideal(e_0*e_1, e_1*e_2)
///

TEST ///
    S = QQ[x_0,x_1]
    E = QQ[e_0, e_1,e_2, SkewCommutative=>true]
    m = matrix{{e_0, e_0*e_1*e_2}}
    D = priddyDifferential(-2, m, S)
    C = priddyComplex(m, S, LengthLimit=>3)
    assert(D == map(E^{{3}, {5}, {7}, {9}},E^{{2}, {4}, {6}},{{e_0, 0, 0}, {e_0*e_1*e_2, e_0, 0}, {0, e_0*e_1*e_2, e_0}, {0, 0, e_0*e_1*e_2}}))
    assert isWellDefined C
    assert isHomogeneous C
    assert(C.dd_(-2) == D)
///


TEST ///
    -- Examples tried
    -- This next example doesn't make sense, because one of the forms is even degree. 
    restart
    needsPackage "ExteriorResolutions"
    S = QQ[x_0,x_1]
    E = QQ[e_0, e_1, e_2, e_3, SkewCommutative=>true]
    --m = matrix{{e_0*e_1, e_2*e_3}}
    m = matrix{{e_0*e_1*e_2, e_1*e_2*e_3}}

    priddyDifferential(-2, m, S)
    priddyComplex(m, S, LengthLimit=>3)
///



end--

restart

needsPackage "ExteriorResolutions"
installPackage "ExteriorResolutions"
viewHelp "ExteriorResolutions"

S = QQ[x_0,x_1]
E = QQ[e_0,e_1,  SkewCommutative=>true]

m = matrix{{e_0^10, e_0*e_1}}

C = priddyComplex(m, S, LengthLimit=>3)
C.dd
prune HH C
