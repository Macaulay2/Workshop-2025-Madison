-- Implementing box product

-- Create the fixed module of the box product *before* modding out by the equivalence relation
genFixedModuleBoxProduct = method()
genFixedModuleBoxProduct (CpMackeyFunctor,CpMackeyFunctor) := Module => (M,N) ->(
    return (getFixedModule(M) ** getFixedModule(N))++(getUnderlyingModule(M) ** getUnderlyingModule(N))
)

syzFixedModuleBoxProduct = method()
syzFixedModuleBoxProduct (CpMackeyFunctor,CpMackeyFunctor) := Module => (M,N) ->(
    return ((getUnderlyingModule(M) ** getFixedModule(N)) ++ (getFixedModule(M) ** getUnderlyingModule(N)) ++ (getUnderlyingModule(M) ** getUnderlyingModule(N)))
)

boxProductRelationOne = method()
boxProductRelationOne (CpMackeyFunctor,CpMackeyFunctor) := Matrix => (M,N) ->(
    outputMapOne := getTransfer(M) ** id_(getFixedModule(N));
    outputMapTwo := id_(getUnderlyingModule(M)) ** getRestriction(N);
    return matrix({{outputMapOne},{-outputMapTwo}})

)

boxProductRelationTwo = method()
boxProductRelationTwo (CpMackeyFunctor,CpMackeyFunctor) := Matrix => (M,N) ->(
    outputMapOne := id_(getFixedModule(M)) ** getTransfer(N);
    outputMapTwo := getRestriction(M)** id_(getUnderlyingModule(N));
    return matrix({{outputMapOne},{-outputMapTwo}})

)

boxProductRelationThree = method()
boxProductRelationThree (CpMackeyFunctor,CpMackeyFunctor) := Matrix => (M,N) ->(

    targetOfZero:= getFixedModule(M) ** getFixedModule(N);
    sourceOfZero:=getUnderlyingModule(M) ** getUnderlyingModule(N);

    outputMapOne := map(targetOfZero, sourceOfZero,0);

    outputMapTwo := getConjugation(M)** getConjugation(N) - id_(getUnderlyingModule(M))** id_(getUnderlyingModule(N));
    return matrix({{outputMapOne},{-outputMapTwo}})
)

totalBoxProductRelation = method()
totalBoxProductRelation (CpMackeyFunctor,CpMackeyFunctor) := Matrix => (M,N) ->(
    return boxProductRelationOne(M,N) | boxProductRelationTwo(M,N) | boxProductRelationThree(M,N)
)

boxProductFixedModule = method()
boxProductFixedModule (CpMackeyFunctor,CpMackeyFunctor) := Matrix => (M,N) ->(
    coker(totalBoxProductRelation(M,N))
)

boxProductUnderlyingModule = method()
boxProductUnderlyingModule (CpMackeyFunctor,CpMackeyFunctor) := Matrix => (M,N) ->(
    getUnderlyingModule(M)**getUnderlyingModule(N)
)

-- Define the transfer map for box products
boxProductTransfer = method()
boxProductTransfer (CpMackeyFunctor,CpMackeyFunctor) := Matrix => (M,N) ->(
    -- Look at the transfer valued in the generators module before modding out by relations
    mapBeforeInducing := map(getFixedModule(M) ** getFixedModule(N), getUnderlyingModule(M) ** getUnderlyingModule(N),0) || id_(getUnderlyingModule(M))** id_(getUnderlyingModule(N));

    return inducedMap(boxProductFixedModule(M,N),boxProductUnderlyingModule(M,N), mapBeforeInducing)
)



boxProductRestriction = method()
boxProductRestriction (CpMackeyFunctor,CpMackeyFunctor) := Matrix => (M,N) ->(

    mapOne:= getRestriction(M) ** getRestriction(N);

    p := getPrimeOrder(M);

    mapTwo:= sum for j from 0 to p-1 list(
        ((getConjugation(M))^j) ** ((getConjugation(N))^j)
    );

    return inducedMap(boxProductUnderlyingModule(M,N),boxProductFixedModule(M,N),mapOne | mapTwo)
)

boxProductConjugation = method()
boxProductConjugation (CpMackeyFunctor,CpMackeyFunctor) := Matrix => (M,N) ->(
    getConjugation(M) ** getConjugation(N)
)


boxProduct = method()
boxProduct (CpMackeyFunctor,CpMackeyFunctor) := CpMackeyFunctor => (M,N) ->(
    if getPrimeOrder(M) != getPrimeOrder(N) then error("-- Primes not the same (incompatable)");
    return makeCpMackeyFunctor(getPrimeOrder(M),boxProductRestriction(M,N), boxProductTransfer(M,N),boxProductConjugation(M,N))
)

CpMackeyFunctor ** CpMackeyFunctor := boxProduct