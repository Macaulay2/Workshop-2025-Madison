
getRandomElementsInModule = method()
getRandomElementsInModule(Module, ZZ) := Matrix => (M, n) -> (
    matrix {for i to n-1 list matrix random M}
)

makeRandomCpMackeyFunctor = method()
makeRandomCpMackeyFunctor(ZZ, List):= CpMackeyFunctor => (p, L) -> (

    if not isPrime p then error " -- p is not prime!";
    if not (length L === 4) then error " -- expected a list of length 4.";

    n := L#0;
    m := L#1;
    k := L#2;
    l := L#3;

    if not (class n === ZZ and class m === ZZ and class k === ZZ and class l === ZZ) then error " -- not a list of integers.";


    A := makeBurnsideMackeyFunctor p;
    B := makeUnderlyingFreeMackeyFunctor p;

    X := directSum ((for i to n-1 list A) | (for j to m-1 list B));

    XUnderlying := getUnderlyingModule X;
    XFixed := getFixedModule X;

    RandomUnderlyingElements := getRandomElementsInModule(XUnderlying, k);
    RandomFixedElements := getRandomElementsInModule(XFixed, l);



    return prune cokernel makeUniversalMap(X, RandomFixedElements, RandomUnderlyingElements)

)
