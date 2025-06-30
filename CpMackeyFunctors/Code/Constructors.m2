------------------------------------------
-- Constructors for common Mackey functors
------------------------------------------

BurnsideMackeyFunctor = method()
BurnsideMackeyFunctor(ZZ) := CpMackeyFunctor => (p) -> (
    C := matrix{{1}};
    T := matrix{{0},{1}};
    R := matrix{{1,p}};

    return makeCpMackeyFunctor(p,R,T,C);
)
