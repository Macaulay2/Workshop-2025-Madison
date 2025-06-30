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

FixedFreeMackeyFunctor = BurnsideMackeyFunctor

UnderlyingFreeMackeyFunctor = method()
UnderlyingFreeMackeyFunctor(ZZ) := CpMackeyFunctor => (p) -> (
    I := mutableMatrix id_(ZZ^p);                                    --declare identity matrix
    C := matrix(rowPermute(I,0,({p-1}|(toList (0..p-2)))));
    T := matrix({for i to p-1 list 1});
    R := matrix(for i to p-1 list {1});

    return makeCpMackeyFunctor(p,R,T,C);
)
