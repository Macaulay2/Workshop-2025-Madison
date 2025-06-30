--We need a previous package 
needsPackage "SRdeformations"

-- Getting the special Monomials (modulo permutations of the variables)

-- First we get the special indexes
IsListSp = (I,n) -> (
    R:=QQ[x_1..x_n];
    aux := 0;
    if (length I) > n then aux = 1 else(
        if (length I) < n then(
            for i from 1 to (n - (length I)) do I = I | {0}
        );
        for i from 1 to (#I - 1) do(
            if abs(I_i - I_(i-1))>1 then(aux=1; break)
        );
        if (last I) > 0 then aux = 1
    );
    if aux == 0 then {true , I} else {false , I}
)

--Generate, if posible, the list of special monomials
ListSpMon=(n,d)->(
    if d>(.5*n*(n-1)) then return "There are no special monomials with that degree";
    R:=QQ[x_1..x_n];
   L:=partitions d;
   M:={};
   for i from 0 to (#L-1) do if (IsListSp(toList L_i,n))_0 then M=M | ShuffMon(vectorToMonomial( vector((IsListSp(toList L_i,n))_1)  , R  ) ,n);
   if #M>0 then toList(set M) else print("There are no special monomials")
)
--Test
ListSpMon(4,1)
ListSpMon(4,8)
ListSpMon(7,8)
ListSpMon(27,8)

-- Shuffle every monomial
ShuffMon=(f,n)->(
    R:=QQ[x_1..x_n];
    P:=(exponents f)_0;
    P= permutations P;
    Mon:={};
    for i from 0 to (#P-1) do(
       Mon = Mon | {vectorToMonomial(vector(P_i), R) } 
    );
    Mon
)
ShuffMon(x_1^2*x_2^3*x_3^1*x_4^2,4)


--Orbit Sum for one monomial
orbSum = (f,G,n) ->(
    R:=QQ[x_1..x_n];
    I:= (exponents f)_0;
    v:= transpose matrix{I};
    g:=0;
    G*v;
    for i from 0 to (#(G*v)-1) do(
        g = g + vectorToMonomial(vector((G*v)_i),R);
    );
    g
)

--Orbit Sums for special Monomials
orbSumList=(G,n,d)->(
    M:=ListSpMon(n,d);
    L:={};
        for i from 0 to (#M-1) do(
         L = L | { orbSum(M_i,G,n)};
        );
        toList(set L)
) 



--Toy Examples 1
m = matrix{{0,0,0,1},{1,0,0,0},{0,1,0,0},{0,0,1,0}}
G = {m , m^2, m^3, m^4}
m2 = matrix{{0,1,0,0},{1,0,0,0},{0,0,0,1},{0,0,1,0}}
G2 = {m2, m2^2}
R = QQ[x_1..x_4]
f=x_1^2*x_2^3*x_3^1*x_4^2
orbSum(f,G,4)
orbSum(f,G2,4)
orbSumList(G,4,3)

--Toy Examples 1
-- Orbit Sum for {s_i}_{i\in\{1,\ldost n\}} where s_i=x_1*x_2*\ldots*x_i
    -- defining s_i
    SPoly=(i,n)->(
    R:=QQ[x_1..x_n];
    f:=1;
    for j from 1 to i do(f=f*x_j);
    f
    )
    SPoly(3,4)

    -- Getting the OrbitSum List for s_i polynomials in n variables
    OrbSumSi=(n,G)->(
        R:=QQ[x_1..x_n];
        L:={};
        for i from 1 to n do( L = L | { orbSum((SPoly(i,n)),G,n)});
        toList(set L)
    )
