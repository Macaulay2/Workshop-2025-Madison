-- other attempt: Generate special Monomials of degree at most d in a constructive way
needsPackage "SRdeformations"

ListSpInd=(n,d)->(
    SpI:={{0}}; --Initializing the list of special powers
    SpDeg:={0}; --Initializing the list of special degrees
    AuxI:={};
    AuxDeg:={};
    for m from 1 to n-1 do(
        for i from 0 to (#SpI)-1 do(
            if (SpDeg_i)+(last (SpI_i))<=d then(AuxI=AuxI | { (SpI_i) | {last (SpI_i)}};
                                                AuxDeg=AuxDeg | { SpDeg_i+ (last (SpI_i))};
                                                if (SpDeg_i)+(last (SpI_i))+1<=d then( AuxI=AuxI | { (SpI_i) | {1 +(last (SpI_i))}};
                                                                                         AuxDeg=AuxDeg | { 1+ SpDeg_i+ (last (SpI_i))}))
                                             else(AuxI=AuxI | { sort ((SpI_i) | {0}) };
                                                  AuxDeg=AuxDeg | { SpDeg_i}  ));
        SpI=AuxI;
        SpDeg=AuxDeg;
        AuxI={};
        AuxDeg:={}
        ); 
    m:= length (SpI_0);
    C:={};
    if m < n then(for i from 1 to (n - m) do( C = C | {0});
                  for i from 0 to (#SpI - 1) do( SpI_i = SpI_i | C));
    SpI   
)
--Test
ListSpInd(2,6)
ListSpMon(2,6)
ListSpInd(6,4)
ListSpMon(4,1)
ListSpMon(4,6)
ListSpMon(4,8)
ListSpMon(7,8)
ListSpMon(8,27)


-- Shuffle every monomial
ShuffMon=(f,n)->(
    R:=QQ[x_1..x_n];
    P:=(exponents f)_0;
    P= permutations P;
    Mon:={};
    for i from 0 to (#P-1) do(
       Mon = Mon | {vectorToMonomial(vector(P_i), R) } 
    );
    toList set Mon
)

ListSpMon=(n,d)->(
    d=min {d,n*(n-1)//2};
    R:=QQ[x_1..x_n];
    SpI:=ListSpInd(n,d);
    SpMon:={};
    for i from 0 to (#SpI - 1) do(SpMon = SpMon| ShuffMon(vectorToMonomial( vector (SpI_i) , R  ) ,n));
    toList set SpMon
)

--Test
ListSpInd(2,6)
ListSpMon(2,6)

--Orbit Sum for one monomial
orbSum = (f,G,n) ->(
    R:=QQ[x_1..x_n];
    I:= (exponents f)_0;
    v:= transpose matrix{I};
    g:=0;
    S=G*v;
    for i from 0 to (#S-1) do(
        g = g + vectorToMonomial(vector(S_i),R);
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
