newPackage("OrbitSum",
    Version => "0.1",
    Date => "July 1st, 2025",
    Authors => {
        {
            Name => "Lilia Alanís-López", 
            Email => "lilia.alanislpz@tec.mx", 
            HomePage => "https://sites.google.com/tec.mx/liliaalanislopez/"
        },
	    {
            Name => "Dalena Vien",
            Email => "",
	        HomePage => ""
        },
        {
            Name => "Indalecio Ruiz",
            Email => "",
            HomePage => ""
        },
        {
            Name => "Francesca Gandini",
            Email => "",
            HomePage => ""
        },
        {   Name => "Summer Strom",
            Email => "",
            HomePage => ""
        } 
    },
    Headline => "Package for permutation groups and orbit sums",
    Keywords => {"Permutation groups","orbit sums","documentation"},
    
    DebuggingMode => true,
    Reload => true
    )

--Description => "This package computes orbit sums of special monomials under the action of a permutation group. It includes functions to generate check special monomials, and compute orbit sums for special monomials. The package is designed for use in invariant theory and related fields."



--We need a previous package
needsPackage "SRdeformations"

-- Getting the special Monomials (modulo permutations of the variables)


listSpInd=(n,d)->( --intakes the number of variables n and the degree d.
    SpI:={{0}}; --Initializing the list of special powers
    SpDeg:={0}; --Initializing the list of special degrees
    AuxI:={};
    AuxDeg:={};
    for m from 1 to n-1 do(
        for i from 0 to (#SpI)-1 do(
            if (SpDeg_i)+(last (SpI_i))<=d then(AuxI=AuxI | { (SpI_i) | {last (SpI_i)}}; -- Checking the degree of the current and last element degrees then appending them if they are leq d.
            AuxDeg=AuxDeg | { SpDeg_i+ (last (SpI_i))}; -- Appending the degree of the new element.
            if (SpDeg_i)+(last (SpI_i))+1<=d then( AuxI=AuxI | { (SpI_i) | {1 +(last (SpI_i))}}; -- Checking if this is now special, then if it is not repeating the process.
                                                        AuxDeg=AuxDeg | { 1+ SpDeg_i+ (last (SpI_i))}))
            else(AuxI=AuxI | { sort ((SpI_i) | {0}) }; -- Extending the index.
                AuxDeg=AuxDeg | { SpDeg_i}  )
        );
        SpI=AuxI;
        SpDeg=AuxDeg;
        AuxI={};
        AuxDeg:={} -- Updating the lists.
        ); 
    SpI  -- Returns the list of special indexes 
)

-- Shuffle every monomial
shuffMon=(f,n)->( -- f is a monomial and n is the number of variables.
    R:=QQ[x_1..x_n];
    P:=(exponents f)_0; -- Takes the exponents of the monomial and makes a list.
    P= permutations P; -- Gives all permutations of the exponent list.
    Mon:={};
    for i from 0 to (#P-1) do(
       Mon = Mon | {vectorToMonomial(vector(P_i), R)} -- Converts the permuted list of exponents back to a list of monomials.
    );
    toList set Mon -- Returns the list of permuted monomials.
)

--Generate, if posible, the list of special monomials
listSpMon=(n,d)->(
    d=min {d,n*(n-1)//2};
    R:=QQ[x_1..x_n];
    SpI:=listSpInd(n,d);
    SpMon:={};
    for i from 0 to (#SpI - 1) do(SpMon = (SpMon| shuffMon(vectorToMonomial( vector (SpI_i) , R  ) ,n)));  -- Takes the special index and converts it to a monomial in R of n variables.
    toList set SpMon --returns the list of special monomials
)

--Test for Special Indexes and Monomials
listSpInd(6,4)
listSpMon(2,6)
listSpMon(4,1)
listSpMon(4,6)



--Orbit Sum for one monomial
orbSum = (f,G,n) ->( -- Intakes a monomial f, a group action G, and the number of variables n.
    R:=QQ[x_1..x_n];
    I:= (exponents f)_0; -- Gets a list of the exponents of f.
    v:= transpose matrix{I}; -- Takes a transpose of matrix made by the list of exponents.
    g:=0;
    S := G*v; -- This takes the group action over the exponents.
    for i from 0 to (#(S)-1) do(  
        g = g + vectorToMonomial(vector(S)_i,R); --This extracts the exponent permuted by the group action and assigns it to a variable.
    );
    g -- Returns the orbit sum of the monomial f.
)
--Orbit Sums for special Monomials
orbSumList=(G,n,d)->( -- Intakes a group action G, the number of variables n, and the degree d.
    M:=listSpMon(n,d); -- This generates the list of special monomials.
    L:={};
        for i from 0 to (#M-1) do(
         L = L | { orbSum(M_i,G,n)}; -- This takes all the special monomials and takes the orbit sum of each.
        );
        toList(set L) -- returns a list of the orbit sums of the special monomials.
) 


export {"listSpInd",
        "listSpMon",
        "shuffMon",
	    "orbSum",
        "orbSumList"
}

beginDocumentation()

load "./Orbitsumdoc.m2"

end
