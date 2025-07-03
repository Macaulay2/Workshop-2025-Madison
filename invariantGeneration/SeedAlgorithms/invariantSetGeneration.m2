-- Created July 3rd, 2025

--  CREDITS  --
--  Advisor: Francesca Gandini
--  Theory: Marcus Cassel & Sumner Strom
--  Code: Gordon Novak
--      * Last Editied - Gordon Novak 07/03/25
--  Documentation: Gordon Novak
--      * Last Documented - 07/03/25
--  Review: Marcus Cassel & Sumner Strom
--      * Last Reviewed - 07/03/25

-- //////////////// --
-- //////////////// --
-- CODE STARTS HERE
-- //////////////// --
-- //////////////// --

needsPackage "InvariantRing"

-- METHOD_NAME: genseeds
-- USAGE: Finds the minimal generating seed invariants for an invariant ring
--      INPUT: 
--          * R     : polynomialRing    => Ring being acted upon
--          * W     : matrix            => Weight matrix representing the group action
--      OUTPUT:
--          * N     : List              => List of generating seeds
genseeds = method();
genseeds(PolynomialRing, Matrix, List) := (R, W, ZList) -> (

    -- Setting n & m for ease of use
    n = numColumns W;
    m = numRows W;

    -- The variable that holds the full rank submatrix
    subVars = matrix{{0}};
    -- We need to keep track of the columns of our full rank matrix
    colList = {};

    -- Keeps track of the indicies of the columns from which we'll form the submatricies
    mList = toList (0..(m-1));

    -- Now, we start going through sets of m x m submatricies to check if the determinant is nonzero.
    for i to n - m do (
        -- Form the submatrix of size m x m
        rod = submatrix(W, mList, toList (i..(i + m - 1)));

        -- If determinant is nonzero, we'll use it:
        if (determinant rod != 0) then (
            subVars = rod;
            colList = toList (i..(i + m - 1));
            break;
        );
    );
    if (subVars == matrix{{0}}) then (
        print ("error: i cant do it.\n");
        return {};
    );

    -- Keeps track of each of the vectors representing the exponents 
    vecExp = {};
    
    for v in toList(set(toList (0..(n-1))) - set(colList)) do (
        wColumns = sort({v} | colList); 
        -- Creates a Census of gerballs
        tempVec = {};
        -- Sets the signFlip to 1 to start
        signFlip = 1;
        -- This for loop allows us to do a cofactor expaction on the chosen sub matrix, and puts it into vector.
        for i in wColumns do (
            -- e is basis vector e_j (e.g. {0,1,0,...,0})
            e = for j from 0 to n-1 list (if j == i then 1 else 0);
            
            -- Makes a square sub matrix that is removing one of the columns
            subM = submatrix(W, mList, sort(toList(set(wColumns) - set({i}))));
            -- scales the unit vector
            tempVec = tempVec | {signFlip * (determinant subM ) * e};
            
            -- F lip the signFlip
            signFlip = signFlip * -1;
        );
        --  Takes every tempVec and adds them together to form a vector that represents exponents.  
        vecExp = vecExp | {fold((a,b) -> (for i to #a - 1 list (a_i + b_i)), tempVec)};
        
    );


    expR = {};
    for vec in vecExp do (
        seed = 1;
        for e to (#vec - 1) do (
            if (vec#e != 0) then (
                seed = seed * (gens R)#e^((ZList#e + vec#e) % ZList#e);
            );
        );
        expR = expR | {seed};
    );

    return expR;

)

-- METHOD_NAME: expandseeds
-- USAGE: Finds the minimal generating seed invariants for an invariant ring
--      INPUT: 
--          * L     : List              => contains generating seeds for expansion.
--          * ZList : List              => contains the list of moduluses of each variable.
--      OUTPUT:
--          * N     : List              => list of minimal generating invariants for the ring.
expandseeds = method();
expandseeds (List, List) := (L, ZList) ->(

    R = ring(L_0);
    gR = gens(ring L_0);
    -- M will be all of the elements in our group and N will include only minimal elements of M

    M = {};
    -- This code accumulates the m' over and over and then mods out by our modDegree
    for m in L do (
        M = M | {m};
        m' = m;
        for p in ZList when (m' != 0) do (

            
            -- Then, we multiply out a power. 
            m' = m' * m;

            for j to (#gR - 1) do (
                while (degree(gR#j, m') > ZList#j) do (
                    m' = lift(m' / ((gR#j)^(ZList#j)), R);
                );
            );

            -- First we add m' to our M
            M = M | {m'};

        );
    );

    -- This for loop goes through the variables and appends the monomials to M
    for i from 0 to (numgens R - 1) do (
        M = M | {(gR#i)^(ZList#i)};
    );

    --remove duplicate monomials
    M = unique(M);
    M = sort(M);


    N = {};

    -- Now, we're going to check if we can make 
    for i when (i < #M) do (

        candidate = M_i;
        for j when (j < #M) do (
            if (i != j and (candidate % M_j) == 0) then (
                
                candidate = candidate // M_j;

                if candidate == 1_R then (
                    break;
                );
            );  
        );

        if (candidate != 1_R) then (
            N = N | {M_i};
        )
    );

    -- Finally, we return M
    return N;
)

invargenset = method();
invargenset(PolynomialRing, Matrix, List) := (R, W, ZList) -> (
    return expandseeds(genseeds(R, W, ZList), ZList);
)

export {
    "genseeds",
    "expandseeds",
    "invargenset"
}

beginDocumentation()

document {
    Key => genseeds, 

    Headline => "Finds the minimal generating seed invariants for an invariant ring.",

    Usage => "genseeds(PolynomialRing, Matrix, List) := (R, W, ZList)",

    Inputs => {
        "R" => PolynomialRing => {"Ring upon which the group action is applied."},
        "W" => Matrix => {"The weight matrix representing the group action."},
        "ZList" => List => {"The list of the dimensions of the weight matrix."}
        },    

    Outputs => {
        "N" => List => {"A list of the generating seeds."},
    },

    PARA {"This function returns the generating seeds that can be used to find all invariants under a group action."},

    SUBSECTION "Examples",
    PARA {"We can use this algorithm in rings of any number of variables. For example, "},
    EXAMPLE {
        "R = QQ[x_0..x_2];",
        "W = matrix{{0,1,1},{1,0,1}};",
        "ZList = {2,2,2};",
        "genseeds(R,W,ZList)"
    }
}

document{
    Key => invargenset, 

    Headline => "Finds the generating set of invariants under a group action.",

    Usage => "invargenset(R, W, ZList)",

    Inputs => {
        "R" => PolynomialRing => {"Ring upon which the group action is applied."},
        "W" => Matrix => {"The weight matrix representing the group action."},
        "ZList" => List => {"The list of the dimensions of the weight matrix."}
        },    

    Outputs => {
        "N" => List => {"A list of the generating seeds."},
    },

    PARA {"This function returns the generating set of invariants under a group action."},

    SUBSECTION "Examples",
    PARA {"We can use this algorithm in rings of any number of variables. For example, "},
    EXAMPLE {
        "R = QQ[x_0..x_2];",
        "W = matrix{{0,1,1},{1,0,1}};",
        "ZList = {2,2,2};",
        "invargenset(R, W, ZList)"
    },
    EXAMPLE {
        "R = QQ[v,w,x,y,z];",
        "ZList = {4,4,4,4,4};",
        "W = matrix{{1,1,1,1,1},{1,1,1,1,0}};",
        "invargenset(R,W,ZList);"
    }
}

document {
    Key => expandseeds, 
    Headline => "Expands a complete set of seed generators to create the minimal generating set of an invariant ring.",
    Usage => "expandseeds(R,W,ZList)",
    Inputs => {
        "L" => List => {"seeds that can be used to generate the ring of invaraints."},
        "ZList" => List => {"dimensions of the weight matrix."}
    },    

    Outputs => {
        "Basis" => List => {"set of generators for the ring of invariants under a group action."},
    },

    PARA {"This function returns the set of generators for a ring of invariants."},

    SUBSECTION "Examples",
    PARA {"We can use this algorithm in rings of any number of variables, but the dimensions of the weight matrix utilized must be known. For automatic utilization, either use the ", TT "invargenset", " method, or combine with the ", TT "genseeds", " method."},
    EXAMPLE {
        "R = QQ[v,w,x,y,z]",
        "ZList = {3,3,3,3,3};",
        "expandseeds({v^2 * y, w^2 * y, x^2 * y}, ZList)"
    },


    SUBSECTION {"Ways to use ", TT "expandseeds", "."},
    UL {
        CODE {"expandseeds(R, W, ZList)"}
    }

}
