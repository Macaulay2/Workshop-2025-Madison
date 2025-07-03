--  CREDITS  --
--  Advisor: Francesca Gandini
--  Theory: Lucas Rizzolo (2021)
--  Initial Code: Lucas Rizzolo (2021)
--      * Last Editied - Gordon Novak 07/02/25
--  Code-Cleanup: Gordon Novak
--      * Last Editied - Gordon Novak 07/03/25
--  Documentation: Gordon Novak & Sasha Arasha
--      * Last Documented - 07/03/25

-- //////////////// --
-- //////////////// --
-- CODE STARTS HERE
-- //////////////// --
-- //////////////// --

needsPackage "InvariantRing"

-- METHOD_NAME: expandseeds
-- USAGE: Finds the minimal generating seed invariants for an invariant ring
--      INPUT: 
--          * L     : List              => contains generating seeds for expansion.
--          * ZList : List              => contains the list of moduluses of each variable.
--      OUTPUT:
--          * N     : List              => list of minimal generating invariants for the ring.
expandseeds = method();
expandseeds (List, List) := (L, ZList) ->(

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

            );  
        );

        if (candidate != 1_R) then (
            N = N | {M_i};
        )
    );

    -- Finally, we return M
    return N;
)

export {"expandseeds"}

beginDocumentation()

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
        "R = QQ[x_1..x_5]"
        "ZList = {3,3,3,3,3};",
        "expandseeds({v^2 * y, w^2 * y, x^2 * y}, ZList)"
    },
    EXAMPLE {
        "R = QQ[v,w,x,y,z];",
        "ZList = {4,4,4,4,4};",
        "W = matrix{{1,1,1,1,1},{1,1,1,1,0}};",
        "expandseeds(genseeds(R,W,ZList), ZList);"
    },


    SUBSECTION {"Ways to use ", TT "expandseeds", "."},
    UL {
        CODE {"expandseeds(R, W, ZList)"}
    }

}
