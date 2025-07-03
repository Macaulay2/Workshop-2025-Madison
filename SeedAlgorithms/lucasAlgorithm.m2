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

-- METHOD_NAME: genseeds
-- USAGE: Finds the minimal generating seed invariants for an invariant ring
--      INPUT: 
--          * R     : polynomialRing    => Ring being acted upon
--          * W     : matrix            => Weight matrix representing the group action
--          * ZList : List              => List of the dimensions of the weight matrix
--      OUTPUT:
--          * N     : List              => List of minimal generating seeds
genseeds = method();
genseeds (PolynomialRing,Matrix,List) := (R,W,ZList) ->(
    
    -- First, we get the weight matrix via the diagonalAction function.
    T = diagonalAction(W,ZList,R);

    -- First we find the invariants under our diagonalAction (a list)
    S1 = invariants T;
    -- Then we sort the list via our lex ordering.
    S1 = sort S1;

    -- First make a shortcut variable for our generators under R
    gR = gens R;
    -- Then, we apply the variables of our polynomialRing to a subring.
    modVars = apply(gR, m -> sub(m, ring T));

    -- Now, we need to find the first element of S1 that isn't a pure power:
    -- First we make a guiding list of all the possible combinations of modded variables.
    comb = subsets(toList(0..(#modVars - 1)), 2);
    m := 0;
    -- Then, we check each possible two-element combination:
    -- Eg.. {0,1}, {0,2}, {1,2}, ect.
    for i in comb do (
        -- Afterwords, we run this combination for each s in our invariants
        for s in S1 do (
            -- If we get a nonzero combination, we've found an element of S1 that isn't a pure power, and we can end the loop.
            if (degree(modVars#(i_0), s) * degree(modVars#(i_1), s) != 0) then (
                m = s;
                break;
            )
        );
        -- This ends the loop.
        if (m != 0) then break;
    );

    -- m' will be our "accumulator" that we constantly multiply by m
    m' := m;

    -- sets the degree that each exponent will be modded by.
    -- this seems to work for Z/p x Z/p but I'm not sure about other groups.
    modDegree := lcm(ZList#0,ZList#1);

    -- M will be all of the elements in our group and N will include only minimal elements of M
    M := {};
    N := {};

    -- This for loop goes through the variables and appends the monomials to M
    for i from 0 to (numgens R - 1) do (
        M = M | {(modVars#i)^(modDegree)};
    );

    for i to (modDegree - 1) when (m' != 0) do (
        -- First we add m' to our M
        M = M | {m'};
        -- Then, we multiply out a power. 
        m' = m'* m;

        for j to (numgens R - 1) do (
            while (degree(modVars#j, m') > modDegree) do (
                m' = lift(m' / ((modVars#j)^(modDegree)), R);
            );
        );
    );

    --remove duplicate monomials
    M = unique(M);
    M = sort(M);

    -- remove nonminimal elements from the set
    for i to (#M - 1) do (
        -- Assume that it is initially minimal
        isMinimal := true;

        -- Then loop through all the elements that aren't i;
        for j to (i-1) when isMinimal do (

            -- Then, we create a variable to check if M#i is always a higher degree than M#j
            higherDegree = true;
            for var in (gens R) when higherDegree do (
                if (degree(var, M#i) < degree(var, M#j)) then (
                    higherDegree = false;
                );
            );
            -- If it is a higher degree, then it's not minimal
            if higherDegree then isMinimal = false;
        );

        -- However, if it is minimal, then we add it to the minimal set, N.
        if (isMinimal) then (N = N | {M#i};);
    );

    -- Finally, we return n
    return N;
)


export {"genseeds"}

beginDocumentation()

document {
    Key => genseeds, 

    Headline => "List the seeds that can be used to generate a complete ring of invariants.",

    Usage => "genseeds(R,W,ZList)",

    Inputs => {
        "R" => PolynomialRing => {"Ring upon which the group action is applied."},
        "W" => Matrix => {"The weight matrix representing the group action."},
        "ZList" => List => {"The list of the dimensions of the weight matrix."}
        },    

    Outputs => {
        "Seeds" => List => {"A list of the generating seeds."},
    },

    EXAMPLE {
        "R = QQ[x_0..x_2];",
        "W = matrix{{0,1,1},{1,0,1}};",
        "ZList = {2,2};",
        "genseeds(R,W,ZList)"
    },

    PARA {"This function returns the generating seeds that can be used to find all invariants under a group action."}

}
