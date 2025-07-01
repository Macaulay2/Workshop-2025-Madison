------------------------------------------------------------------------
--  Documentation for getActionMatrix
------------------------------------------------------------------------

document
{
     Key       => getActionMatrix,
     Headline  => "construct the action matrix of a ring element from the elimination template",

     Usage     => "A = getActionMatrix(actVar, mp, M)",

     Inputs    =>
     {
         "actVar" => {"action variable whose action on the quotient ring one wishes to encode."},

         "mp"     => {"a triple (E,R,B) of lists of monomials produced by a border–basis routine; here\n",
                      "      *   E : \"excessive\" monomials to be eliminated,\n",
                      "      *   R : \"reducible\" monomials expressible through the basis,\n",
                      "      *   B : ordered list of \"basic\" monomials spanning the quotient."},

         "M"      => {"Elimination template matrix of size (|E|+|R|) × (|R|+|B|),\n",
                      "      where the rows correspond to the monomials in E and R (in that order),\n",
                      "      and the columns correspond to the monomials in R and B (in that order).\n",
                      "      The matrix must have full numerical rank."}
     },

     Outputs   =>
     {
         "A"      => {"a real square matrix whose size equals |B|;  `A_(i,j)` is the coefficient of the j-th basic monomial in the expansion of *actVar*·B_i in the quotient ring.  Columns corresponding to basis monomials **not divisible** by *actVar* form an appended identity block, ensuring that the column order matches `B` exactly."}
     },

     Description =>
     {
         TEXT
         "The routine realises the standard linear-algebraic construction of an\n",
         "*action matrix* used in eigenvalue / border-basis root-finding algorithms (see, e.g., Sommese & Wampler, *The Numerical Solution of Systems of Polynomials*).\n\n",

         "Algorithmic steps\n",
         "-----------------\n",
         "1.  **Row elimination of excessive monomials** – the sub-matrix of *M* corresponding to *E* is factored by an LU-decomposition (package *NumericalLinearAlgebra*) and removed from the system via suitable row operations; this is numerically cheaper than a full RREF.\n",
         "2.  **Resolution of reducible monomials** – for the remaining rows the linear system `Mr * A = Mb` is solved, expressing each reducible monomial in *R* as a linear combination of the basic monomials *B*.\n",
         "3.  **Identity block for “extra” basis elements** – basic monomials not divisible by `actVar` contribute an identity block on the right, giving the final action matrix the correct shape and ordering.\n\n",

         "The method assumes that the block `Mr` is nonsingular; if the input data stem from a border basis of a *zero-dimensional* ideal this condition is automatically satisfied in exact arithmetic.  For badly-conditioned numerical data you may wish to add pivot-tolerance tests."
     },

     Caveat    =>
     {
         "• The matrix *M* must have full numerical rank; no rank-revealing pivoting beyond the LU factorisation is performed.\n",
         "• All arithmetic is carried out in double precision (ring RR).  Replace calls to `LUdecomposition` and `solve` with their exact counterparts if exact arithmetic is required."
     },

     SeeAlso   => {LUdecomposition, solve},

     Examples  => PRE
     "
     needsPackage \"NumericalLinearAlgebra\"

     -- a toy system ---------------------------------------------------
     R  = QQ[x,y];
     E  = {x^3, y^3};
     Rm = {x^2*y, x*y^2};
     B  = {1, x, y, x^2, y^2};
     mp = (E, Rm, B);

     -- random dense template of size (|E|+|Rm|) × (|Rm|+|B|)
     M  = random(RR^(#E+#Rm), RR^(#Rm+#B));

     A  = getActionMatrix(x, mp, M);

     -- check: size of A equals |B| × |B|
     assert(numrows A == #B and numcols A == #B)
     "
}
