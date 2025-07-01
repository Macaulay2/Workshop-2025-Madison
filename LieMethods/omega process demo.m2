-- ***********************************************
-- Macaulay2 SLreynolds via the Omega process
-- (Based on the Singular)
-- ***********************************************

-- Set the size of the matrices.
n = 3;

------------------------------------------------------------
-- Step 1. Define the ring P in the x–variables.
P = QQ[x_1..x_n];
gens P

-- Step 2. Define the ring in the g–variables.

nn = n^2;

------------------------------------------------------------
-- Step 3. Construct the ring R = QQ[x_1,…,x_n, g_11,…,g_nn].
R = QQ[x_1..x_n, g_1..g_nn];
gensR = gens R

-- Split the generators into the x–variables and the g–variables.
xVars = gensR_{0..(n-1)};
gVars = gensR_{n..(n+n^2-1)};

------------------------------------------------------------
-- Step 4. Build the generic n×n matrix M for the SL_n–representation.

M = genericMatrix(R, g_1, n,n);
M


------------------------------------------------------------
-- Step 5. Define the group ideal for SLₙ: ideal(det(M) - 1).
gideal = ideal(det M - 1);
gideal

------------------------------------------------------------
-- Step 6. Define the substitution map act : P --> r.
-- First form the column vector v = (x₁, …, xₙ)^T.
vx = genericMatrix(R, x_1, n, 1);
vx
-- w = (transpose(vx))*M, shape 1xn
w = (transpose(vx)) * M;
w
-- Define the ring map act : P --> R.
act = map(R, P, w);
act

------------------------------------------------------------
-- Step 7. Define the “projection” map gtozero: R --> P forgetting g variables.
gtozero = map(P, R, join(gens P, nn:0));

------------------------------------------------------------
-- Step 8. Omega process
-- omegaProcess(h, q) is the result of differentiating h
-- with respect to g_(q), g_(q-1), …, g_(1)
omegaProcess = (h, q) -> (
    hNew := h;
    for i from 1 to q do (
       hNew = diff(g_(q-i+1),hNew);
    );
    hNew
);


------------------------------------------------------------
-- Step 9. Define the SLreynolds function.
SLreynolds = f -> (
    h = act(f)%gideal; --****fix this!!! need to kill gs
    --print h;
    q = n - 1;
    output = 0_P;  -- initialize output (in P, since gtozero maps R to P)
    c = 1;         -- normalization constant
    iVal = 1;      -- loop counter
    
    while (h != 0) do (
       output = output + gtozero(h)/c;
       h = omegaProcess(h, q);
       -- update c by the product of iVal, iVal+1, …, iVal+q-1.
       for j from iVal to (iVal+q-1) do (
          c = c * j;
       );
       iVal = iVal + 1;
    );
    output
);

------------------------------------------------------------
-- Step 10. Examples
use P;
f1 = (x_1^2 + x_2^2 + x_3^2)_P; 
result1 = SLreynolds(f1);
print result1;

f2 = ((x_1 * x_2)^2)_P;
result2 = SLreynolds(f2);
print result2;
