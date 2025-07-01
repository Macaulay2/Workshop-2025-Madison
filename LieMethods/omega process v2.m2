needsPackage "InvariantRing";

SLreynolds = method()
SLreynolds (RingElement) := f -> (

    if not isHomogeneous f then
        error "f must be homogeneous";

    P = ring f;

    n = numgens P;

    xVars = {x_1..x_n};
    gVars = {g_(1,1)..g_(n,n)};

    R = QQ[join(xVars, gVars)];

    M = genericMatrix(R, g_(1,1), n,n);

    gideal = ideal(det M - 1);

    vx = genericMatrix(R, x_1, n, 1);

    w = (transpose(vx)) * M;

    act = map(R, P, w);

    gtozero = map(P, R, join(gens P, n^2:0));

    omegaProcess = (h, q) -> (
        hnew := h;
        -- for k=1..q, differentiate w.r.t. g_{ row = q-k+1, col = k }
        for k from 0 to q-1 do (
            row := q - k;
            print(row);
            col := k + 1;
            print(idx);
            hnew = diff(g_(row, col), hnew);
            print(hnew);
        );
        hnew
    );
    

    h = act(f)%gideal;
    print(h);
    q = n - 1;
    output = 0_P;  -- initialize output (in P, since gtozero maps R to P)
    c = 1;         -- normalization constant
    iVal = 1;      -- loop counter
    
    print("Starting loop");

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


n = 3;

P = QQ[x_1..x_n];

f = x_1*x_2*x_3^4 + x_1^2*x_2^3*x_3 + x_1^2*x_2^2*x_3^2 + x_1^4*x_2*x_3;

SLreynolds(f)

g = x_1 + x_3^4;

SLreynolds(g) --should raise error: not homogeneous!



