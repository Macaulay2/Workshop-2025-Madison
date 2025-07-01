--basis for Lie algebra A1
X1 = matrix{{1,0},{0,-1}};
X2 = matrix{{0,1},{0,0}};
X3 = matrix{{0,0},{1,0}};

--Killing form for A1
n = 1;
Killing = pair -> (A = pair_0; B = pair_1; pair = 2*(n+1)*trace(A*B));

--dual basis
Y1 = (1/8)*X1;
Y2 = (1/4)*X3;
Y3 = (1/4)*X2;

--commutator bracket
bracket = pair -> (
X = pair_0; Y = pair_1;
pair = X*Y - Y*X
);

--x vars for 2x2 matrices; a vars for coefficients of polys
R = QQ[x_(0,0) .. x_(1,1)];

--map  R -> R of Lie algebra commutator action on 2x2 matrices
Act = X -> (
M = transpose genericMatrix(R,2,2);
M = bracket(X,M);
X = map(R,R,{M_(0,0),M_(0,1),M_(1,0),M_(1,1)})
);

--degree 2 homogeneous polys over R
f = x_(0,0)*x_(1,1) - x_(0,1)*x_(1,0);


Cas = f -> (
    f1 = 0; f2 = 0; f3 = 0;

    A1 = Act(X1); B1 = Act(Y1);
    for i to 3 do (for j from i to 3 do (
        f1 = f1 + contract(R_i*R_j,f)*((A1 B1 R_i)*R_j + R_i*(A1 B1 R_j) + (A1 R_i)*(B1 R_j) + (B1 R_i)*(A1 R_j)) 
    ));

    A2 = Act(X2); B2 = Act(Y2);
    for i to 3 do (for j from i to 3 do (
        f2 = f2 + contract(R_i*R_j,f)*((A2 B2 R_i)*R_j + R_i*(A2 B2 R_j) + (A2 R_i)*(B2 R_j) + (B2 R_i)*(A2 R_j)) 
    ));

    A3 = Act(X3); B3 = Act(Y3);
    for i to 3 do (for j from i to 3 do (
        f3 = f3 + contract(R_i*R_j,f)*((A3 B3 R_i)*R_j + R_i*(A3 B3 R_j) + (A3 R_i)*(B3 R_j) + (B3 R_i)*(A3 R_j)) 
    ));

    f = f1 + f2 + f3
);

-- Reynolds demo 
V = QQ^2

v1 = matrix{{0},{2}}
v2 = matrix{{1},{4}}




