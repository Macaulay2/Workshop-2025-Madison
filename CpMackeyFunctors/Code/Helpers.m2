-- Inputs a matrix M and a number n, and outputs the matrix M^n
matrixPower = method()
matrixPower(Matrix,ZZ) := Matrix => (M,n) -> (
    if n < 0 then error "Exponent must be non-negative";
    if n < 2 then error "Exponent is too small";
    result := M;
    for i from 2 to n do (
        result = result * M;
    );
    return result
)