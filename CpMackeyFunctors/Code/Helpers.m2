matrixPower = method()
matrixPower(Matrix,ZZ) := Matrix => (M,n) -> (
    if n < 0 then error "Exponent must be non-negative";
    result := M;
    for i from 2 to n do (
        result = result * M;
    );
    return result
)