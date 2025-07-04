-- Input: A matrix
-- Output: Boolean that gives whether the matrix defines a nondegenerate symmetric bilinear form over a field of characteristic not 2

isWellDefinedGWOverAlgebra = method()
isWellDefinedGWOverAlgebra Matrix := Boolean => M -> (
    
    -- Return false if the matrix isn't square and symmetric
    if not isSquareAndSymmetric M then return false;

    -- Return false if the matrix represents a degenerate form
    if isDegenerate M then return false;

    -- Return false if the matrix isn't defined over a field
    if not (isField ring M or (instance(ring M, QuotientRing) and isField coefficientRing ring M and dim ring M == 0)) then return false;
    
    -- Returns false if the matrix is defined over a field of characteristic 2
    if char(ring M) == 2 then return false;

    -- Otherwise, return true
    true
    )

-- We define GrothendieckWittClass to be a new type, meant to represent the isomorphism class 
-- of a nondegenerate symmetric bilinear form over a field of characteristic not 2

GrothendieckWittClassOverAlgebra = new Type of HashTable
GrothendieckWittClassOverAlgebra.synonym = "Grothendieck-Witt Class over algebra"

-- Input: A matrix M representing a nondegenerate symmetric bilinear form over a field of characteristic not 2
-- Output: The GrothendieckWittClass representing the symmetric bilinear form determined by M

makeGWClassOverAlgebra = method()
makeGWClassOverAlgebra Matrix := GrothendieckWittClassOverAlgebra => M -> (
   if isWellDefinedGWOverAlgebra M then (
        new GrothendieckWittClassOverAlgebra from {
            symbol matrix => M,
            symbol cache => new CacheTable,
            }
        )
    else (
        error "makeGWClass called on a matrix that does not represent a nondegenerate symmetric bilinear form over a field of characteristic not 2";
	)
    )

getRing = method()
getRing GrothendieckWittClassOverAlgebra := Ring => beta -> (
    ring getMatrix beta
    )

getBaseField GrothendieckWittClassOverAlgebra := Ring => beta -> (
    toField ring getMatrix beta
    )

-- Input: A GrothendieckWittClass representing a symmetric bilinear form determined by a matrix M
-- Output: The matrix M

getMatrix GrothendieckWittClassOverAlgebra := Matrix => alpha -> (
    alpha.matrix
    )

-- Input: A Grothendieck-Witt class beta over QQ, RR, CC, or a finite field of characteristic not 2
-- Output: A diagonalized form of beta, with squarefree entries on the diagonal

getDiagonalClassOverAlgebra = method()
getDiagonalClassOverAlgebra GrothendieckWittClassOverAlgebra := GrothendieckWittClassOverAlgebra => beta -> (

    -- Check if the diagonal class has already been computed; if so, recall it from the cache
    if beta.cache.?getDiagonalClassOverAlgebra then return beta.cache.getDiagonalClassOverAlgebra;

    getDiagonalClassOfBetaMatrix := diagonalizeViaCongruenceOverAlgebra getMatrix beta;

    -- The computed diagonal class gets stored in the cache
    beta.cache.getDiagonalClass = makeGWClassOverAlgebra getDiagonalClassOfBetaMatrix;
    makeGWClassOverAlgebra getDiagonalClassOfBetaMatrix
    )

diagonalizeViaCongruenceOverAlgebra = method()
diagonalizeViaCongruenceOverAlgebra Matrix := Matrix => AnonMut -> (
    k := ring AnonMut;
    --if not isField k then error "expected matrix over a field";
    kk := toField k;

    --Warning: this might work too generally!

    if not isSquareAndSymmetric AnonMut then
	error "matrix is not symmetric";
    
    -- If the matrix is already diagonal, then return it
    if isDiagonal AnonMut then return AnonMut;
    
    -- Otherwise, we iterate through positions below the diagonal, performing row operations followed by the corresponding
    -- column operations in order to obtain a diagonal matrix congruent to the original
    B := sub(AnonMut, kk);
    A := mutableMatrix B;

    n := numRows A;
    for col from 0 to n - 1 do (
	-- If diagonal entry in column "col" is zero
        if A_(col,col) == 0 then (
            for row from col + 1 to n - 1 do ( 
		-- Scan for nonzero entries in column "col" below the diagonal entry
                if A_(row,col) != 0 then (
                    if A_(row,row) == 0 then (
		        -- Row operation to make A_(col,col) nonzero
                        rowAdd(A, col,1,row);
		        -- Column operation to keep reduced matrix congruent to original matrix
                        columnAdd(A, col,1,row);
                        )
                    else (
		        -- Row and column swaps to make A_(col,col) nonzero
                        rowSwap(A, col,row);
                        columnSwap(A, col,row);
                        );
                    break;
                    );
                );
            );
        -- Now A_(col,col) != 0 unless there was a zero row/column; we use it to clear column "col" below this entry
        if A_(col,col) != 0 then (
            for row from col + 1 to n - 1 do (
                temp := A_(row,col);
                -- Row operation to make A_(row,col) zero
                rowAdd(A, row, -temp/A_(col,col), col);
	        -- Column operation to keep reduced matrix congruent
                columnAdd(A, row, -temp/A_(col,col), col);
                );
            );
        );
    sub(matrix A,k) 
    )

getDiagonalEntriesOverAlgebra = method()
getDiagonalEntriesOverAlgebra GrothendieckWittClassOverAlgebra := List => beta -> (
    
    M := diagonalizeViaCongruenceOverAlgebra getMatrix beta;
    n := numRows M;
    L := {};
    
    for i from 0 to n - 1 do
	L = append(L, M_(i,i));
    L
    )
