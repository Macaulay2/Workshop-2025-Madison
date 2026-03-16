restart
-*
Problem: given the following data

  (1) Two lines L_1, L_2 in PP^3 (projective space) which intersect
  (2) Two matching lines l_1, l_2 in PP^2 (projective plane)
      (L_1, l_1) is a match
      (L_2, l_2) is a match

Goal: recover a linear projection P : PP^3 --> PP^2
      such that P(L_i) = l_i for i = 1, 2,
        with matrix representation
	  P = diag(f, f, 1) * (R | 0),
	    R R^T = I, det(R) = 1

(This describes a rotating camera with variable zoom.)
4 DOF
  3 DOF for R
  1 DOF for f
*-

-- create a synthetic problem-solution pair
-- point where 3D lines intersect
FF = QQ
a = random(FF^3, FF^1) || matrix{{1}}
-- other 2 points defining two 3D lines
-- L_1 = <a, b1>, L_2 = <a, b2>
b1 = random(FF^3, FF^1) || matrix{{1}}
b2 = random(FF^3, FF^1) || matrix{{1}}
-- this method function creates the quaternion parametrization of SO3
-*
Q : PP^3 --> P(C^(3\times 3)) = PP^8

  [w:x:y:z] -> 3 x 3 matrix defined in function below

Claim: The Zariski closure of the image of Q is the variety of
       3 x 3 rotation matrices, i.e.
       for R = Q([w:x:y:z]), we have
       R R^T = I, det (R) = 1,
       and conversely "most points" satisfying these equations
       come from some choice of [w:x:y:z]

Q is the quaternion parametrization of SO(3)

Note: The map Q is (generically) 1-to-1, but usually people lift
 to the sphere S^3, which is a double-cover of P^3
*-
Q2R = method(Options => {Normalized=>false, FF=>FF})
Q2R (Thing, Thing, Thing, Thing) := o -> (w, x, y, z) -> (
    M := matrix{
	{w^2+x^2-y^2-z^2, 2*x*y-2*w*z, 2*w*y+2*x*z},
	{2*x*y+2*w*z, w^2-x^2+y^2-z^2, -2*w*x+2*y*z},
	{-2*w*y+2*x*z, 2*w*x+2*y*z, w^2-x^2-y^2+z^2}
	};
    if o.Normalized then (1/(w^2+x^2+y^2+z^2)) * M else M
    )
Q2R List := o -> L -> Q2R(L#0, L#1, L#2, L#3, o)

-- unknown camera parameters: goal is to recover these
(w0, x0, y0, z0, f0) = (random FF, random FF, random FF, random FF, random FF)
R0 = Q2R(w0, x0, y0, z0, Normalized=>true)
-- check this a rotation matrix
R0 * transpose R0, det R0

-- generate "ground truth" camera matrix
P0 = diagonalMatrix{f0, f0, 1} * (R0 | matrix{{0},{0},{0}})
-- implicit equations of 2D lines obtained by projection
l1 = gens ker transpose(P0 * (a | b1))
l2 = gens ker transpose(P0 * (a | b2))
-*
Goal: recover P0 just from the data of (L1,L2,l1,l2)
*-

S = FF[w..z,f]
R = Q2R(w,x,y,z) -- _SCALED_ rotation matrix

P = diagonalMatrix{f,f,1} * (R | matrix{{0},{0},{0}})
-- 5 equations in 5 unknowns f, w, x, y, z
I = ideal(
    transpose l1 * P * a,
    transpose l1 * P * b1,
    transpose l2 * P * a,
    transpose l2 * P * b2,
    w^2+x^2+y^2+z^2-1
    )
dim I, degree I, radical I == I
needsPackage "EigenSolver"
minimalProblemSolutions = zeroDimSolve I
-- did we recover the ground-truth solution
groundTruthSolution = matrix{(1/sqrt(w0^2+x0^2+y0^2+z0^2)*{w0,x0,y0,z0})|{f0}}
-- yes! some ways to verify this below
select(minimalProblemSolutions, x -> norm(matrix x- groundTruthSolution) < 1e-10)
position(minimalProblemSolutions, x -> norm(matrix x- groundTruthSolution) < 1e-10)
-*
IDEA:

 Add this example to paper: Start in Intro, End in Section 4?

 Intro: Explain problem and set up notation

 Section 4: Explain solution ingredients
   (1) getTemplate (solve for initial data (L1, L2, l1, l2)
   (2) copyTemplate (solve for new initial data, "with less pain") 
*-
