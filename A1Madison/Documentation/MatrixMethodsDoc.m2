doc ///
	Key
		diagonalizeViaCongruence
		(diagonalizeViaCongruence, Matrix)
		[diagonalizeViaCongruence, linearTolerance]
	Headline
		diagonalizes a symmetric matrix via congruence
	Usage
		diagonalizeViaCongruence M
	Inputs
		M : Matrix
			a symmetric matrix over any field or finite étale algebras over a field
		linearTolerance => RR
			a positive number specifying the tolerance to which entries are considered zero
	Outputs
		: Matrix
			a diagonal matrix congruent to @TT("M")@
	Description
		Text
			Given a symmetric matrix @TT("M")@ over any field or finite étale algebra over a field, this command gives a diagonal matrix congruent to @TT("M")@. Note that the order in which the diagonal terms appear is not specified. 
		Example
			R = QQ[x]/(x^2 - 1)
			M = matrix(R, {{1,2},{2,x}});
			diagonalizeViaCongruence M
	Caveat
		When computing over inexact fields such as $\mathbb{R}$ or $\mathbb{C}$, the @TT("linearTolerance")@ option specifies the tolerance to which entries are considered zero. The default tolerance is $10^{-12}$.
	SeeAlso
		getDiagonalClass
		getDiagonalEntries
///