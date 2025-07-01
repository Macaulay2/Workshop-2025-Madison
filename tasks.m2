-- Tasks and projects:

-- 1. Rational Polytopes
-- Discussion on Package Design, names of functions
-- Additions to Normaliz package, how should they be incorporated? Pull request to Noraliz
-- Examples, Documentation
-- Add function toricVariety, toricIdeal of a rational polytope

-- 2. Weighted Ehrhart Theory
-- Examples
-- Implement Theorem 1.1 from https://arxiv.org/pdf/2402.11328 [Sums of Weighted Lattice Points of Polytopes]
-- Use Normaliz to compute weighted Ehrhart series, see [Sections 5.2.12, 7.9, Normaliz Documentation https://github.com/Normaliz/Normaliz/blob/master/doc/Normaliz.pdf]
-- Work out how to do the most general version of weighted Ehrhart theory

-- 3. Equivariant Ehrhart Theory
-- Goal: compute the Equivariant Ehrhart Series for Polytope under cyclic group action
-- Rough outline:
-- > Check if a polytope is invariant under action of matrix
-- > Get the fixed Polytopes and their Ehrhart Series
-- > Compute the Character Table for the group
-- > put result into the right form
-- Package the equivariant Ehrhart file into a package
