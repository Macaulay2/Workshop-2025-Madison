-*
   Copyright 2023, ??????
    
   You may redistribute this file under the terms of the GNU General Public
   License as published by the Free Software Foundation, either version 2 of
   the License, or any later version.
*-

beginDocumentation()

doc ///
 Node
  Key
   Permutation Group
  Headline
   a package to deal with permutation groups and orbit sums
  Description
   Text
     is a package to give basic ring permutations and orbit sums
   Example
     needsPackage("PermutationGroup")
   Text
     In particular, this package allows you to do computations with permutation groups and orbit sums.
  Subnodes
    ListSpInd
    ListSpMon
    ShuffMon
	orbSum
    orbSumList
///

doc ///
    Key
        ListSpInd
    Headline
        List the special monomials that are related to a degree.
    Usage
        ListSpInd(n,d)
    Inputs
        n: Number
            the number of variables in the polynomial ring
        d: Number
            the degree of the polynomial ring
    Outputs
        M: List
            a list of special monomials
    Description
        Text
            This function returns all the special monomials of degree d in n variables.
        Example

///

doc ///
    Key
        ListSpMon
    Headline
        List the special monomials that are related to a degree.
    Usage
        ListSpMon(n,d)
    Inputs
        n: Number
            the number of variables in the polynomial ring
        d: Number
            the degree of the polynomial ring
    Outputs
        M: List
            a list of special monomials
    Description
        Text
           This function returns all the special monomials of degree d in n variables.
        Example

///

doc ///
  Key
   ShuffMon
  Headline
    Permutes monomial.
  Usage
    ShuffMon(f,n)
  Inputs
   f: Monomial
       a monomial in the polynomial ring
    n: Number
       the number of variables in the polynomial ring
  Outputs
    Mon: List
         a list of monomials
  Description
   Text
    This function takes a monomial and permutes all the variables of the monomial and puts all permutations in a list.
   Example
    R = QQ[x,y];
    g = {x^2, x*y, y^2};
    S = subring g;
    numgens presentationRing S
///

doc ///
  Key
   orbSum
  Headline
   Computes the orbit sum of a monomial.
  Usage
    orbSum(f,G,n)
  Inputs
    f: Monomial
         a monomial in the polynomial ring
     G: Group
         a permutation group
     n: Number
         the number of variables in the polynomial ring
  Outputs
    g: Monomial
         the orbit sum of the monomial f
  Description
   Text
    This function computes the orbit sum of a monomial f under the action of a permutation group G.
   Example

///

doc ///
  Key
   orbSumList
  Headline
    Computes the orbit sums of a list of monomials.
  Usage
    orbSumList(G,n,d)
  Inputs
    G: Group
         a permutation group
    n: Number
        the number of variables in the polynomial ring
    d: Number
        the degree of the polynomial ring
  Outputs
    L: List
        a list of orbit sums of special monomials.
  Description
   Text
    This function computes the orbit sums of a list of special monomials of degree d in n variables under the action of a permutation group G.
   Example

///
