-*
   Copyright 2025, ??????
    
   You may redistribute this file under the terms of the GNU General Public
   License as published by the Free Software Foundation, either version 2 of
   the License, or any later version.
*-

beginDocumentation()

///document {
 Node
  Key => {PermutationGroup}
  Headline => "a package to deal with permutation groups and orbit sums",
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
    shuffMon
	orbSum
    orbSumList
}

document {
    Key => {ListSpInd},
    Headline => "List the special monomials that are related to a degree.",
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

}

document {
    Key => {ListSpMon}
    Headline => "List the special monomials that are related to a degree.",
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

}



document {
  Key => {orbSum}
  Headline => "Computes the orbit sum of a monomial.",
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
}

document {
  Key => {orbSumList}
  Headline => "Computes the orbit sums of a list of monomials.",
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

}///

document {
  Key => {shuffMon},
  Headline => "Permutes monomial.",
  Usage => "shuffMon(f,n)",
  Inputs => {
      "f" => Monomial => {"a monomial in the polynomial ring"},
      "n" => Number => {"the number of variables in the polynomial ring"},
      },

  Outputs => {
      "Mon" => List => {"a list of monomials"},
      },

  PARA {
    "This function takes a monomial and permutes all the variables of the monomial and puts all permutations in a list."
    },

   Example {
    "R = QQ[x,y]",
    "g = {x^2, x*y, y^2}",
    "S = subring g",
    "numgens presentationRing S"
    },
}

document {
	Key => {action, (action, RingOfInvariants)},
	
	Headline => "the group action that produced a ring of invariants",
	
	Usage => "action S",
	
	Inputs => {
	    	"S" => RingOfInvariants => {"of the group action being returned"},
		},
	
	Outputs => {
		GroupAction => {"the action that produced the ring of invariants in the input"}
		},
	"This function is provided by the package ", TO InvariantRing,".",
	
	PARA {
	    "This example shows how to recover the action of a
	    torus that produced a certain ring of invariants.
	    Note other action types are possible and may produce
	    a different looking output."
	    },
    	
	EXAMPLE {
		"R = QQ[x_1..x_4]",
		"T = diagonalAction(matrix {{0,1,-1,1},{1,0,-1,-1}}, R)",
		"S = R^T",
		"action S"
		},
	    }