-*
   Copyright 2025, ??????
    
   You may redistribute this file under the terms of the GNU General Public
   License as published by the Free Software Foundation, either version 2 of
   the License, or any later version.
*-

document {
  Key       => OrbitSum,
  Headline  => "Invariance of permutation groups.",
  BOLD "OrbitSum", " is a package dedicated to permutation groups and orbit sums.",
  PARA {"Current algorithms include:"},
  UL {
      TO listSpInd,
      TO listSpMon,
      TO orbSum,
      TO orbSumList,
      TO shuffMon
  },
  Contributors => "Gordon Novak & Sasha Arasha, for the documentation.",
  Acknowledgement => "A special thanks to the organizers of the 2025 Madison workshop",
  Caveat => {"warning"}
}

document {
  Key => listSpInd, 

  Headline => "List the special indicies that are related to a degree.",

  Usage => "listSpInd(n,d)",

  Inputs => {
    "n" => Number => {"the number of variables in the polynomial ring"},
        
    "d" => Number => {"the degree of the polynomial ring"}
    },    

  Outputs => {
    "M" => List => {"a list of special indicies"},
  },

  EXAMPLE {
    "listSpInd(3,4)",
    "listSpInd(1,2)",
    "listSpInd(4,5)"
  },

  "This function returns all the special incidies of degree ", CODE "d", " in ", CODE "n", " variables."

}

document {
  Key => listSpMon,

  Headline => "Lists the special monomials that are related to a degree",

  Usage => "listSpMon(n,d)", 
      
  Inputs => {
      "n" => Number => {"the number of variables in the polynomial ring"},
      "d" => Number => {"the degree of the polynomial ring"},
  },
  Outputs => {
      List => {"a list of special monomials"},       
  },

  EXAMPLE {
    "listSpMon(3,4)",
    "listSpMon(1,2)",
    "listSpMon(4,2)"
  },

  "This function returns all the special monomials of degree d in n variables."
}

document {
  Key => orbSumList,

  Headline => "Computes the orbit sums of a list of monomials.",

  Usage => "orbSumList(G,n,d)",

  Inputs => {
      "G" => {"a permutation group"},
         
      "n" => Number => {"the number of variables in the polynomial ring"},
        
      "d" => Number => {"the degree of the polynomial ring"}, 
  },
  Outputs => {
      List  => {"a list of orbit sums of special monomials"},
  },
  "This function computes the orbit sums of a list of special monomials of degree d in n variables under the action of a permutation group G.",

  EXAMPLE {
    "needsPackage(\"SpechtModule\")",
    "G = generatePermutationGroup {{1,0,2,3},{1,2,3,0}};",
    "n = 3",
    "d = 3",
    "orbSumList(G,n,d)"
  }
}

document {
  Key => orbSum, 

  Headline => "Computes the orbit sum of a monomial",

  Usage => "orbSum(r,G,n)",

  Inputs => {
    "r" => RingElement => {"a monomial in the polynomial ring"},
        
    "G" => {"a permutation group"},     

    "n" => Number => {"the number of variables in the polynomial ring"},
  },    

  Outputs => {
    "g" => RingElement => {"the orbit sum of the monomial f"},
  },

  PARA {
    "This function computes the orbit sum of a monomial f under the action of a permutation group G."
  },
  
  EXAMPLE {
    "needsPackage(\"SpechtModule\")",
    "R = QQ[x..z]",
    "r = x^2*y",
    "G = generatePermutationGroup {{1,0,2,3},{1,2,3,0}};",
    "orbSum(r,G,numgens R)"
  }
   
}

document {
  Key => shuffMon, 

  Headline => "Permutes monomial.",

  Usage =>  "shuffMon(r,n)",

  Inputs => {
    "r" => RingElement => {"a monomial in the polynomial ring"},
      
    "n" => Number => {"the number of variables in the polynomial ring"}, 
  }, 

  Outputs => {
    "Mon" => List => {"a list of monomials"},
  }, 

  PARA {"This function takes a monomial and permutes all the variables of the monomial and puts all permutations in a list."}, 

  EXAMPLE {
    "R = QQ[x..z];",
    "r = x^2*y;",
    "shuffMon(r,numgens R)"
  }
}


