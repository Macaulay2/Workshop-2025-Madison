document {
    Key => {transferGW, (transferGW, GrothendieckWittClass)},
	Headline => "the transfer of Grothendieck-Witt from an extension field to a base field",
	Usage => "transferGW(beta)",
	Inputs => {
	    GrothendieckWittClass => "beta" => {"Grothendieck-Witt class over a field of characteristic not 2"}
	    },
	Outputs => {
	    GrothendieckWittClass => {"the image of the Grothendieck-Witt class ", TT "beta", "in", TEX///$\text{GW}(k)$///, " under the canonical transfer map."}
	    },
	PARA {"Given a field extension ", TEX///$L/k$///, " and a Grothendieck-Witt class, ", TT "beta", " in ", TEX///$\text{GW}(L)$///, " computes the image of ", TT"beta", " under the canonical map ", TEX///$\text{GW}(L)\to\text{GW}(k)$///, " for a field extension ", TEX///$L/k$///, "."},
	EXAMPLE lines ///
		 beta = makeGWuClass matrix(QQ, {{0,2},{2,0}});
		 getScalar beta
	 	 ///,
    SeeAlso => {"GrothendieckWittClass", "algebraicTrace"}
        }