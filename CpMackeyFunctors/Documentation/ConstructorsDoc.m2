document{
    Key => {makeBurnsideMackeyFunctor, (makeBurnsideMackeyFunctor, ZZ)},
    Headline => "constructs the Burnside Mackey functor for the group",
    Usage => "makeBurnsideMackeyFunctor(p)",
    Inputs => {
	ZZ => "p" => {"a prime number", TEX///$p$///},
	},
    Outputs => {
	CpMackeyFunctor => {"the ", TEX///$C_p$///, "-Mackey functor ", TEX///$\underline{A}$///, "."},
	},
    PARA{
        "The ", EM "Burnside Mackey functor", " of a group ", TEX///$G$///, " is defined by sending ", TEX///$G/H$///, " to the Burnside ring ", TEX///$A(H)$///, ", with transfer and restriction coming from transfer and restriction of ", TEX///$G$///,"-sets. When ", TEX///$G$///, " is a cyclic group of prime order, this admits a nice form. The underlying module is given by ", TEX///$\mathbb{Z}$///, " with trivial conjugation action, while the fixed module is ", TEX///$\mathbb{Z} \{1,t\}$///, ". Restriction is defined as",
        TEX///$$\text{res}_e^{C_p} \colon \mathbb{Z}\{1,t\} \to \mathbb{Z}$$///,
        "by sending ", TEX///$t\mapsto p$///, ". Transfer is of the form ",
        TEX///$$\text{tr}_e^{C_p} \colon \mathbb{Z} \to \mathbb{Z}\{1,t\}$$///,
        "by sending ", TEX///$x\mapsto xt$///, "."
        },
    EXAMPLE lines///                                                                         
    makeBurnsideMackeyFunctor(5)
    ///,                  
    PARA{EM "Citations:"},
    UL{
	{"[S73] J.P. Serre, ", EM "A course in arithmetic,", " Springer-Verlag, 1973."},
	{"[MH73] Milnor and Husemoller, ", EM "Symmetric bilinear forms,", " Springer-Verlag, 1973."},
    },
    --SeeAlso => {"getHilbertSymbolReal", "getHasseWittInvariant"}
}