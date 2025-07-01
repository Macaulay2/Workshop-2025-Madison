document{
    Key => {makeBurnsideMackeyFunctor, (makeBurnsideMackeyFunctor, ZZ)},
    Headline => "constructs the Burnside Mackey functor",
    Usage => "makeBurnsideMackeyFunctor(p)",
    Inputs => {
	    ZZ => "p" => {"a prime number ", TEX///$p$///},
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
    PARA{
        "The Burnside Mackey functor represents the functor ", TEX///$\text{Mack}_{C_p} \to \text{Ab}$///, " which sends a Mackey functor ", TEX///$M$///, " to its fixed point level ", TEX///$M(C_p/C_p)$///, ". In other words, we have a natural isomorphism",
        TEX///$$\text{Hom}_{\text{Mack}_{C_p}}\left(\underline{A},M \right) \cong M(C_p/C_p).$$///
    },
    SeeAlso => {"makeFixedFreeMackeyFunctor"}
}


document{
    Key => {makeUnderlyingFreeMackeyFunctor, (makeUnderlyingFreeMackeyFunctor, ZZ)},
    Headline => "constructs the free Mackey functor on an underlying generator",
    Usage => "makeUnderlyingFreeMackeyFunctor(p)",
    Inputs => {
	    ZZ => "p" => {"a prime number ", TEX///$p$///},
	},
    Outputs => {
	    CpMackeyFunctor => {"the free Mackey functor on an underlying generator."},
	},
    PARA{
        "The free ", TEX///$C_p$///,"-Mackey functor on an underlying generator represents the functor ", TEX///$\text{Mack}_{C_p} \to \text{Ab}$///, " which sends a Mackey functor ", TEX///$M$///, " to its underlying level ", TEX///$M(C_p/e)$///, ". This means there is a natural isomorphism",
        TEX///$$\text{Hom}_{\text{Mack}_{C_p}}\left(\underline{B},M \right) \cong M(C_p/e).$$///, 
        "In components, the underlying module is the free module on the ", TEX///$C_p$///, "-set ", 
        TEX///$C_p/e=\{1,\gamma,\gamma^2,\ldots,\gamma^{p-1}\}$///, " with conjugation induced by the ", TEX///$C_p$///, "-action of left multiplication. The fixed module is the module ", TEX///$\ZZ$///, ". Restriction is the map",
        TEX///$$\text{res}_e^{C_p} \colon \ZZ \to \ZZ\{1,\gamma,\ldots,\gamma^{p-1}\}$$///,
        "by sending ", TEX///$1\mapsto 1+\gamma+\cdots+\gamma^{p-1}$///, ". The transfer is the map ",
        TEX///$$\text{tr}_e^{C_p} \colon \ZZ\{1,\gamma,\ldots,\gamma^{p-1}\} \to \ZZ$$///,
        "by sending ", TEX///$\gamma^i\mapsto 1$///, " for all ", TEX///$i$///,"."
    },
    EXAMPLE lines///
    makeUnderlyingFreeMackeyFunctor(5)
    ///
}


document{
    Key => {makeFixedPointMackeyFunctor, (makeFixedPointMackeyFunctor, ZZ, Matrix)},
    Headline => "constructs the fixed-point Mackey functor of a C_p-module",
    Usage => "makeFixedPointMackeyFunctor(p,C)",
    Inputs => {
	    ZZ => "p" => {"a prime number ", TEX///$p$///},
        Matrix => "C" => {"an order ", TEX///$p$///," automorphism ", TEX///$\text{c}$///},
	},
    Outputs => {
	    CpMackeyFunctor => {"the fixed-point Mackey functor of the ", TEX///$C_p$///, "-module specified by ", TEX///$\text{c}$///},
	},
    PARA{
        "Given a module ", TEX///$M$///, " with ", TEX///$C_p$///, "-action specified by an automorphism ", TEX///$\text{c}\colon M\to M$///, " there is a naturally associated ", EM"fixed-point Mackey functor ", TEX///$\text{FP}(M)$///, ". The underlying module is the module ", TEX///$M$///," with conjugation given by ", TEX///$\text{c}$///,". The fixed module is the module ", TEX///$M^{C_p}$///, " of fixed points of the action, represented as the kernel of ", TEX///$1-\text{c}$///,". The restriction is the map ",
        TEX///$$\text{res}_e^{C_p} \colon M^{C_p} \to M$$///,
        "by induced by the inclusion of fixed-points. The transfer is the map ",
        TEX///$$\text{tr}_e^{C_p} \colon M\to M^{C_p}$$///,
        "given by ", TEX///$1+\text{c}+\text{c}^2+\cdots+\text{c}^{p-1}$///, "."
    },
    EXAMPLE lines///
    makeFixedPointMackeyFunctor(2,matrix{{0,1,0},{1,0,0},{0,0,1}})
    ///
}


document{
    Key => {makeOrbitMackeyFunctor, (makeOrbitMackeyFunctor, ZZ, Matrix)},
    Headline => "constructs the orbit Mackey functor of a C_p-module",
    Usage => "makeOrbitMackeyFunctor(p,C)",
    Inputs => {
	    ZZ => "p" => {"a prime number ", TEX///$p$///},
        Matrix => "C" => {"an order ", TEX///$p$///," automorphism ", TEX///$\text{c}$///},
	},
    Outputs => {
	    CpMackeyFunctor => {"the orbit Mackey functor of the ", TEX///$C_p$///, "-module specified by ", TEX///$\text{c}$///},
	},
    PARA{
        "Given a module ", TEX///$M$///, " with ", TEX///$C_p$///, "-action specified by an automorphism ", TEX///$\text{c}\colon M\to M$///, " there is a naturally associated ", EM"orbit Mackey functor ", TEX///$\text{O}(M)$///, ". The underlying module is the module ", TEX///$M$///," with conjugation given by ", TEX///$\text{c}$///,". The fixed module is the quotient module ", TEX///$M_{C_p}$///, " of the action, represented as the cokernel of ", TEX///$1-\text{c}$///,". The restriction is the map ",
        TEX///$$\text{res}_e^{C_p} \colon M_{C_p} \to M$$///,
        "given by ", TEX///$1+\text{c}+\text{c}^2+\cdots+\text{c}^{p-1}$///, ". The transfer is the map",
        TEX///$$\text{tr}_e^{C_p} \colon M\to M_{C_p}$$///,
        "induced by the projection to the quotient."
    },
    EXAMPLE lines///
    makeOrbitMackeyFunctor(2,matrix{{0,1,0},{1,0,0},{0,0,1}})
    ///
}


document{
    Key => {makeZeroMackeyFunctor, (makeZeroMackeyFunctor, ZZ)},
    Headline => "constructs the zero Mackey functor for the group",
    Usage => "makeZeroMackeyFunctor(p)",
    Inputs => {
	    ZZ => "p" => {"a prime number ", TEX///$p$///},
	},
    Outputs => {
	    CpMackeyFunctor => {"the zero Mackey functor."},
	},
    PARA{
        "Perhaps the easiest Mackey functor is the ", EM "zero Mackey functor", " which has the zero-module as both the underlying and fixed modules. This is the zero object in the abelian category of Mackey functors, and is an important object to have for homological algebra."
    },
    EXAMPLE lines///
    makeZeroMackeyFunctor(2)
    ///
}

document{
    Key => {makeComplexRepresentationMackeyFunctor, (makeComplexRepresentationMackeyFunctor, ZZ)},
    Headline => "constructs the complex representation Mackey functor",
    Usage => "makeComplexRepresentationMackeyFunctor(p)",
    Inputs => {
	    ZZ => "p" => {"a prime number ", TEX///$p$///},
	},
    Outputs => {
	    CpMackeyFunctor => {"the complex representation Mackey functor for the group ", TEX///$C_p.$///},
	},
    PARA{
        "The ", EM "complex representation", " Mackey functor of a group ", TEX///$G$///, " is defined by sending ", TEX///$G/H$///, " to the Grothendieck group of complex ", TEX///$G$/// , "-representations. The transfer and restriction come from induction and restriction of ", TEX///$G$///,"-representations. When ", TEX///$G$///, " is a cyclic group of prime order, this admits a nice form. The underlying module is given by ", TEX///$\mathbb{Z}$///, " with trivial conjugation action, while the fixed module is ", TEX///$\mathbb{Z}\{\lambda_{0},\lambda_1,\dots,\lambda_{p-1}\}$///, ". The element ", TEX///$\lambda_i$///,  " represents the one dimensional complex representation given by multiplication by ",TEX///$e^{2\pi i/p}.$///, " Restriction is defined as",
        TEX///\[\text{res}_e^{C_p} \colon \mathbb{Z}\{\lambda_{0},\lambda_1,\dots,\lambda_{p-1}\} \to \mathbb{Z}\]///,
        "by sending ", TEX///$\lambda_i\mapsto 1$///, ". The transfer ",
        TEX///\[\text{tr}_e^{C_p} \colon \mathbb{Z} \to \mathbb{Z}\{\lambda_{0},\lambda_1,\dots,\lambda_{p-1}\}\]///,
        "by sending ", TEX///$x\mapsto x\cdot \left(\sum^{p-1}_{0} \lambda_i\right)$///, "."
        },
    EXAMPLE lines///                                                                         
    makeComplexRepresentationMackeyFunctor(5)
    ///
}

document{
    Key => {makeRealRepresentationMackeyFunctor, (makeRealRepresentationMackeyFunctor, ZZ)},
    Headline => "constructs the real representation Mackey functor",
    Usage => "makeRealRepresentationMackeyFunctor(p)",
    Inputs => {
	ZZ => "p" => {"a prime number ", TEX///$p$///},
	},
    Outputs => {
	CpMackeyFunctor => {"the real representation Mackey functor for the group ", TEX///$C_p.$///},
	},
    PARA{
        "The ", EM "real representation", " Mackey functor of a group ", TEX///$G$///, " is defined by sending ", TEX///$G/H$///, " to the Grothendieck group of real orthogonal ", TEX///$G$/// , "-representations. The transfer and restriction come from induction and restriction of ", TEX///$G$///,"-representations. When ", TEX///$G$///, " is a cyclic group of prime order p, with p odd, this admits a nice form. The underlying module is given by ", TEX///$\mathbb{Z}$///, " with trivial conjugation action, while the fixed module is ", TEX///$\mathbb{Z}\{\lambda_{0},\lambda_1,\dots,\lambda_{(p-1)/2}\}$///, ". The element ", TEX///$\lambda_i$///, " for ",TEX///$i>0$/// ,
        " represents the two dimensional real representation given by rotation by ",TEX///$(2\pi i)/p$///, " radians. The element ", TEX///$\lambda_0$///," represents the trivial one dimensional representation.  Restriction is defined as",
        TEX///\[\text{res}_e^{C_p} \colon \mathbb{Z}\{\lambda_{0},\lambda_1,\dots,\lambda_{(p-1)/2}\} \to \mathbb{Z}\]///,
        "by sending ", TEX///\[\lambda_i\mapsto \begin{cases} 2 & i>0 \\ 1 & i=0. \end{cases} \]///, "Transfer is of the form ",
        TEX///\[\text{tr}_e^{C_p} \colon \mathbb{Z} \to \mathbb{Z}\{\lambda_{0},\lambda_1,\dots,\lambda_{(p-1)/2}\}\]///,
        "by sending ", TEX///$x\mapsto x\cdot \left(\sum^{(p-1)/2}_{0} \lambda_i\right)$///, "."
        },
    EXAMPLE lines///                                                                         
    makeRealRepresentationMackeyFunctor(5)
    ///
}