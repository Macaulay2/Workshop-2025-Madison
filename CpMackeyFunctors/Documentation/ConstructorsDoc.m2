doc ///
    Key
        "constructing examples of Mackey functors"
    Headline
        info about how to construct explicit Mackey functors
    Description
        Text
            New Mackey functors can be built using @TO "makeCpMackeyFunctor"@, which constructs a Mackey functor out of a restriction, transfer, and conjugation @TO "matrix"@, in that order. For instance:
            
            todo - example here
            
            There are other constructors for various types of Mackey functors, including the @TO2(makeBurnsideMackeyFunctor,"Burnside Mackey functor")@, the @TO2(makeRealRepresentationMackeyFunctor,"real")@ and @TO2(makeComplexRepresentationMackeyFunctor,"complex representation")@ Mackey functors, @TO2(makeFixedPointMackeyFunctor,"fixed point")@ and @TO2(makeOrbitMackeyFunctor,"orbit")@ Mackey functors, {\em free} Mackey functors on an @TO2(makeUnderlyingFreeMackeyFunctor,"underlying element")@, and of course the @TO2(makeZeroMackeyFunctor,"zero Mackey functor")@. 

        Example
            makeZeroMackeyFunctor(5)

        Text
            Furthermore we have {\bf random constructors} which allow us to build a @TO2(makeRandomCpMackeyFunctor,"random Mackey functor")@ over the group $C_p$.

    SeeAlso
        "background on Mackey functors"
        "the abelian category of Mackey functors"
        "explicit applications of the CpMackeyFunctors package"
///

doc ///
    Key 
        makeBurnsideMackeyFunctor
        (makeBurnsideMackeyFunctor, ZZ)
        makeFixedFreeMackeyFunctor
    Headline
        constructs the Burnside Mackey functor
    Usage
        makeBurnsideMackeyFunctor(p)
        makeFixedFreeMackeyFunctor(p)
    Inputs 
        p : ZZ
            a prime number=> {
    Outputs
        : CpMackeyFunctor
            the $C_p$-Mackey functor $\underline{A}$
    Description
        Text
            The {\em Burnside Mackey functor} of a group $G$ is defined by sending $G/H$ to the Burnside ring $A(H)$, with transfer and restriction coming from transfer and restriction of finite $G$-sets. When $G$ is a cyclic group of prime order, these maps admit a particularly nice form. The underlying module is given by $\mathbb{Z}$, with trivial conjugation actin=on, while the fixed module is $\mathbb{Z}\{1,t\}$. Restriction is defined as $$\text{res}_e^{C_p} \colon \mathbb{Z}\{1,t\} \to \mathbb{Z}$$ by sending $t\mapsto p$. Transfer is of the form $$\text{tr}_e^{C_p} \colon \mathbb{Z} \to \mathbb{Z}\{1,t\}$$ by sending $x\mapsto xt$.
        Example
            makeBurnsideMackeyFunctor(5)
///                                                                       

doc ///
    Key
        makeUnderlyingFreeMackeyFunctor
        (makeUnderlyingFreeMackeyFunctor, ZZ)
    Headline
        constructs the free Mackey functor on an underlying generator
    Usage
        makeUnderlyingFreeMackeyFunctor(p)
    Inputs
        p : ZZ
            a prime number $p$
    Outputs
        : CpMackeyFunctor
            the free Mackey functor on an underlying generator
    Description
        Text
            The free $C_p$-Mackey functor on an underlying generator represents the functor $\text{Mack}_{C_p} \to \text{Ab}$ which sends a Mackey functor $M$ to its underlying level $M(C_p/e)$. This means there is a natural isomorphism $$\text{Hom}_{\text{Mack}_{C_p}}\left(\underline{B},M \right) \cong M(C_p/e).$$
        Text
            In components, the underlying module is the free module on the $C_p$-set  $C_p/e=\{1,\gamma,\gamma^2,\ldots,\gamma^{p-1}\}$, with conjugation induced by the $C_p$-action of left multiplication. The fixed module is the module $\ZZ$. Restriction is the map $$\text{res}_e^{C_p} \colon \ZZ \to \ZZ\{1,\gamma,\ldots,\gamma^{p-1}\}$$ by sending $1\mapsto 1+\gamma+\cdots+\gamma^{p-1}$. The transfer is the map $$\text{tr}_e^{C_p} \colon \ZZ\{1,\gamma,\ldots,\gamma^{p-1}\} \to \ZZ$$ by sending $\gamma^i\mapsto 1$ for all $i$.
        Example
            makeUnderlyingFreeMackeyFunctor(5)
///

doc ///
    Key
        makeFixedPointMackeyFunctor
        (makeFixedPointMackeyFunctor, ZZ, Matrix)
    Headline
        constructs the fixed-point Mackey functor of a C_p-module
    Usage
        makeFixedPointMackeyFunctor(p,C)
    Inputs
        p : ZZ
            a prime number $p$
        C : Matrix
            an order $p$ automorphism $c$
    Outputs
        : CpMackeyFunctor
            the fixed-point Mackey functor of the $C_p$-module specified by $\text{c}$.
    Description
        Text
            Given a module $M$ with $C_p$-action specified by an automorphism $\text{c}\colon M\to M$, there is a naturally associated {\em fixed-point Mackey functor} $\text{FP}(M)$. The underlying module is the module $M$ with conjugation given by $\text{c}$. The fixed module is the module $M^{C_p}$ of fixed points of the action, represented as the kernel of $1-\text{c}$. The restriction is the map $$\text{res}_e^{C_p} \colon M^{C_p} \to M$$ induced by the inclusion of fixed-points. The transfer is the map $$\text{tr}_e^{C_p} \colon M\to M^{C_p}$$ given by $1+\text{c}+\text{c}^2+\cdots+\text{c}^{p-1}$.
        Example
            makeFixedPointMackeyFunctor(2,matrix{{0,1,0},{1,0,0},{0,0,1}})
///

doc ///
    Key
        makeOrbitMackeyFunctor
        (makeOrbitMackeyFunctor, ZZ, Matrix)
    Headline
        constructs the orbit Mackey functor of a C_p-module
    Usage
        makeOrbitMackeyFunctor(p,C)
    Inputs
        p : ZZ
            a prime number $p$
        C : Matrix
            an order $p$ automorphism $\text{c}$
    Outputs
        : CpMackeyFunctor
            the orbit Mackey functor of the $C_p$-module specified by $\text{c}$
    Description
        Text
            Given a module $M$ with $C_p$-action specified by an automorphism $\text{c}\colon M\to M$, there is a naturally associated {\em orbit Mackey functor} $\text{O}(M)$. The underlying module is the module $M$ with conjugation given by $\text{c}$. The fixed module is the quotient module $M_{C_p}$ of the action, represented as the cokernel of $1-\text{c}$. The restriction is the map $$\text{res}_e^{C_p} \colon M_{C_p} \to M$$ given by $1+\text{c}+\text{c}^2+\cdots+\text{c}^{p-1}$. The transfer is the map $$\text{tr}_e^{C_p} \colon M\to M_{C_p}$$ induced by the projection to the quotient.
        Example
            makeOrbitMackeyFunctor(2,matrix{{0,1,0},{1,0,0},{0,0,1}})
///

doc ///
    Key
        makeZeroMackeyFunctor
        (makeZeroMackeyFunctor, ZZ)
    Headline
        constructs the zero Mackey functor for the group
    Usage
        makeZeroMackeyFunctor(p)
    Inputs
        p : ZZ
            a prime number $p$
    Outputs
        : CpMackeyFunctor
            the zero Mackey functor for $C_p$
    Description
        Text
            Perhaps the easiest Mackey functor is the ", EM "zero Mackey functor", " which has the zero-module as both the underlying and fixed modules. This is the zero object in the abelian category of Mackey functors, and is an important object to have for homological algebra.
        Example
            makeZeroMackeyFunctor(2)
///

doc ///
    Key
        makeComplexRepresentationMackeyFunctor
        (makeComplexRepresentationMackeyFunctor, ZZ) 
    Headline
        constructs the complex representation Mackey functor
    Usage
        makeComplexRepresentationMackeyFunctor(p)
    Inputs
        p : ZZ
            a prime number $p$
    Outputs
        : CpMackeyFunctor
            the complex representation Mackey functor for the group $C_p$.
    Description
        Text
            The {\em complex representation} Mackey functor of a group $G$ is defined by sending $G/H$ to the Grothendieck group of complex $G$-representations. The transfer and restriction come from induction and restriction of $G$-representations. When $G$ is a cyclic group of prime order, this admits a nice form. The underlying module is given by $\ZZ$ with trivial conjugation action, while the fixed module is $\mathbb{Z}\{\lambda_{0},\lambda_1,\dots,\lambda_{p-1}\}$. The element $\lambda_i$ represents the one dimensional complex representation given by multiplication by $e^{2\pi i/p}$. Restriction is defined as",
            \[\text{res}_e^{C_p} \colon \mathbb{Z}\{\lambda_{0},\lambda_1,\dots,\lambda_{p-1}\} \to \mathbb{Z}\]
            by sending $\lambda_i\mapsto 1$. The transfer
            \[\text{tr}_e^{C_p} \colon \mathbb{Z} \to \mathbb{Z}\{\lambda_{0},\lambda_1,\dots,\lambda_{p-1}\}\]
            sends $x\mapsto x\cdot \left(\sum^{p-1}_{0} \lambda_i\right)$.
        Example
            makeComplexRepresentationMackeyFunctor(5)
    SeeAlso
        makeRealRepresentationMackeyFunctor
///

doc ///
    Key
        makeRealRepresentationMackeyFunctor
        (makeRealRepresentationMackeyFunctor, ZZ)
    Headline
        constructs the real representation Mackey functor
    Usage
        makeRealRepresentationMackeyFunctor(p)
    Inputs
        p : ZZ
            a prime number $p$
    Outputs
        : CpMackeyFunctor
            the real representation Mackey functor for the group $C_p$
    
    Description
        Text
            The {\em real representation} Mackey functor of a group $G$ is defined by sending $G/H$ to the Grothendieck group of real orthogonal $G$-representations. The transfer and restriction come from induction and restriction of $G$-representations. When $G$ is a cyclic group of prime order p, with p odd, this admits a nice form. The underlying module is given by $\ZZ$ with trivial conjugation action, while the fixed module is $\mathbb{Z}\{\lambda_{0},\lambda_1,\dots,\lambda_{(p-1)/2}\}$. The element $\lambda_i$ for $i>0$ represents the two dimensional real representation given by rotation by $(2\pi i)/p$ radians. The element $\lambda_0$ represents the trivial one dimensional representation.  Restriction is defined as
            \[\text{res}_e^{C_p} \colon \mathbb{Z}\{\lambda_{0},\lambda_1,\dots,\lambda_{(p-1)/2}\} \to \mathbb{Z}\]
            by sending
            \[\lambda_i\mapsto \begin{cases} 2 & i>0 \\ 1 & i=0. \end{cases} \]
            Transfer is of the form
            \[\text{tr}_e^{C_p} \colon \mathbb{Z} \to \mathbb{Z}\{\lambda_{0},\lambda_1,\dots,\lambda_{(p-1)/2}\}\]
            by sending $x\mapsto x\cdot \left(\sum^{(p-1)/2}_{0} \lambda_i\right)$.
        Example
            makeRealRepresentationMackeyFunctor(5)
    SeeAlso
        makeComplexRepresentationMackeyFunctor
///
