doc ///
    Key
        "list of common Mackey functors"
    Headline
        a list of common Cp-Mackey functors
    Description
        Text
            Here we include a table of common Mackey functors for the group $C_p$. The Lewis diagrams have to be formatted in @TT{"\\substack"}@ in order to compile, apologies for their appearance.
            
            @TABLE{
                ((BOLD "  Ravenel symbol  "),(BOLD "  Lewis diagram  "), (BOLD "  Description  "), (BOLD "  Implementation  ")),
                ((TEX "$\\square$"), (TEX "$$\\substack{\\ZZ \\\\ p\\uparrow\\downarrow 1\\\\ \\ZZ \\\\ \\circlearrowleft \\\\ 1}$$"),("fixed point Mackey functor of ", TEX"$\\ZZ$"," as a trivial ", TEX"$C_p$", "-module"), (TT "makeFixedPointMackeyFunctor(p,id_(ZZ^1))")),
                ((TEX "$\\circ$"), (TEX "$$\\substack{\\ZZ/p \\\\ 0\\uparrow\\downarrow 0\\\\ 0 \\\\ \\circlearrowleft \\\\ 0}$$"), ("todo"), ("todo")),
                ((TEX "$\\overset{\\slash}{\\square}$"), (TEX "$$\\substack{\\ZZ \\\\ 1\\uparrow\\downarrow p\\\\ \\ZZ \\\\ \\circlearrowleft \\\\ 1}$$"), ("todo"), ("todo")),
                ((TEX "$\\ominus^n$"), (TEX "$$\\substack{\\ZZ / (q^n-1) \\\\ 1\\uparrow\\downarrow \\sum_{j=0}^{p-1} q^{jn} \\\\ \\ZZ/ (q^{pn} - 1) \\\\ \\circlearrowleft \\\\ 1}$$"), ("todo"), ("todo"))
            }@
///