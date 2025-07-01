newEliminationTemplate = method(Options => {})
newEliminationTemplate (RingElement, Ideal) := o -> (aVar, J) -> (
    (sh, mp) := getTemplate(aVar, basis(R/J), J);
    M := getTemplateMatrix(shifts, monomialPartition, J);
    new EliminationTemplate from {
	shifts => sh,
        monomialPartition => mp,
        templateMatrix => M,
        actionVariable => aVar,
        ideal => J
    }
    )

