newPackage(
    "EquivariantEhrhart",
    Version => "0.1",
    Date => "June 30, 2025",
    Headline => "compute equivariant Ehrhart series of lattice polytopes",
    Authors => {
        {
            Name => "Oliver Clarke",
            Email => "oliver.clarke@durham.ac.uk",
            HomePage => "https://oliverclarkemath.com"
        },
        {
            Name => "Mike Cummings",
            Email => "mike.cummings@uwaterloo.ca",
            HomePage => "https://mikecummings.ca"
        },
        {
            Name => "Sean Grate",
            Email => "sean.grate@auburn.edu",
            HomePage => "https://seangrate.com"
        }
    },
    Keywords => {"Combinatorics", "Convex Geometry"},
    PackageImports => {"BettiCharacters", "Permutations", "RationalPolytopes"},
    AuxiliaryFiles => true,
    DebuggingMode => false
)

export {
    -- methods
    "conjugacyClasses",
    "cycleTypeRepresentatives",
    "equivariantEhrhartSeries",
    "fixedPolytope",
    "generateGroup",
    "isSymmetric",
    "orbitPolytope",
    "representationRing",

    -- options
    "OnlyListRepresentatives",
    "ReturnHStarList",
    "ReturnPartitionList",
    "ReturnTable"
}

-* Code section *-
load "./EquivariantEhrhart/code.m2"

-* Documentation section *-
--beginDocumentation()
--load "./EquivariantEhrhart/docs.m2"

-* Test section *-
load "./EquivariantEhrhart/tests.m2"

end--

-* Development section *-
restart
debug needsPackage "EquivariantEhrhart"
check "EquivariantEhrhart"

uninstallPackage "EquivariantEhrhart"
restart
installPackage "EquivariantEhrhart"
viewHelp "EquivariantEhrhart"
