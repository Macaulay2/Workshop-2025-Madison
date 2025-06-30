--CpMackeyFunctors.m2
newPackage(
    "CpMackeyFunctors",
    Version=>"1.0",
    Date=>"June 24, 2025",
    Authors=>{
        {Name=>"Thomas Brazelton",
            Email=>"brazelton@math.harvard.edu",
            HomePage=>"https://tbrazel.github.io/"},
        {Name=>"Ben Spitz",
            Email=>"benspitz@virginia.edu",
            HomePage=>"https://benspitz.com/"},
        {Name=>"Michael R. Zeng",
            Email=>"zengrf@uw.edu",
            HomePage=>"https://sites.google.com/uw.edu/zengrf/"},
        {Name=>"Chenglu Wang",
            Email=>"cwang@math.harvard.edu",
            HomePage=>"url"}
        {Name=>"David Chan",
            Email=>"chandav2@msu.edu",
            HomePage=>"https://sites.google.com/view/davidchanmath"}
        {Name=>"Add Yourself",
            Email=>"email",
            HomePage=>"url"}
    },
    Headline=>"for working with Mackey functors for cyclic groups of prime order",
    PackageImports=>{},
    PackageExports=>{},
    AuxiliaryFiles => true,
    DebuggingMode=>true
)

load "./CpMackeyFunctors/Code/Helpers.m2"
-- helper commands are internal only, and not exported

-- Code for the CpMackeyFunctor type
load "./CpMackeyFunctors/Code/MackeyFunctor.m2"
export{
    -- CpMackeyFunctors
    "CpMackeyFunctor",
    "isWellDefinedCpMackeyFunctor",
    "makeCpMackeyFunctor",
    "primeorder",
    "Underlying",
    "Fixed",
    "Res",
    "Tr",
    "Conj",
    "getPrimeOrder",
    "getUnderlyingModule",
    "getFixedModule",
    "getRestriction",
    "getTransfer",
    "getConjugation"
}


----------------------------
----------------------------
-- Testing
----------------------------
----------------------------

-- load "./tests/HelpersTest.m2"
-- load "./tests/MackeyFunctorsTest.m2"
