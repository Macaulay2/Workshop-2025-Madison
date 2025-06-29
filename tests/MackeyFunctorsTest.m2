-- Note: to load the necessary packages from the repository, either:
-- first call installPackage in testbot.m2:
--  installPackage("Elimination", FileName => "/home/macaulay/Elimination.m2")
-- or specify the path to needsPackage:
--  needsPackage("Elimination", FileName => "/home/macaulay/Elimination.m2")

-- test code and assertions

needsPackage("CpMackeyFunctors")


U:=ZZ^0;
F:=ZZ^0;
R:=matrix({});
C:=R;
T:=R;
M:=makeCpMackeyFunctor({U,F},R,T,C)
assert(isWellDefinedCpMackeyFunctor M)