-- New versions of Normaliz (tested for Normaliz 3.10.2)
-- have the ability to compute with rational polytopes
-- however, we have to make a small change to doWriteNmzData function
-- below is a working example

needsPackage "Normaliz"
debug Normaliz

-- writes the given data in a normaliz input file
doWriteNmzData = method()
-- writes several matrices in a normaliz input file
doWriteNmzData List := matrices -> (
    checkNmzFile("doWriteNmzData");
    outf := nmzFile | ".in" << "";
    for p in matrices do (
	sgr := p#0;
	nmzMode := p#1;
	outf << numRows sgr << endl;
	outf << numColumns sgr << endl;
	--

	if ring sgr =!= ZZ and ring sgr =!= QQ then error("matrix with non-rational entries");

	for i from 0 to numRows sgr - 1 do (
	    s := "";
	    for j from 0 to numColumns sgr - 1
	    do s = s | toString(sgr_(i,j)) | " "; -- this handles ZZ and QQ entries
	    outf << s << endl;
	    );

	-- Until version 3.9.4, input type normal_toric_ideal was called lattice_ideal
	if normalizProgram#"version" < "3.10" and nmzMode == "normal_toric_ideal" then nmzMode = "lattice_ideal";
	outf << nmzMode << endl);
    outf << close)

-- The following function doesn't need to be redefined necessarily;
-- this version does not delete temporary files (may be useful for debugging)
runNormaliz = method(Options => options normaliz)
runNormaliz(Matrix, String) := opts -> (sgr, nmzMode) -> runNormaliz({(sgr, nmzMode)}, opts)
runNormaliz List := opts -> s -> (
    print opts;

    setNmzFile();

    if 0 < #opts.grading then s = append(s, (matrix {opts.grading}, "grading"));
    doWriteNmzData s;

    dir := select(".*/", nmzFile);
    runDir := if dir != {} then dir#0 else null;
    runProgram(normalizProgram, getNmzExec(), collectNmzOptions() | baseFilename nmzFile,
	RunDirectory => runDir, Verbose => debugLevel > 0);

    -- return nothing if .gen is not generated
    if not nmzGen then ( if nmzFilename == "" then rmNmzFiles(); return );

    if not opts.allComputations then (
	nmzData := readNmzData "gen";
	rc := new RationalCone from { "gen" => nmzData, "inv" => getNumInvs() };
	-- if nmzFilename == "" then rmNmzFiles();
	return rc);

    -- read all files written
    files := { "inv" => getNumInvs() };
    suffixes := { "gen","egn","esp","tri","typ","ht1","ext","tgn" };
    for s in suffixes do if fileExists(nmzFile | "." | s) then files = append(files, s => readNmzData s);

    L := readMultipleNmzData "cst";
    files = append(files, "sup" => L#0);
    files = append(files, "equ" => L#1);
    files = append(files, "cgr" => L#2);

    C := new RationalCone from files;

    -- if nmzFilename == "" then rmNmzFiles(); -- removing used normaliz files
    C)


end

restart
load "normalizWithRationalPolytopes.m2"

M = matrix {{1,5}}
(runNormaliz(M,"polytope"))#"gen"

M=transpose matrix "1/2, -1/3, 1/4; 1/2, -1/3, -1/2"
C = (normaliz(M,"polytope"))

needsPackage "RationalPolytopes"
P = convexHull transpose M
H = ehrhartSeries P

denom = C#"inv"#"hilbert series denom"
use ring H
D = product for i in denom list (1-t^i)
ring D === ring H

D*H
C#"inv"#"hilbert series num"



M=transpose matrix "1, -1/3, 1/4; 1/2, -1/3, -2/2"
(normaliz(M,"polytope"))#"gen"

M = matrix "1,0;-1,0;0,1/20;-1,11/20"
elapsedTime (normaliz(M,"polytope"))
