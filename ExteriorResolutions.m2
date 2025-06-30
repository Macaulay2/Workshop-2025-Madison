newPackage("ExteriorResolutions",
    Version => "1.1",
    Date => "24 June 2025",
    Headline => "Injective resolutions over exterior algebras",
    Authors => {
	{Name => "Penelope Beall",           Email => "pbeall@ucdavis.edu",	       HomePage => "https://pbeall.github.io" },
    {Name => "Michael K. Brown",         Email => "mkb0096@auburn.edu",	       HomePage => "http://webhome.auburn.edu/~mkb0096/" },
	{Name => "Caitlin M. Davis",	     Email => "cmdavis22@wisc.edu",	       HomePage => "https://sites.google.com/wisc.edu/caitlindavis/home" },
	{Name => "Andrew Karsten",	     Email => "akk0071@auburn.edu",	       HomePage => "https://www.auburn.edu/cosam/departments/math/students/grad/graduate-students.htm" },
	{Name => "Jiaheng Li",		     Email => "henryli@gatech.edu",	       HomePage => "" },
	{Name => "Jianuo Zhou",		     Email => "jzhou632@gatech.edu",	       HomePage => "https://math.gatech.edu/people?field_job_type_tid=14"},
	{Name => "Boyana Martinova",	     Email => "boyana.martinova@gmail.com",    HomePage => "https://sites.google.com/view/bmartinova/home"},
	{Name => "Mahrud Sayrafi",	     Email => "mahrud@umn.edu",		       HomePage => "https://mahrud.github.io/" },
	{Name => "Gregory G. Smith",	     Email => "ggsmith@mast.queensu.ca",       HomePage => "https://mast.queensu.ca/~ggsmith/" },
	{Name => "Sreehari Suresh Babu",     Email => "sreeharisbabu183@gmail.com",    HomePage => "https://sreehari183.github.io/" }
    },
    PackageExports => {"Complexes"},
    Keywords => {"Commutative Algebra"}
  )

export {
    --methods
    "injRes",
	"priddyComplex"
    }

--------------------------------------------------
--- Injective resolutions
--------------------------------------------------

--Input: 
--Output:
injRes = method();
injRes(Complex) := (C) -> (
)

priddyComplex = method(TypicalValue=>Complex)
priddyComplex(Matrix) := (m) -> (
	
	
	
	d0 := m;
	d1 := syz d0;
	d2 := syz d1;
	d3 := syz d2;
	
	--complex {d0, d1, d2, d3}
	
	complex apply(10, i -> priddyDifferential(i, m))
)
priddyComplex(List) := (l) -> (
	priddyComplex matrix {l}
)


priddyDifferential = method(TypicalValue=>Matrix)
priddyDifferential(ZZ, Matrix, Ring) := (i, m, S) -> (
	
    monsSrc = basis(-i, S);
	monsTgt = basis(-i+1, S);
	
	--src =  -- P_{-i}\
	--tgt = 
)


beginDocumentation()

--------------------------------------------------
--- Differential Modules
--------------------------------------------------

doc ///
   Key 
      ExteriorResolutions
   Headline 
      Package that implements homological constructions over the exterior algebra
   Description
      Text
   SeeAlso
   References
       Text
///


--- TESTS

TEST /// 
    (1+1 == 2) == true
///


