-- -*- coding: utf-8 -*-
-*
`TestPackage` is an abbreviated version of `FirstPackage`.
Create your own package using `packageTemplate`.
*-

newPackage(
    "TestPackage",
    Version => "0.1",
    Date => "April 1, 2025",
    Authors => {{
	    Name => "First Lastname",
	    Email => "lastname@email.edu",
	    HomePage => "https://github.com/Macaulay2/FirstLastname"
	    }},
    Headline => "a Macaulay2 package to be tested"
    )

export {"firstFunction"}

firstFunction = method(TypicalValue => String)
firstFunction ZZ := String => n -> if n == 1 then "Hello World!" else "D'oh!"

TEST ///
    assert ( firstFunction 2 == "D'oh!" )
///
