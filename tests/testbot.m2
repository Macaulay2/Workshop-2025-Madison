-* This is the default testbot.m2 provided for this workshop.

   It contains Macaulay2 code and is automatically run every time you
   either open a pull request
   or push a new commit to a project branch on the GitHub respository,
   unless you add the string "[skip ci]" to the commit message.
*-

-- To add your project directories containing Macaulay2 source code files
-- to the load path, uncomment and edit the following line:    
--   path = { currentDirectory() | "pathToMyFiles/", "absolutePathToMyFiles/", "etc.../" } | path 
-- Terminate each directory name with a "/".

-- Edit the following lines to preload and check your package 
needsPackage ("TestPackage", FileName => "tests/TestPackage.m2")
check TestPackage

-- and/or run a series of examples either with
load "tests/example.m2"
-- or with
capture get "tests/example.m2"
-- at every push on GitHub.
