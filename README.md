# dissol240

This repository contains all sources for the dissolFoam package.
The package has been tested with OpenFOAM v2.4.0

Installation:

1) Apply patches and recompile OpenFOAM v2.x.x
2) Compile applications - libraries, solvers, utilities - with Allwmake
3) Run cases (approximate execution times in docs/SI.pdf)

Notes:
    1) surfRoughGen may require paths to fftw3 library and header files
        edit applications/utilities/surfRoughGen/Make/options
    2) Some of the blockMeshDict files require v2.4.0 (for mesh grading)
    3) dissolFoam requires both the 0 and current directories to restart
    4) More information about the package can be found in docs/SI.pdf

