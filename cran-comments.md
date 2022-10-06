## 2019.03.28 - Initial submission
We have revised the vignette to substantially reduce check time.

Regarding the note, we have intentionally left the connector word (and) in lowercase to accurately reflect the acronym of our package name. The title meets the requirement of title case.

## Test environments
* local OS X install, R 3.5.3
* ubuntu 14.04 (on travis-ci), R 3.5.3
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 1 note

* This is a new release.


## 2019.04.02 - Updated version 0.2.1
Cleaned up unused dependencies and imports.
Fixed expected variable inputs to sqrt() and round() with explicit casts in C++ code to remove Solaris build warning.


## 2020.01.17 - Updated version 1.0.0
Based on internal testing and journal reviewer comments, we have added several updates in a new major release of the package.

## Test environments
* Oracle Solaris 10, x86, 32 bit, R-patched (experimental)
* local Mac OS X install, R 3.6.2
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Debian Linux, R-devel, GCC ASAN/UBSAN
* Fedora Linux, R-devel, clang, gfortran
* Ubuntu 14.04 (on travis-ci), R 3.5.3, 3.6.2
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 note


## 2020.03.22 - Updated version 1.1.0
Bug fixes:
- transpositional error in cellular automata dispersal function
- incorrect modification of transition matrices based on approach to carrying capacity
  in competition density function
Performance improvements:
- carrying capacity raster is generated once per timestep if function is used.
  
## Test environments
* Oracle Solaris 10, x86, 32 bit, R-patched (experimental)
* local Mac OS X install, R 3.6.2
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Debian Linux, R-devel, GCC ASAN/UBSAN
* Fedora Linux, R-devel, clang, gfortran

## R CMD check results
0 errors | 0 warnings | 0 note


## 2021.02.25 - Updated version 1.2.0
Updates:
- Revised transition matrix stochasticity documentation for clarity
- Fixed error in growth transition functions documentation
- Added functionality and documentation for two-sex transition matrices
Bug fixes:
- Fixed density_dependence_stages length mismatch in dispersal function
Performance improvements:
- revised testing code to speed up test runtimes

## Test environments
* Oracle Solaris 10, x86, 32 bit, R-release, Oracle Developer Studio 12.6
* macOS 10.13.6 High Sierra, R-release, brew
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Ubuntu Linux 20.04.1, R-release, GCC
* Debian Linux, R-devel, GCC ASAN/UBSAN
* Fedora Linux, R-devel, clang, gfortran

## R CMD check results
0 errors | 0 warnings | 2 notes
checking sizes of PDF files under 'inst/doc' ... NOTE
  Unable to find GhostScript executable to run checks on size reduction
checking top-level files ... NOTE
Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.


## 2021.03.03 - Minor updated version 1.2.1
Updates:
- Revised deprecated C++ function
- Reduced installation size in vignettes

## R CMD check results
0 errors | 0 warnings | 0 notes


## 2022.10.05 - Updated version 1.3.0
Updates:
- Re-wrote plotting and visualisation functions and updated vignettes accordingly
- Added support for two-sex simulations
- Added minimum distance parameter to cellular automata dispersal function
Bug fixes:
- Fixed neighbouring cell indexing error in kernel-based dispersal function
Misc:
- Removed reference to deprecated function in future package
- Removed empty NEWS file

## Test environments
* macOS 11.5.2 Big Sur, R-release, brew
* Windows Server 2022, R-devel, 64 bit
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Debian Linux, R-devel, GCC ASAN/UBSAN
* Fedora Linux, R-devel, clang, gfortran

## R CMD check results
0 errors | 1 warning | 1 note
  WARNING
  'qpdf' is needed for checks on size reduction of PDFs
  NOTE
  Found the following files/directories:
    'lastMiKTeXException'