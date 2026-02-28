## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Package description

lmmprobe implements a partitioned Empirical Bayes ECM algorithm for sparse
high-dimensional linear mixed modeling, as described in Zgodic et al. (2025)
<doi:10.1007/s11222-025-10649-z>.

## Test environments

* Local: Arch Linux, R 4.5.2, GCC 15.2.1
* GitHub Actions (all passed R CMD check with 0 errors, 0 warnings):
  - ubuntu-latest, R release
  - ubuntu-latest, R devel
  - ubuntu-latest, R oldrel-1
  - macOS-latest (ARM64), R release
  - windows-latest, R release

## Notes

* The only NOTE on CRAN incoming checks is "New submission" and the
  maintainer address, which is expected for a first submission.

* The package includes a subset (500 of 15,378 probes) of a publicly available
  gene expression dataset (GEO accession GSE65391, Banchereau et al. 2016,
  Cell 165(3):551-565). The full dataset is available from the Gene Expression
  Omnibus.

* The package uses OpenMP for parallelization in C++ code. The OpenMP header
  is conditionally included (#ifdef _OPENMP) so the package compiles and runs
  correctly on platforms without OpenMP support (e.g., macOS with Apple Clang),
  falling back to single-threaded execution.
