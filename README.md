# SDPA

The SDPA (SemiDefinite Programming Algorithm) is a software package for solving semidefinite program (SDP). It is based on a Mehrotra-type predictor-corrector infeasible primal-dual interior-point method. The SDPA handles the standard form SDP and its dual. It is implemented in C++ language utilizing the LAPACK(BLAS) and MUMPS https://mumps-solver.org for matrix computation. The SDPA incorporates dynamic memory allocation and deallocation. So, the maximum size of an SDP that can be solved depends on the size of computational memory which user’s computer loads. 

Older versions of SDPA and manuals are available at
https://sdpa.sourceforge.net/

The latest version 7.4.4 supports the following features

1: ILP64 (int = 64bit) support. It can handle huge problems!

2: Parallel computation (multi-threading) for faster search direction computation, etc.

ILP64 support can solve huge problems, but it consumes more memory than LP64 and slows down the speed a little.
First select LP64 (int = 32bit) or ILP (int = 64bit).

See the INSTALL file for information on how to create binaries.

Last updated January 1, 2024
